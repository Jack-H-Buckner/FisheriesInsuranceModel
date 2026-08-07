#####################################################################
#####################################################################
### Demand metric 3: probability of bankruptcy.
###
### Ported from analysis/table_4.jl, the production version of the
### base model's calculation (src/calculate_bankruptcy_rate.jl is a
### broken stub: N = 2, T = 2, an unconditional throw(), and a
### data.frame(df) R-ism).  The method is unchanged:
###
###   N paired trajectories of T periods, simulated forward from a
###   stationary biological state at a fixed level of savings, with and
###   without access to index insurance;
###   Poisson MLE hazard  λ̂ = (# bankruptcies) / (total exposure time);
###   a paired nonparametric bootstrap over the N trajectories;
###   reported as the 10 year probability P10 = 1 - exp(-10 λ̂).
###
### Three things are done differently, all forced by the port rather
### than chosen:
###
### 1. Each trajectory builds its OWN MCRandomVariable.  `simulate`
###    resamples the shock nodes in place (X() in the package's
###    analysis.jl), so the single Xmc that table_4.jl shares across
###    Threads.@threads is written by every thread at once.  That is a
###    data race, and it also makes the per iteration Random.seed!
###    ineffective, so the base model's numbers are not in fact
###    reproducible run to run.  A fresh Xmc per trajectory is cheap
###    and fixes both.
###
### 2. `paired = true` (the default) re-seeds before each of the two
###    simulation calls, so the with and without insurance runs see the
###    SAME shock sequence.  Mortality and payouts evolve independently
###    of savings, so the two runs then share an identical biological
###    and index path and differ only through insurance.  table_4.jl
###    calls the design paired but runs the two simulations back to
###    back off one stream, so they actually get different shocks and
###    are paired only through the shared starting state.  Pass
###    `paired = false` to reproduce that behaviour.
###
### 3. Grids and solutions come from `rebuild_problem`, so the value
###    function loaded is the one this grid was solved on.  table_4.jl
###    hardcodes its own sizes, which disagree with src/run.jl's.
###
### Jack H. Buckner, Oregon State University, 2026
### Generated with Claude Code.
#####################################################################
#####################################################################

# Guarded: the driver may already have included setup.jl through another metric,
# and including it twice in one session re-runs the BASE_PARAMETERS loop and
# warns on the const.
@isdefined(rebuild_problem) || include(joinpath(@__DIR__, "setup.jl"))

using Random, Statistics


# Defaults, all from analysis/table_4.jl.  `burnin` is the length of the
# biological trajectory simulated to reach a stationary stock state before the
# financial simulation starts; table_4.jl uses T for it, so it is defaulted to
# T here as well rather than being introduced as a new free parameter.
const BANKRUPTCY_SETTINGS = (
    N          = 60,     # trajectories per scenario
    T          = 250,    # periods per trajectory (the censoring horizon)
    B          = 2000,   # bootstrap resamples
    horizon    = 10,     # years, for the reported P10
    seed       = 123,    # trajectory seeds are hash((seed, i))
    boot_seed  = 456,
    s0_savings = 0.0,    # savings at the start of each trajectory
)


# ------------------------------------------------------------------
# One trajectory
# ------------------------------------------------------------------

# A fresh shock sampler.  `Nquad` sizes the Gauss Hermite quadrature object that
# init_ranom_variables also builds and that is discarded here: the Bellman step
# uses prob.X, built by rebuild_problem at the solved Nquad, while the realised
# shocks come from this MCRandomVariable, whose sampler draws the recruitment
# deviation directly.  So a small Nquad costs nothing and changes nothing.
new_shocks(p) = init_ranom_variables(p; Nquad = 2)[2]


"""
    stationary_start(p, Xmc, burnin; savings = 0.0)

A starting state for the financial simulation: the biological state reached
after `burnin` periods of the population model, with savings held at `savings`
and the firm solvent.

The savings state is `log(savings - s̄)`, so the default `savings = 0.0` gives
`log(-s̄)` — the firm starts with nothing saved and its full borrowing capacity
available, matching `analysis/table_4.jl:109`.
"""
function stationary_start(p, Xmc, burnin::Int; savings = 0.0)
    savings > p.s̄ || throw(ArgumentError(
        "savings = $savings is at or below the borrowing floor s̄ = $(p.s̄); " *
        "the trajectory would start bankrupt"))
    sim, = simulate_trajectory(p, Xmc, burnin)
    return vcat([1.0, log(savings - p.s̄)], sim[:, end])
end


# Time to bankruptcy from one simulated path.  `Ibt` is the bankruptcy indicator
# over the T+1 states, 1.0 solvent and 2.0 bankrupt; index 1 is the initial
# state, so a first hit at index k means k-1 periods elapsed before bankruptcy.
# A path that never goes bankrupt is censored at T.
function time_to_bankruptcy(Ibt, T::Int)
    k = findfirst(>(1.0), Ibt)
    k === nothing && return (exposure = float(T), event = false)
    return (exposure = float(k - 1), event = true)
end


# ------------------------------------------------------------------
# The N paired trajectories
# ------------------------------------------------------------------

"""
    bankruptcy_paths(prob_ins, prob_no; kwargs...)

Simulate `N` paired trajectories with and without index insurance and return the
exposure time and bankruptcy indicator of each.

Returns `(exposure, event, exposure_no, event_no, s0)`, the per trajectory
survival data that `hazard_rates` turns into a rate.

Keywords default to `BANKRUPTCY_SETTINGS`; `paired` controls whether the two
runs share a shock stream (see the header of this file).
"""
function bankruptcy_paths(prob_ins, prob_no;
                          N::Int          = BANKRUPTCY_SETTINGS.N,
                          T::Int          = BANKRUPTCY_SETTINGS.T,
                          burnin::Int     = T,
                          seed            = BANKRUPTCY_SETTINGS.seed,
                          s0_savings      = BANKRUPTCY_SETTINGS.s0_savings,
                          paired::Bool    = true,
                          verbose::Bool   = true)

    assert_comparable(prob_ins, prob_no)
    p = prob_ins.p

    # Vector{Bool}, not the BitVector that falses(N) gives: BitVector packs 64
    # flags into one word, so two threads writing different trajectories can
    # write the same word and lose one of them.  Vector{Bool} is byte addressed
    # and safe.  (analysis/table_4.jl uses falses(N) here.)
    exposure    = zeros(N);  event    = zeros(Bool, N)
    exposure_no = zeros(N);  event_no = zeros(Bool, N)
    s0_all      = zeros(5, N)

    verbose && println("simulating ", N, " paired trajectories of ", T,
                       " periods on ", Threads.nthreads(), " threads",
                       paired ? " (common random numbers)" : "")

    # Per trajectory seeding.  Random.seed! seeds the TASK local default RNG,
    # and Threads.@threads runs each chunk of iterations as one task, so seeding
    # at the top of the body makes trajectory i reproducible independently of
    # the thread count and the scheduling.  Hashing (seed, i) avoids the
    # correlations between adjacent seeds that seed + i would leave.
    Threads.@threads for i in 1:N
        Random.seed!(hash((seed, i)) % UInt32)

        # One sampler per trajectory: simulate mutates it in place, so a shared
        # one would be a data race across threads.
        Xmc = new_shocks(p)
        s0 = stationary_start(p, Xmc, burnin; savings = s0_savings)
        s0_all[:, i] .= s0

        # Drawn from the seeded stream, so it is a deterministic function of i.
        sim_seed = rand(UInt32)

        Random.seed!(sim_seed)
        x = simulation(prob_ins, El, Xmc, T; with_policy = false, s0 = s0)

        # Common random numbers: re-seeding puts the second run on the same
        # shock stream as the first.  Xmc can be reused because X() overwrites
        # every sampled row of its nodes, so no state carries over.
        paired && Random.seed!(sim_seed)
        x_no = simulation(prob_no, El, Xmc, T; with_policy = false, s0 = s0)

        a = time_to_bankruptcy(x.Ibt, T)
        b = time_to_bankruptcy(x_no.Ibt, T)
        exposure[i], event[i]       = a.exposure, a.event
        exposure_no[i], event_no[i] = b.exposure, b.event
    end

    return (exposure = exposure, event = event,
            exposure_no = exposure_no, event_no = event_no, s0 = s0_all)
end


# ------------------------------------------------------------------
# Hazard rate, bootstrap, and the reported probability
# ------------------------------------------------------------------

# Poisson MLE hazard λ̂ = (# events) / (total exposure time).  Under a constant
# hazard / exponential time to bankruptcy model this is both the MLE and an
# unbiased rate estimator.  Zero exposure cannot arise from bankruptcy_paths
# (every trajectory is solvent at t = 0 and so contributes at least one period
# once it fails at t >= 1), but a bootstrap resample of an all zero column
# would, so it is guarded rather than left to produce a silent NaN.
hazard(event, exposure) = sum(exposure) > 0 ? sum(event) / sum(exposure) : NaN

# Rate ratio, with the degenerate case made explicit: if the no insurance arm
# records no bankruptcies at all the ratio is undefined, not zero or infinite.
rate_ratio(λ, λ_no) = (isnan(λ) || isnan(λ_no) || λ_no == 0) ? NaN : λ / λ_no

# 10 year bankruptcy probability under the constant hazard model.  Monotone in
# λ, so a CI on λ transforms endpoint by endpoint into a valid CI on P10.  P10
# is reported as the interpretable marginal probability; ratios stay on the
# instantaneous λ scale, which is horizon invariant where P10 ratios are not.
prob_within(λ, horizon) = isnan(λ) ? NaN : 1 - exp(-horizon * λ)


"""
    hazard_rates(paths; B, horizon, boot_seed)

Turn the per trajectory survival data from `bankruptcy_paths` into hazard rates,
10 year bankruptcy probabilities and bootstrap confidence intervals.

The bootstrap resamples TRAJECTORIES, taking the with and without insurance run
of the same trajectory together.  Pairing preserves the correlation the common
random numbers induce between the two arms, so the CI on the ratio is not
inflated by variation that is common to both.
"""
function hazard_rates(paths;
                      B::Int     = BANKRUPTCY_SETTINGS.B,
                      horizon    = BANKRUPTCY_SETTINGS.horizon,
                      boot_seed  = BANKRUPTCY_SETTINGS.boot_seed)

    event, exposure       = paths.event,    paths.exposure
    event_no, exposure_no = paths.event_no, paths.exposure_no
    N = length(event)

    λ    = hazard(event, exposure)
    λ_no = hazard(event_no, exposure_no)

    Random.seed!(boot_seed)
    λ_boot        = zeros(B)
    λ_no_boot     = zeros(B)
    logratio_boot = zeros(B)
    reduction_boot = zeros(B)
    for b in 1:B
        idx = rand(1:N, N)
        λb    = hazard(event[idx],    exposure[idx])
        λb_no = hazard(event_no[idx], exposure_no[idx])
        λ_boot[b]        = λb
        λ_no_boot[b]     = λb_no
        # 0 events, 0 no insurance events, or either NaN all give a non finite
        # log ratio; those resamples are dropped below rather than quantiled.
        logratio_boot[b] = log(rate_ratio(λb, λb_no))
        # The risk DIFFERENCE, which stays estimable where the ratio does not.
        # When insurance removes bankruptcy entirely the ratio is 0 with an
        # undefined interval - every resample has a zero numerator - and yet the
        # reduction is the largest it can be.  That is the common case here, so
        # the difference is reported alongside the base model's ratio rather
        # than instead of it.
        reduction_boot[b] = prob_within(λb_no, horizon) - prob_within(λb, horizon)
    end

    # 95% percentile intervals.  With few events per scenario a resample can
    # contain none at all, so the λ draws are filtered too, not just the ratio.
    qs(v) = (u = filter(isfinite, v); isempty(u) ? (NaN, NaN) : Tuple(quantile(u, [0.025, 0.975])))
    λ_lo,    λ_hi    = qs(λ_boot)
    λ_no_lo, λ_no_hi = qs(λ_no_boot)
    lr_lo,   lr_hi   = qs(logratio_boot)
    red_lo,  red_hi  = qs(reduction_boot)
    ratio_lo, ratio_hi = isnan(lr_lo) ? (NaN, NaN) : (exp(lr_lo), exp(lr_hi))

    return (
        N              = N,
        n_events       = sum(event),
        n_events_no_I  = sum(event_no),
        exposure       = sum(exposure),
        exposure_no_I  = sum(exposure_no),
        lambda         = λ,
        lambda_lo      = λ_lo,
        lambda_hi      = λ_hi,
        lambda_no_I    = λ_no,
        lambda_no_I_lo = λ_no_lo,
        lambda_no_I_hi = λ_no_hi,
        lambda_ratio   = rate_ratio(λ, λ_no),
        ratio_lo       = ratio_lo,
        ratio_hi       = ratio_hi,
        horizon        = horizon,
        P10            = prob_within(λ,    horizon),
        P10_lo         = prob_within(λ_lo, horizon),
        P10_hi         = prob_within(λ_hi, horizon),
        P10_no_I       = prob_within(λ_no,    horizon),
        P10_no_I_lo    = prob_within(λ_no_lo, horizon),
        P10_no_I_hi    = prob_within(λ_no_hi, horizon),
        # Reduction in the 10 year bankruptcy probability attributable to
        # insurance.  Always estimable, unlike the ratio, and on the same
        # interpretable scale as P10 itself.  Its interval comes from the paired
        # resample, so the correlation the common random numbers induce between
        # the two arms is carried through rather than being thrown away.
        P10_reduction    = prob_within(λ_no, horizon) - prob_within(λ, horizon),
        P10_reduction_lo = red_lo,
        P10_reduction_hi = red_hi,
        # Share of trajectories that fail inside the simulated window.  Unlike
        # P10 this makes no constant hazard assumption, so a large gap between
        # the two is a sign that assumption is doing real work.
        frac_bankrupt      = mean(event),
        frac_bankrupt_no_I = mean(event_no),
        boot_dropped       = B - count(isfinite, logratio_boot),
        B                  = B,
    )
end


"""
    bankruptcy_rate(prob_ins, prob_no; kwargs...)

Probability of bankruptcy with and without access to index insurance.

Simulates the paired trajectories and returns the hazard rates, the 10 year
bankruptcy probabilities and their bootstrap intervals as a single NamedTuple.
Keywords are those of `bankruptcy_paths` and `hazard_rates`.
"""
function bankruptcy_rate(prob_ins, prob_no;
                         N::Int        = BANKRUPTCY_SETTINGS.N,
                         T::Int        = BANKRUPTCY_SETTINGS.T,
                         burnin::Int   = T,
                         B::Int        = BANKRUPTCY_SETTINGS.B,
                         horizon       = BANKRUPTCY_SETTINGS.horizon,
                         seed          = BANKRUPTCY_SETTINGS.seed,
                         boot_seed     = BANKRUPTCY_SETTINGS.boot_seed,
                         s0_savings    = BANKRUPTCY_SETTINGS.s0_savings,
                         paired::Bool  = true,
                         verbose::Bool = true)

    t0 = time()
    paths = bankruptcy_paths(prob_ins, prob_no; N = N, T = T, burnin = burnin,
                             seed = seed, s0_savings = s0_savings,
                             paired = paired, verbose = verbose)
    res = hazard_rates(paths; B = B, horizon = horizon, boot_seed = boot_seed)
    res = merge(res, (T = T, seed = seed, paired = paired,
                      s0_savings = s0_savings,
                      minutes = (time() - t0) / 60))
    verbose && print_bankruptcy(res)
    return res
end


"""
    bankruptcy_rate(scenario; kwargs...)

Convenience method: rebuild `scenario` and its base's shared no insurance
solution from `results_index/`, then run the calculation.
"""
function bankruptcy_rate(scenario::AbstractString; kwargs...)
    r  = rebuild_problem(scenario)
    r0 = rebuild_problem(scenario; insurance = false)
    return bankruptcy_rate(r.prob, r0.prob; kwargs...)
end


# ------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------

fmt(x; digits = 4) = isnan(x) ? "NA" : string(round(x, digits = digits))
ci(x, lo, hi; digits = 4) = "$(fmt(x; digits = digits)) ($(fmt(lo; digits = digits)), $(fmt(hi; digits = digits)))"

function print_bankruptcy(res; io = stdout)
    println(io, "events (insurance / none): ", res.n_events, " / ",
            res.n_events_no_I, " out of ", res.N,
            "   exposure ", Int(res.exposure), " / ", Int(res.exposure_no_I), " periods")
    println(io, "P", res.horizon, "        = ", ci(res.P10, res.P10_lo, res.P10_hi))
    println(io, "P", res.horizon, "_no_I   = ", ci(res.P10_no_I, res.P10_no_I_lo, res.P10_no_I_hi))
    println(io, "reduction  = ", ci(res.P10_reduction, res.P10_reduction_lo, res.P10_reduction_hi))
    println(io, "λ ratio    = ", ci(res.lambda_ratio, res.ratio_lo, res.ratio_hi))
    res.boot_dropped > 0 && println(io, "note: ", res.boot_dropped, " of ", res.B,
        " bootstrap resamples had no events in one arm, so the ratio CI is based on ",
        res.B - res.boot_dropped, " draws; read the reduction instead")
    return nothing
end
