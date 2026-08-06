#####################################################################
#####################################################################
### Demand metric 2: change in value from access to insurance.
###
### The companion to coverage.jl.  Where coverage_by_budget evaluates
### the POLICY function at each simulated state, this evaluates the
### VALUE function of the insurance and no insurance solutions at the
### same state and differences them:
###
###   value_diff_by_budget   ΔV and the implied change in the
###                          certainty equivalent, integrated over the
###                          stationary distribution of biological
###                          states, as a function of the budget
###                          available above the borrowing floor
###
### Both solutions are built on the same state grid (build_grids leaves
### the state grid untouched when insurance = false), so the two value
### functions can be evaluated at identical state vectors.  Only the
### action grid differs, which is what the comparison is about.
### Jack H. Buckner, Oregon State University, 2026
### Generated with claude code. 
#####################################################################
#####################################################################

# Guarded: coverage.jl includes setup.jl too, and including it twice in one
# session re-runs the BASE_PARAMETERS loop and warns on the const.
@isdefined(rebuild_problem) || include(joinpath(@__DIR__, "setup.jl"))


# ------------------------------------------------------------------
# Value to consumption units
# ------------------------------------------------------------------
#
# V is a discounted sum of per period utility, so V*(1-δ) is the constant per
# period utility that would deliver the same value, and inverting the utility
# function turns that into consumption units.
#
# CAUTION, and the reason two mappings are reported below.  `utility`
# (src/income_model.jl:44) is  u(κ) = (κ^(1-γ)-1)/(1-γ) + 1,  but
# `inverse_utility` (:52) is  (U(1-γ)+1)^(1/(1-γ)),  which inverts the CRRA
# WITHOUT the +1.  The two disagree by that constant:
#
#   γ != 0 (RA sets, γ = 1.001): the offset enters the CE multiplicatively and
#           cancels in a ratio, so pct_increase_CE is unaffected.
#   γ == 0 (no_RA sets):         u(κ) = κ, so ce = E[κ] + 1.  The +1 inflates
#           the denominator and DAMPS the reported percentage change.
#
# `pct_increase_CE` reproduces src/analysis.jl:291-299 exactly, so the index
# numbers stay comparable with the base model's published table.
# `pct_increase_consumption` uses the true inverse of `utility` and is the
# number to quote if the no_RA percentages are read as consumption changes.
# For γ != 0 the two agree to rounding; at γ = 0 they do not, by construction.

# Base model's mapping, src/analysis.jl:201.
certainty_equivalent(V, δ, γ) = inverse_utility(V * (1 - δ), γ)

# Inverse of `utility` as actually defined: solve (κ^(1-γ)-1)/(1-γ) + 1 = u.
consumption_equivalent(V, δ, γ) = ((V * (1 - δ) - 1) * (1 - γ) + 1)^(1 / (1 - γ))


# ------------------------------------------------------------------
# ΔV over the stationary distribution of the stock
# ------------------------------------------------------------------

# Same log spaced span as coverage.jl's default_budgets, repeated here so this
# file can be included on its own without pulling in the plotting code.
budget_grid(prob; n::Int = 15) =
    exp.(range(extrema(prob.V.states[2, prob.V.states[1, :] .== 1])..., length = n))


# The two solutions have to be the same problem apart from the premium grid,
# otherwise the difference is meaningless.  Cheap to check, and it catches the
# two mistakes that actually happen: passing two insurance solves, and pairing a
# scenario with another base's no insurance solve.
function assert_comparable(prob_ins, prob_no)
    probe = [1.0, 0.0, 0.0, 0.0, 1.0]
    all(prob_no.u(probe, prob_no.p)[2, :] .== 0.0) || throw(ArgumentError(
        "prob_no is not a no insurance solution: its premium grid is not [0.0]"))
    any(prob_ins.u(probe, prob_ins.p)[2, :] .> 0.0) || throw(ArgumentError(
        "prob_ins has a degenerate premium grid; it looks like a no insurance solve"))
    prob_ins.V.Bslines[1].grid == prob_no.V.Bslines[1].grid || throw(ArgumentError(
        "the two solutions are on different state grids, so V cannot be " *
        "differenced state by state"))
    prob_ins.δ ≈ prob_no.δ || throw(ArgumentError(
        "different discount rates: $(prob_ins.δ) vs $(prob_no.δ)"))
    prob_ins.p.γ ≈ prob_no.p.γ || throw(ArgumentError(
        "different risk aversion: $(prob_ins.p.γ) vs $(prob_no.p.γ)"))
    return nothing
end


"""
    value_diff_by_budget(prob_ins, prob_no, sim_states, budgets; nmax = 2000)

Change in value from access to index insurance, as a function of the budget
available above the borrowing floor, integrated over the stationary
distribution of biological states.

For each budget `s - s̄`, the savings state is held at `log(s - s̄)` while
`(log b, Δ log b, M)` ranges over `sim_states`, and both value functions are
evaluated at that state.  Unlike `coverage_by_budget` this only interpolates
`V`, it does not re-solve the one period problem, so it is cheap and `nmax`
defaults to the whole stationary sample.

Returns one NamedTuple per budget, with

  mean_value_diff          E[V_insurance - V_no_insurance], in utils.  The
                           +1 shift in `utility` is common to both and cancels,
                           so this is the offset free measure.
  pct_increase_CE          100*(E[ce_I] - E[ce_N])/E[ce_N], reproducing
                           src/analysis.jl:291-299 for comparability with the
                           base model's table
  pct_increase_consumption the same ratio using the true inverse of `utility`;
                           see the note at the top of this file, and prefer this
                           one when reading no_RA (γ = 0) results as consumption
  frac_better              fraction of states where insurance is weakly better,
                           which must be 1.0 up to interpolation error: access
                           to a larger action set cannot lower V

Pair the arguments with `rebuild_problem(scenario)` and
`rebuild_problem(scenario; insurance = false)`; the latter resolves to the
base's shared no insurance directory.
"""
function value_diff_by_budget(prob_ins, prob_no, sim_states,
                              budgets = budget_grid(prob_ins);
                              nmax::Int = 2000)
    assert_comparable(prob_ins, prob_no)

    states = thin_states(sim_states, nmax)
    n = size(states, 2)
    p = prob_ins.p
    δ, γ = prob_ins.δ, p.γ

    out = NamedTuple[]
    for s in budgets
        st = log(s)                      # the savings state is log(s - s̄)

        Vi = Vector{Float64}(undef, n)
        Vn = Vector{Float64}(undef, n)
        for j in 1:n
            state = [1.0, st, states[1, j], states[2, j], states[3, j]]
            Vi[j] = prob_ins.V(state)
            Vn[j] = prob_no.V(state)
        end

        ΔV = Vi .- Vn
        ce_i, ce_n = mean(certainty_equivalent.(Vi, δ, γ)), mean(certainty_equivalent.(Vn, δ, γ))
        c_i,  c_n  = mean(consumption_equivalent.(Vi, δ, γ)), mean(consumption_equivalent.(Vn, δ, γ))

        push!(out, (
            budget                   = s,
            savings                  = s + p.s̄,   # savings on the model's own scale
            mean_value_ins           = mean(Vi),
            mean_value_no            = mean(Vn),
            mean_value_diff          = mean(ΔV),
            q10_value_diff           = quantile(ΔV, 0.1),
            q90_value_diff           = quantile(ΔV, 0.9),
            mean_ce_ins              = ce_i,
            mean_ce_no               = ce_n,
            pct_increase_CE          = 100 * (ce_i - ce_n) / ce_n,
            pct_increase_consumption = 100 * (c_i - c_n) / c_n,
            frac_better              = mean(ΔV .>= -1e-8),
            n_states                 = n,
        ))
    end
    return out
end


"""
    value_diff_by_budget(scenario; kwargs...)

Convenience method: rebuild both solutions for `scenario` and its base's no
insurance counterpart, simulate the stationary states, and run the sweep.
"""
function value_diff_by_budget(scenario::AbstractString; n_sim::Int = 2000,
                              seed::Int = 123, budgets = nothing, kwargs...)
    r  = rebuild_problem(scenario)
    r0 = rebuild_problem(scenario; insurance = false)
    sim = stationary_states(r.p, r.Xmc; n = n_sim, seed = seed)
    b = something(budgets, budget_grid(r.prob))
    return value_diff_by_budget(r.prob, r0.prob, sim, b; kwargs...)
end


"""
    plot_value_diff_by_budget(rows; field = :pct_increase_CE, label = "")

Change in welfare against budget.  `field` selects the measure:
`:pct_increase_CE`, `:pct_increase_consumption` or `:mean_value_diff` (the last
is shown with a 10-90% band over the stationary distribution).
"""
function plot_value_diff_by_budget(rows; field::Symbol = :pct_increase_CE,
                                   label = "", plt = nothing, kwargs...)
    b = [r.budget for r in rows]
    y = [getfield(r, field) for r in rows]
    f = plt === nothing ? Plots.plot : Plots.plot!
    ylab = field === :mean_value_diff ? "Δ value (utils)" : "Increase in welfare (%)"

    if field === :mean_value_diff
        lo = [r.q10_value_diff for r in rows]
        hi = [r.q90_value_diff for r in rows]
        return f(b, y; ribbon = (y .- lo, hi .- y), label = label,
                 xlabel = "Budget (savings above the floor)", ylabel = ylab,
                 xscale = :log10, linewidth = 2, fillalpha = 0.15, kwargs...)
    end
    return f(b, y; label = label, xlabel = "Budget (savings above the floor)",
             ylabel = ylab, xscale = :log10, linewidth = 2, kwargs...)
end
