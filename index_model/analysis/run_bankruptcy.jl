#####################################################################
#####################################################################
### Driver: bankruptcy rates across the index insurance sweep.
###
### Runs demand metric 3 for every solved scenario and writes one tidy
### row per (scenario, arm) to manuscript_index/index_bankruptcy.csv.
###
### Run from the repository root, AFTER `make sync_index`:
###
###     julia --project=. --threads=8 \
###           index_model/analysis/run_bankruptcy.jl
###
###     ... one scenario only
###     julia --project=. --threads=8 \
###           index_model/analysis/run_bankruptcy.jl no_RA_int_pr60_re100
###
###     ... a fast smoke run, to check the plumbing before committing
###     ... to the full calculation
###     julia --project=. --threads=8 \
###           index_model/analysis/run_bankruptcy.jl --quick
###
### Options (all optional, in any order, before or after the scenario
### names):
###
###   --all            include the RA_* scenarios too; by default only
###                    the no_RA_* sets are run, as in the plan
###   --N=60           trajectories per scenario
###   --T=250          periods per trajectory
###   --B=2000         bootstrap resamples
###   --seed=123       trajectory seeds are hash((seed, i))
###   --savings=0.0    savings at the start of each trajectory
###   --unpaired       give the two arms independent shocks, which is
###                    what analysis/table_4.jl does
###   --quick          N=8, T=60, B=200; for testing, not for results
###   --out=PATH       output CSV
###
### Threads matter: the N trajectories are the parallel unit, and each
### re-solves the one period problem at every step, so this is the
### expensive part of the index analysis.  Give it as many threads as
### the machine has.
### Jack H. Buckner, Oregon State University, 2026
### Generated with Claude Code.
#####################################################################
#####################################################################

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

ENV["GKSwstype"] = "100"   # headless GR: setup.jl pulls in the plotting code
using Plots, LaTeXStrings

include(joinpath(@__DIR__, "bankruptcy.jl"))


#####################################################################
### Arguments
#####################################################################

const DEFAULT_OUT = joinpath(normpath(joinpath(@__DIR__, "..", "..")),
                             "manuscript_index", "index_bankruptcy.csv")

const QUICK = (N = 8, T = 60, B = 200)

function parse_args(args)
    opts = Dict{Symbol,Any}(
        :N => BANKRUPTCY_SETTINGS.N, :T => BANKRUPTCY_SETTINGS.T,
        :B => BANKRUPTCY_SETTINGS.B, :seed => BANKRUPTCY_SETTINGS.seed,
        :savings => BANKRUPTCY_SETTINGS.s0_savings,
        :paired => true, :all => false, :out => DEFAULT_OUT, :quick => false)
    scenarios = String[]

    for a in args
        if !startswith(a, "--")
            push!(scenarios, String(a))
        elseif a == "--all";      opts[:all] = true
        elseif a == "--unpaired"; opts[:paired] = false
        elseif a == "--quick";    opts[:quick] = true
        elseif startswith(a, "--N=");       opts[:N] = parse(Int, a[5:end])
        elseif startswith(a, "--T=");       opts[:T] = parse(Int, a[5:end])
        elseif startswith(a, "--B=");       opts[:B] = parse(Int, a[5:end])
        elseif startswith(a, "--seed=");    opts[:seed] = parse(Int, a[8:end])
        elseif startswith(a, "--savings="); opts[:savings] = parse(Float64, a[11:end])
        elseif startswith(a, "--out=");     opts[:out] = String(a[7:end])
        else
            throw(ArgumentError("unknown option $a; see the header of " *
                                "index_model/analysis/run_bankruptcy.jl"))
        end
    end

    # --quick overrides the sizes but not an explicit --N / --T / --B, so
    # `--quick --T=120` does what it looks like.
    if opts[:quick]
        for (k, v) in pairs(QUICK)
            any(startswith(a, "--$k=") for a in args) || (opts[k] = v)
        end
    end

    # Fail on a typo now rather than after the first scenario has run.
    foreach(scenario_row, scenarios)
    return (scenarios = scenarios, opts = opts)
end


"""
    sweep_scenarios(; include_RA = false)

Scenario names to run, in sweep order.  Defaults to the `no_RA_*` sets: the
plan's demand metrics are computed on the risk neutral base, where demand comes
from the borrowing constraint and bankruptcy risk rather than from risk
aversion.  `--all` adds the `RA_*` sets.
"""
function sweep_scenarios(; include_RA::Bool = false)
    rows = unique_scenarios()
    include_RA || (rows = [r for r in rows if startswith(r.base_tag, "no_RA")])
    return [r.name for r in rows]
end


# A scenario can only be run once both its own solution and its base's shared no
# insurance solution are on disk.  Missing solutions are reported and skipped
# rather than raising, so a partially synced sweep still produces a table.
function solution_status(scenario)
    ins = joinpath(output_dir(scenario; insurance = true),  "sol.jld2")
    non = joinpath(output_dir(scenario; insurance = false), "sol.jld2")
    isfile(ins) || return (ok = false, why = "no insurance solve at $(relpath(ins))")
    isfile(non) || return (ok = false, why = "no no-insurance solve at $(relpath(non))")
    return (ok = true, why = "")
end


#####################################################################
### One row per (scenario, arm)
#####################################################################

# Contract design and derived index properties, so both arms can be plotted
# against the shared precision ladder straight from the CSV without re-deriving
# anything.  p_low and p_stay come from the solved problem's own parameters, not
# from a constant here.
function design_columns(scenario, r, p)
    mc = marginal_chain(p.T)
    m = index_moments(r.pr, r.re, mc.p_low, mc.p_stay)
    return (
        scenario      = scenario,
        base          = r.base,
        cost_tag      = isempty(r.cost_tag) ? "fair" : r.cost_tag,
        pr            = r.pr,
        re            = r.re,
        p_low         = mc.p_low,
        p_stay        = mc.p_stay,
        p_index       = m.p_index,
        return_period = m.return_period,
        fpr           = m.fpr,
        b0            = m.b0,
        cf            = p.cf,
        cv            = p.cv,
        eta_bar       = price_per_exposure(m.p_index, p.cf, p.cv),
    )
end

const RESULT_COLUMNS = (
    :n_events, :n_events_no_I, :N, :exposure, :exposure_no_I,
    :frac_bankrupt, :frac_bankrupt_no_I,
    :P10, :P10_lo, :P10_hi, :P10_no_I, :P10_no_I_lo, :P10_no_I_hi,
    :P10_reduction, :P10_reduction_lo, :P10_reduction_hi,
    :lambda, :lambda_lo, :lambda_hi,
    :lambda_no_I, :lambda_no_I_lo, :lambda_no_I_hi,
    :lambda_ratio, :ratio_lo, :ratio_hi,
    :horizon, :T, :seed, :paired, :s0_savings, :B, :boot_dropped, :minutes,
)

function scenario_rows(scenario, opts)
    r = scenario_row(scenario)
    ri  = rebuild_problem(scenario)
    ri0 = rebuild_problem(scenario; insurance = false)

    res = bankruptcy_rate(ri.prob, ri0.prob;
                          N = opts[:N], T = opts[:T], B = opts[:B],
                          seed = opts[:seed], s0_savings = opts[:savings],
                          paired = opts[:paired])

    design = design_columns(scenario, r, ri.p)
    stats  = NamedTuple{RESULT_COLUMNS}(Tuple(getproperty(res, k) for k in RESULT_COLUMNS))

    # The perfect index (pr = re = 1) is solved once but belongs to both arms, so
    # it is emitted under each label; that gives every curve its endpoint without
    # a second solve.  Every other scenario has exactly one arm.
    return [merge((arm = arm,), design, stats) for arm in arms_for(scenario)]
end


#####################################################################
### Main
#####################################################################

function main(args)
    parsed = parse_args(args)
    opts = parsed.opts
    scenarios = isempty(parsed.scenarios) ?
        sweep_scenarios(include_RA = opts[:all]) : parsed.scenarios

    println("=" ^ 70)
    println("bankruptcy rates for the index insurance sweep")
    println("scenarios  ", length(scenarios))
    println("settings   N = ", opts[:N], "  T = ", opts[:T], "  B = ", opts[:B],
            "  seed = ", opts[:seed], "  savings = ", opts[:savings],
            opts[:paired] ? "" : "  (unpaired)")
    println("threads    ", Threads.nthreads())
    println("output     ", opts[:out])
    opts[:quick] && println("QUICK MODE - these numbers are a plumbing check, not a result")
    println("=" ^ 70)

    rows = NamedTuple[]
    skipped = String[]

    for (i, scenario) in enumerate(scenarios)
        st = solution_status(scenario)
        if !st.ok
            @warn "skipping $scenario: $(st.why)"
            push!(skipped, scenario)
            continue
        end

        println("\n[", i, "/", length(scenarios), "] ", scenario)
        append!(rows, scenario_rows(scenario, opts))

        # Written after every scenario, not once at the end: the full sweep is
        # hours of simulation and a failure on the last scenario must not throw
        # away the ones already done.
        write_csv(opts[:out], rows)
    end

    println("\n", "=" ^ 70)
    if isempty(rows)
        println("no scenarios could be run.  Solutions live in results_index/; ",
                "run `make sync_index` first.")
        return 1
    end

    println("wrote ", length(rows), " rows for ",
            length(unique(r.scenario for r in rows)), " scenarios to ", opts[:out])
    isempty(skipped) || println("skipped (unsolved): ", join(skipped, ", "))

    # A compact version of the table in the log, so a run is readable without
    # opening the CSV.
    println()
    w = maximum(length.([r.scenario for r in rows])) + 2
    println(rpad("scenario", w), rpad("arm", 16),
            join(lpad.(["pr", "re", "P10", "P10_no_I", "reduction", "ratio", "events"], 11)))
    for r in rows
        println(rpad(r.scenario, w), rpad(r.arm, 16),
                join(lpad.([fmt(r.pr; digits = 2), fmt(r.re; digits = 2),
                            fmt(r.P10), fmt(r.P10_no_I), fmt(r.P10_reduction),
                            fmt(r.lambda_ratio),
                            "$(r.n_events)/$(r.n_events_no_I)"], 11)))
    end
    return 0
end


if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
