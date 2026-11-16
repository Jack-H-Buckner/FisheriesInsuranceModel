#####################################################################
#####################################################################
### Change in utility from access to index insurance, across the
### precision ladder.  The welfare half of index_model/analysis.ipynb,
### as a script.
###
###     julia --project=. --threads=8 \
###           index_model/analysis/run_utility_figures.jl
###
### Options:
###   --base=RA_int    scenario prefix.  The notebook used the RISK
###                    AVERSE sets, and it has to: the welfare gain is a
###                    certainty equivalent, and no_RA has γ = 0, so
###                    there is no risk premium to recover there
###   --nsim=2000      stationary draws the value functions average over
###   --burnin=500     draws discarded before the average
###   --tag=_cheaper   read a preliminary results tree and suffix every
###                    output file to match
###   --outdir=DIR     default index_model/figures
###   --dry-run        print the scenario table and exit
###   --figures-only   redraw the figure from the CSV already on disk,
###                    recomputing nothing
###
### Everything is plotted against F1 = 2*pr*re/(pr+re), the harmonic
### mean of precision and recall, so the two arms share one accuracy
### axis and this figure sits on the same axis as the bankruptcy
### panels.  On the iso premium arm F1 reduces to pr.
###
### Coverage (metric 1) was in the notebook and in this script's first
### version; it lives in analysis/coverage.jl and is still tested, but
### no driver here calls it - this file is the welfare metric only, and
### welfare only interpolates V, so it is cheap.
###
### FOUR DEFECTS IN THE NOTEBOOK THIS CORRECTS:
###
###  1. Cells 10-13 read `RA_int_pr*_cv50`, which is no longer in the
###     sweep - the cost variant was renamed cv33 (cv = 0.3333, a 50%
###     loading factor).  Those directories still exist under
###     results_index/ but `scenario_row` now throws for them, so those
###     cells cannot run at all.  The ladder here comes from
###     `scenario_table()`, so it follows the sweep instead of going
###     stale behind it.
###  2. The variables holding those results are named `value_diff_cv50`,
###     `pr_cv50`, `value_diff_100_cv50` while actually holding cv33
###     numbers - a rename waiting to be misread.
###  3. Cell 9 computes `value_diff_100` / `value_diff_50` on the no_RA
###     base and never uses them; cells 11 and 13 then reuse the same
###     names for RA results.
###  4. The stationary draws are taken once here, explicitly, from the
###     base's own parameters.  The notebook rebuilt them mid-way from
###     whichever `problem` happened to be left in scope by the previous
###     loop.
###
### Jack H. Buckner, Oregon State University, 2026
### Generated with Claude Code.
#####################################################################
#####################################################################

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

ENV["GKSwstype"] = "100"
using Plots, LaTeXStrings

include(joinpath(@__DIR__, "value_diffs.jl"))     # pulls in setup.jl
include(joinpath(@__DIR__, "figure_support.jl"))


#####################################################################
### The scenario ladder
#####################################################################
#
# PR_LADDER, COST_VARIANTS, scenario_table and tag100 all come from
# parameters/generate_scenarios.jl, which setup.jl already loads.  Reusing them
# rather than retyping the ladder is what keeps this script in step with the
# sweep: change a rung or a cost variant there and this figure follows.

cost_suffix(cost_tag) = isempty(cost_tag) ? "" : "_" * cost_tag
series_key(arm, pr, cost_tag) =
    (arm == "iso_premium" ? "iso_" : "full_") * tag100(pr) * cost_suffix(cost_tag)

# key is `<arm>_<pr>[_cv33]`, e.g. iso_75_cv33 / full_50.  The perfect index
# (pr = re = 1) is ONE scenario belonging to both arms: it is stored under the
# iso key and both arms' series reference it, so it is never evaluated twice.
function build_ladder(base_tag)
    seen = Set{String}()
    out = NamedTuple[]
    for r in scenario_table()
        r.base_tag == base_tag || continue
        r.name in seen && continue          # the perfect index appears under both arms
        push!(seen, r.name)
        push!(out, (key = series_key(r.arm, r.pr, r.cost_tag), scenario = r.name,
                    arm = r.arm, pr = r.pr, re = r.re,
                    # `nothing` means the scenario inherits the base set's cost
                    cv = something(r.cv, BASE_PARAMETERS[r.base].cv)))
    end
    return out
end

iso_keys(cost_tag)  = [series_key("iso_premium", pr, cost_tag) for pr in PR_LADDER]
full_keys(cost_tag) = vcat(series_key("iso_premium", 1.0, cost_tag),   # shared endpoint
                           [series_key("full_coverage", pr, cost_tag)
                            for pr in PR_LADDER if pr < 1.0])


#####################################################################
### Arguments
#####################################################################

const FIGURE_DIR = joinpath(normpath(joinpath(@__DIR__, "..")), "figures")

function parse_args(args)
    opts = Dict{Symbol,Any}(:base => "RA_int", :nsim => 2000, :burnin => 500,
                            :tag => "", :outdir => FIGURE_DIR, :dryrun => false,
                            :figonly => false)
    for a in args
        if a == "--dry-run";                opts[:dryrun] = true
        elseif a == "--figures-only";       opts[:figonly] = true
        elseif startswith(a, "--base=");    opts[:base] = String(a[8:end])
        elseif startswith(a, "--nsim=");    opts[:nsim] = parse(Int, a[8:end])
        elseif startswith(a, "--burnin=");  opts[:burnin] = parse(Int, a[10:end])
        elseif startswith(a, "--tag=");     opts[:tag] = String(a[7:end])
        elseif startswith(a, "--outdir=");  opts[:outdir] = String(a[10:end])
        else
            throw(ArgumentError("unknown option $a; see the header of " *
                                "index_model/analysis/run_utility_figures.jl"))
        end
    end
    return opts
end

out(opts, name, ext) = tagged_path(opts[:outdir], name, opts[:tag], ext)


#####################################################################
### Run the ladder
#####################################################################

"""
    run_ladder(scenarios, opts) -> (rows, skipped)

One row per scenario: its design, its F1 score, and every statistic
`value_diff_by_budget` returns at a budget of 1.0 above the borrowing floor,
integrated over the stationary distribution of biological states.

The stationary states are drawn ONCE and reused for every scenario, which is
what makes the numbers comparable down the ladder.

Welfare needs the no insurance solve as well as the scenario's own, so both must
be present; a scenario missing either is reported and skipped rather than
raising, so a partial sweep still produces the panel it can support.
"""
function run_ladder(scenarios, opts)
    rows = NamedTuple[]
    skipped = String[]

    first_ok = findfirst(s -> solved(s.scenario; tag = opts[:tag],
                                     needs_no_insurance = true), scenarios)
    first_ok === nothing && return rows, [s.scenario for s in scenarios]

    seed_problem = rebuild_problem(scenarios[first_ok].scenario; tag = opts[:tag])
    sim_states = stationary_states(seed_problem.p, seed_problem.Xmc;
                                   n = opts[:nsim], burnin = opts[:burnin])
    println("stationary states: ", size(sim_states, 2), " draws (burnin ",
            opts[:burnin], ") from ", scenarios[first_ok].scenario)

    for (i, s) in enumerate(scenarios)
        if !solved(s.scenario; tag = opts[:tag], needs_no_insurance = true)
            @warn "skipping $(s.scenario): the scenario or its no-insurance solve is missing"
            push!(skipped, s.scenario)
            continue
        end

        println("\n[", i, "/", length(scenarios), "] ", s.scenario,
                "   arm=", s.arm, "  pr=", s.pr, "  re=", s.re, "  cv=", s.cv,
                "  F1=", round(f1_score(s.pr, s.re), digits = 3))

        r  = rebuild_problem(s.scenario; tag = opts[:tag])
        r0 = rebuild_problem(s.scenario; insurance = false, tag = opts[:tag])
        v  = value_diff_by_budget(r.prob, r0.prob, sim_states, [1.0])[1]

        println("  ΔCE ", round(v.pct_increase_CE, digits = 4), "%",
                "   Δconsumption ", round(v.pct_increase_consumption, digits = 4), "%",
                "   better in ", round(100 * v.frac_better, digits = 1), "% of states")

        push!(rows, (key = s.key, scenario = s.scenario, arm = s.arm,
                     pr = s.pr, re = s.re, cv = s.cv, f1 = f1_score(s.pr, s.re),
                     results_tag = isempty(opts[:tag]) ? "production" : opts[:tag],
                     budget = 1.0, n_states = v.n_states,
                     pct_increase_CE = v.pct_increase_CE,
                     pct_increase_consumption = v.pct_increase_consumption,
                     mean_value_ins = v.mean_value_ins,
                     mean_value_no = v.mean_value_no,
                     mean_value_diff = v.mean_value_diff,
                     q10_value_diff = v.q10_value_diff,
                     q90_value_diff = v.q90_value_diff,
                     mean_ce_ins = v.mean_ce_ins, mean_ce_no = v.mean_ce_no,
                     frac_better = v.frac_better))

        write_csv(out(opts, "utility_by_f1", ".csv"), rows)
    end
    return rows, skipped
end


#####################################################################
### Figure
#####################################################################

const UTILITY_LABEL = "Increase in utility with insurance \n (certainty equivalents)"

# One arm, both cost structures, against F1.  This is panel A of the composite,
# and the only figure analysis.ipynb saved, so it keeps its filename.
#
# Marker convention, carried over from the hand made composite:
#   filled grey = 0% loading,  open white = 50% loading.
function figure_change_in_utility(rows, path)
    fair = ladder(rows, iso_keys(""))
    cv33 = ladder(rows, iso_keys("cv33"))
    (isempty(fair) && isempty(cv33)) && return nothing

    plt = nothing
    for (series, label, colour) in ((fair, "loading factor 0%",  "grey"),
                                    (cv33, "loading factor 50%", "white"))
        isempty(series) && continue
        x = [r.f1 for r in series]
        y = [r.pct_increase_CE for r in series]
        plt = plt === nothing ? Plots.plot(x, y, color = "grey", label = "") :
                                Plots.plot!(x, y, color = "grey", label = "")
        Plots.scatter!(x, y, color = colour, label = label)
    end
    Plots.plot!(xlab = F1_LABEL, ylab = UTILITY_LABEL; PANEL...)
    savefig_dpi(plt, path; dpi = FIGURE_DPI)
    println("  wrote ", relpath(path))
    return plt
end


#####################################################################
### Main
#####################################################################

function main(args)
    opts = parse_args(args)
    scenarios = build_ladder(opts[:base])

    println("=" ^ 70)
    println("change in utility across the precision ladder (was analysis.ipynb)")
    println("base       ", opts[:base])
    println("scenarios  ", length(scenarios))
    println("settings   nsim = ", opts[:nsim], "  burnin = ", opts[:burnin])
    println("results    ", isempty(opts[:tag]) ? "results_index/<scenario>/  (production)" :
                           "results_index/<scenario>$(opts[:tag])/  (preliminary)")
    println("output     ", opts[:outdir])
    println("threads    ", Threads.nthreads())
    println("=" ^ 70)

    if opts[:dryrun]
        print_scenario_table(scenarios; tag = opts[:tag], needs_no_insurance = true)
        return 0
    end

    mkpath(opts[:outdir])

    if opts[:figonly]
        csv = out(opts, "utility_by_f1", ".csv")
        isfile(csv) || (println("no table at ", relpath(csv),
                                "; run without --figures-only first"); return 1)
        rows = read_csv(csv)
        println("redrawing from ", relpath(csv), " (", length(rows), " rows)")
        figure_change_in_utility(rows, out(opts, "change_in_utility", ".png"))
        return 0
    end

    rows, skipped = run_ladder(scenarios, opts)

    println("\n", "=" ^ 70)
    if isempty(rows)
        println("no scenarios could be run; solve them first, or check --base / --tag")
        return 1
    end

    csv = out(opts, "utility_by_f1", ".csv")
    println("wrote ", length(rows), " rows to ", relpath(csv))
    isempty(skipped) || println("skipped (unsolved): ", join(skipped, ", "))

    println("\nfigures")
    figure_change_in_utility(rows, out(opts, "change_in_utility", ".png"))

    println()
    w = maximum(length(r.scenario) for r in rows) + 2
    println(rpad("scenario", w), rpad("arm", 16),
            join(lpad.(["pr", "re", "cv", "F1", "ΔCE %", "Δcons %"], 11)))
    for r in rows
        println(rpad(r.scenario, w), rpad(r.arm, 16),
                join(lpad.([r.pr, r.re, r.cv, round(r.f1, digits = 3),
                            fmt(r.pct_increase_CE; digits = 3),
                            fmt(r.pct_increase_consumption; digits = 3)], 11)))
    end
    return 0
end


if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
