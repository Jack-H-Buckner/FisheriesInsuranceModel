#####################################################################
#####################################################################
### Bankruptcy rates across the precision ladder, and the figures
### built from them.  A script version of
### index_model/bankruptcy_analysis.ipynb.
###
### Runs every scenario the notebook ran, collects ALL of the returned
### statistics into one table, writes it as a CSV, and reproduces the
### four saved figures from that table.
###
###     julia --project=. --threads=10 \
###           index_model/analysis/run_bankruptcy_figures.jl
###
### Options:
###   --N=100        trajectories per scenario (the notebook used 100)
###   --tag=_cheaper read a preliminary results tree instead of
###                  production, and suffix every output file to match
###   --outdir=DIR   where the CSV and PNGs go (default index_model/figures)
###   --quick        N=8, T=60, B=200; plumbing check, not results
###   --dry-run      print the scenario table and exit without solving
###   --figures-only redraw the panels from the CSV already on disk,
###                  simulating nothing.  This is the whole point of
###                  writing the table first: re-cutting a figure (a
###                  style change, a new dpi) must not cost hours
###
### The table is the deliverable; every figure is derived from it and
### nothing is recomputed per panel.  So a figure can be re-cut from
### the CSV without re-running six hours of simulation.
###
### THREE DIFFERENCES FROM THE NOTEBOOK, all deliberate:
###
###  1. `no_RA_int_pr100_re100_cv33` was run twice there - N=100 in
###     cell 3 and N=25 in cell 10 - and the N=25 copy is what cells 12
###     and 13 use as the shared perfect index endpoint AND as the
###     denominator for both curves.  Here it is run once, at the same
###     N as everything else.
###  2. The cell 11 table mislabelled the full coverage arm as
###     `pr75_re75_100_cv33` etc. (they are re = 100, not re = 75), and
###     its last two rows repeated pr60 and pr50 under the names
###     pr60_re42 / pr50_re33, so the pr42 and pr33 results were never
###     written out at all.  Arm, pr and re are carried as their own
###     columns here rather than encoded in a hand written string.
###  3. F1 is computed once as 2*pr*re/(pr+re) instead of the two
###     separately derived expressions in cell 13.  Same numbers - on
###     the iso premium arm 2p^2/2p reduces to p - one definition.
###
### Jack H. Buckner, Oregon State University, 2026
### Generated with Claude Code.
#####################################################################
#####################################################################

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

ENV["GKSwstype"] = "100"   # headless GR, so this runs over ssh / in a job
using Plots, LaTeXStrings

include(joinpath(@__DIR__, "bankruptcy.jl"))
include(joinpath(@__DIR__, "figure_support.jl"))


#####################################################################
### The scenarios, exactly the set the notebook ran
#####################################################################
#
# `key` is what the figure code looks rows up by, so a figure never depends on
# the order of this list.  `arm` follows generate_scenarios.jl: iso_premium is
# re = pr (cost held fixed, accuracy falls), full_coverage is re = 1 (every
# event covered, precision bought with false positives).

# Derived from scenario_table() rather than hand written, so a rung added to the
# sweep (or the iso arm's pr42 / pr33 solves landing) is picked up here without
# this list being edited.  Filtered to the risk neutral base: the bankruptcy
# channel is what drives demand when gamma = 0.
cost_suffix(cost_tag) = isempty(cost_tag) ? "" : "_" * cost_tag
series_key(arm, pr, cost_tag) =
    (arm == "iso_premium" ? "iso_" : "full_") * tag100(pr) * cost_suffix(cost_tag)

function build_ladder(base_tag)
    seen = Set{String}()
    out = NamedTuple[]
    for r in scenario_table()
        r.base_tag == base_tag || continue
        r.name in seen && continue      # the perfect index appears under both arms
        push!(seen, r.name)
        push!(out, (key = series_key(r.arm, r.pr, r.cost_tag), scenario = r.name,
                    arm = r.arm, pr = r.pr, re = r.re,
                    cv = something(r.cv, BASE_PARAMETERS[r.base].cv)))
    end
    return out
end

iso_keys(cost_tag)  = [series_key("iso_premium", pr, cost_tag) for pr in PR_LADDER]
full_keys(cost_tag) = vcat(series_key("iso_premium", 1.0, cost_tag),
                           [series_key("full_coverage", pr, cost_tag)
                            for pr in PR_LADDER if pr < 1.0])

# f1_score, row_for, rows_for, ladder, tagged_path, solved and the panel
# styling constants come from figure_support.jl, shared with run_demand_figures.jl.


#####################################################################
### Arguments
#####################################################################

const FIGURE_DIR = joinpath(normpath(joinpath(@__DIR__, "..")), "figures")

function parse_args(args)
    opts = Dict{Symbol,Any}(:N => 100, :T => BANKRUPTCY_SETTINGS.T,
                            :B => BANKRUPTCY_SETTINGS.B, :tag => "",
                            :base => "no_RA_int", :outdir => FIGURE_DIR,
                            :quick => false, :dryrun => false, :figonly => false)
    for a in args
        if a == "--quick";       opts[:quick] = true
        elseif a == "--dry-run"; opts[:dryrun] = true
        elseif a == "--figures-only"; opts[:figonly] = true
        elseif startswith(a, "--base=");   opts[:base] = String(a[8:end])
        elseif startswith(a, "--N=");      opts[:N] = parse(Int, a[5:end])
        elseif startswith(a, "--T=");      opts[:T] = parse(Int, a[5:end])
        elseif startswith(a, "--B=");      opts[:B] = parse(Int, a[5:end])
        elseif startswith(a, "--tag=");    opts[:tag] = String(a[7:end])
        elseif startswith(a, "--outdir="); opts[:outdir] = String(a[10:end])
        else
            throw(ArgumentError("unknown option $a; see the header of " *
                                "index_model/analysis/run_bankruptcy_figures.jl"))
        end
    end
    if opts[:quick]
        any(startswith(a, "--N=") for a in args) || (opts[:N] = 8)
        any(startswith(a, "--T=") for a in args) || (opts[:T] = 60)
        any(startswith(a, "--B=") for a in args) || (opts[:B] = 200)
    end
    return opts
end

# Output names carry the tag, so a preliminary run can never overwrite the
# production figures.  The "probabiltiy" misspelling is preserved deliberately:
# the files are referenced by name elsewhere, and renaming them here would
# silently break whatever includes them.
out_path(opts, name, ext) = tagged_path(opts[:outdir], name, opts[:tag], ext)


#####################################################################
### Run the scenarios
#####################################################################

"""
    run_scenarios(opts) -> Vector{NamedTuple}

One row per scenario: its design, every statistic `bankruptcy_rate` returns,
and the derived percentage reduction.  This vector IS the table - pass it to
`DataFrame(rows)` if DataFrames is installed, or to `write_csv`.

A scenario whose solution is missing is reported and skipped rather than
raising, so a partial sweep still produces a table and the figures it supports.
"""
function run_scenarios(scenarios, opts)
    rows = NamedTuple[]
    skipped = String[]

    for (i, s) in enumerate(scenarios)
        # bankruptcy_rate contrasts the scenario against the shared no-insurance
        # solve, so both have to be present.
        if !solved(s.scenario; tag = opts[:tag], needs_no_insurance = true)
            @warn "skipping $(s.scenario): the scenario or its no-insurance solve is missing"
            push!(skipped, s.scenario)
            continue
        end

        println("\n[", i, "/", length(scenarios), "] ", s.scenario,
                "   arm=", s.arm, "  pr=", s.pr, "  re=", s.re, "  cv=", s.cv)
        res = bankruptcy_rate(s.scenario; tag = opts[:tag],
                              N = opts[:N], T = opts[:T], B = opts[:B])

        push!(rows, merge((key = s.key, scenario = s.scenario, arm = s.arm,
                           pr = s.pr, re = s.re, cv = s.cv,
                           f1 = f1_score(s.pr, s.re),
                           results_tag = isempty(opts[:tag]) ? "production" : opts[:tag]),
                          res,
                          # The figure quantity, derived once here so every panel
                          # and the CSV cannot disagree about it.
                          (pct_reduction = 100 * (1 - res.P10 / res.P10_no_I),)))

        # Written after every scenario: this is hours of simulation and a
        # failure on the last one must not throw away the rest.
        write_csv(out_path(opts, "bankruptcy_rates", ".csv"), rows)
    end
    return rows, skipped
end


# The no insurance rate is the denominator of every percentage reduction.  It
# comes from the SAME no insurance solution for every scenario on this base, and
# under common draws every scenario simulates the identical path through it, so
# these must all agree.  Averaging and checking the spread turns that into a
# free diagnostic rather than an assumption - the notebook picked one scenario's
# value arbitrarily, which hid any disagreement.
function baseline_no_insurance(rows)
    vals = [r.P10_no_I for r in rows]
    lo, hi = extrema(vals)
    hi - lo > 1e-8 && @warn """
        P10_no_I differs across scenarios, but every scenario shares one no
        insurance solution and, under common draws, one realised path - so
        these should be identical.  The percentage reductions below use the
        mean.  A spread here means the runs were not made with common draws.""" lo hi
    return sum(vals) / length(vals)
end


#####################################################################
### Figures
#####################################################################


# Bar chart of the 10 year bankruptcy probability down the precision ladder,
# with the no insurance rate as the right hand reference bar.
function figure_levels(rows, keys, no_I, path)
    rs = rows_for(rows, keys)
    isempty(rs) && return nothing
    labels = vcat([string(round(Int, 100 * r.pr)) for r in rs], "No ins.")
    rates  = vcat([r.P10 for r in rs], no_I)

    plt = Plots.bar(labels, rates, color = "lightgrey",
                    xlab = PRECISION_LABEL,
                    ylab = "10 years bankruptcy rate", label = "")
    Plots.plot!(; PANEL...)
    savefig_dpi(plt, path; dpi = FIGURE_DPI)
    println("  wrote ", relpath(path))
    return plt
end


# Percentage reduction against precision, one series per cost structure.
function figure_reduction_by_precision(rows, path)
    fair = ladder(rows, iso_keys(""))
    cv33 = ladder(rows, iso_keys("cv33"))
    (isempty(fair) || isempty(cv33)) && return nothing

    # Categorical x axis, ascending precision, as in the notebook.
    xf = [string(round(Int, 100 * r.pr)) for r in fair]
    xc = [string(round(Int, 100 * r.pr)) for r in cv33]

    plt = Plots.plot(xf, [r.pct_reduction for r in fair], color = "grey", label = "")
    Plots.scatter!(xf, [r.pct_reduction for r in fair], color = "black",
                   label = "loading factor 0%")
    Plots.plot!(xc, [r.pct_reduction for r in cv33], color = "grey", label = "")
    Plots.scatter!(xc, [r.pct_reduction for r in cv33], color = "white",
                   label = "loading factor 50%")
    Plots.plot!(xlab = PRECISION_LABEL, ylab = REDUCTION_LABEL, ylim = (0, 100); PANEL...)
    savefig_dpi(plt, path; dpi = FIGURE_DPI)
    println("  wrote ", relpath(path))
    return plt
end


# Percentage reduction against F1, one series per arm.  Plotting against F1 is
# what makes the two arms comparable: they run down different precision ladders
# but F1 puts both on one accuracy axis.
function figure_reduction_by_f1(rows, path)
    full = sort(rows_for(rows, full_keys("cv33")), by = r -> r.f1)
    iso  = sort(rows_for(rows, iso_keys("cv33")),  by = r -> r.f1)
    (isempty(full) || isempty(iso)) && return nothing

    plt = Plots.plot([r.f1 for r in full], [r.pct_reduction for r in full],
                     color = "grey", label = "")
    Plots.scatter!([r.f1 for r in full], [r.pct_reduction for r in full],
                   color = "grey", label = "recall = 100%")
    Plots.plot!([r.f1 for r in iso], [r.pct_reduction for r in iso],
                color = "grey", label = "")
    Plots.scatter!([r.f1 for r in iso], [r.pct_reduction for r in iso],
                   color = "white", label = "recall = precision")
    Plots.plot!(xlab = F1_LABEL, ylab = REDUCTION_LABEL,
                ylim = (0, 100); PANEL...)
    savefig_dpi(plt, path; dpi = FIGURE_DPI)
    println("  wrote ", relpath(path))
    return plt
end


"""
    make_figures(rows, opts)

Rebuild every saved figure from an existing table.  Separate from
`run_scenarios` on purpose: re-cutting a panel should not mean re-running the
simulation, so this can be called on rows read back from the CSV.
"""
function make_figures(rows, opts)
    isempty(rows) && return nothing
    mkpath(opts[:outdir])
    no_I = baseline_no_insurance(rows)
    println("\nfigures (no insurance baseline P10 = ", round(no_I, digits = 4), ")")

    figure_levels(rows, iso_keys(""), no_I,
                  out_path(opts, "bankruptcy_probabiltiy", ".png"))
    figure_levels(rows, iso_keys("cv33"), no_I,
                  out_path(opts, "bankruptcy_probabiltiy_cv33", ".png"))
    figure_reduction_by_precision(rows, out_path(opts, "change_in_bankruptcy_rate", ".png"))
    figure_reduction_by_f1(rows, out_path(opts, "change_in_bankruptcy_rate_f1", ".png"))
    return nothing
end


#####################################################################
### Main
#####################################################################

function main(args)
    opts = parse_args(args)
    scenarios = build_ladder(opts[:base])

    println("=" ^ 70)
    println("bankruptcy rates and figures, from bankruptcy_analysis.ipynb")
    println("base       ", opts[:base])
    println("scenarios  ", length(scenarios))
    println("settings   N = ", opts[:N], "  T = ", opts[:T], "  B = ", opts[:B])
    println("results    ", isempty(opts[:tag]) ? "results_index/<scenario>/  (production)" :
                           "results_index/<scenario>$(opts[:tag])/  (preliminary)")
    println("output     ", opts[:outdir])
    println("threads    ", Threads.nthreads())
    opts[:quick] && println("QUICK MODE - these numbers are a plumbing check, not a result")
    println("=" ^ 70)

    if opts[:dryrun]
        print_scenario_table(scenarios; tag = opts[:tag], needs_no_insurance = true)
        return 0
    end

    mkpath(opts[:outdir])

    # Redraw only: the table on disk IS the result, so a style or resolution
    # change should never mean re-running the simulation that produced it.
    if opts[:figonly]
        csv = out_path(opts, "bankruptcy_rates", ".csv")
        isfile(csv) || (println("no table at ", relpath(csv),
                                "; run without --figures-only first"); return 1)
        rows = read_csv(csv)
        println("redrawing from ", relpath(csv), " (", length(rows), " rows)")
        make_figures(rows, opts)
        return 0
    end

    rows, skipped = run_scenarios(scenarios, opts)

    println("\n", "=" ^ 70)
    if isempty(rows)
        println("no scenarios could be run; solve them first, or check --tag")
        return 1
    end

    csv = out_path(opts, "bankruptcy_rates", ".csv")
    println("wrote ", length(rows), " rows to ", relpath(csv))
    isempty(skipped) || println("skipped (unsolved): ", join(skipped, ", "))

    make_figures(rows, opts)

    # The table, in the log, so a run is readable without opening the CSV.
    println()
    w = maximum(length(r.scenario) for r in rows) + 2
    println(rpad("scenario", w), rpad("arm", 16),
            join(lpad.(["pr", "re", "cv", "F1", "P10", "P10_no_I", "reduction"], 11)))
    for r in rows
        println(rpad(r.scenario, w), rpad(r.arm, 16),
                join(lpad.([r.pr, r.re, r.cv, round(r.f1, digits = 3),
                            fmt(r.P10), fmt(r.P10_no_I),
                            fmt(r.pct_reduction; digits = 1) * "%"], 11)))
    end
    return 0
end


if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
