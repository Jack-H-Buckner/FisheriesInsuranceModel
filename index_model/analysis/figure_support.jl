#####################################################################
#####################################################################
### Shared helpers for the figure scripts.
###
### run_bankruptcy_figures.jl and run_demand_figures.jl both walk the
### same precision ladder, look results up by key, and write tag
### suffixed output.  That machinery lives here so the two cannot
### disagree about what F1 means or where a preliminary run's files go.
###
### Deliberately NOT included: the scenario lists themselves.  Each
### script names the scenarios it runs explicitly, so what a figure is
### built from is readable in the file that builds it.
### Jack H. Buckner, Oregon State University, 2026
### Generated with Claude Code.
#####################################################################
#####################################################################


# Harmonic mean of precision and recall.  On the iso premium arm re = pr and
# this reduces to pr; on the full coverage arm re = 1 and it is 2pr/(pr+1).
# One definition, because the two arms are plotted against a shared F1 axis and
# a discrepancy there would silently shift one curve relative to the other.
f1_score(pr, re) = 2 * pr * re / (pr + re)


# Look a row up by its key.  `nothing` when that scenario was skipped, which the
# figure code filters on so an unsolved scenario drops a point rather than
# taking the whole panel down.
row_for(rows, key) = begin
    i = findfirst(r -> r.key == String(key), rows)
    i === nothing ? nothing : rows[i]
end

rows_for(rows, keys) = filter(!isnothing, [row_for(rows, k) for k in keys])

# Ladder series, ascending precision, ready to plot.
ladder(rows, keys) = sort(rows_for(rows, keys), by = r -> r.pr)


# Output paths carry the results tag, so a preliminary run can never overwrite a
# production figure or table.
tagged_path(outdir, name, tag, ext) = joinpath(outdir, name * tag * ext)


# Scenarios whose solution has not landed yet are reported and skipped rather
# than raising: a partially synced sweep should still produce the panels it can
# support.  Returns true when both the insurance solve and the shared no
# insurance solve for its base are present.
function solved(scenario; tag::AbstractString = "", needs_no_insurance::Bool = false)
    isfile(joinpath(output_dir(scenario; insurance = true, tag = tag), "sol.jld2")) || return false
    needs_no_insurance || return true
    return isfile(joinpath(output_dir(scenario; insurance = false, tag = tag), "sol.jld2"))
end


# One line per scenario: what would run, and whether it can.
function print_scenario_table(scenarios; tag::AbstractString = "",
                              needs_no_insurance::Bool = false, io = stdout)
    w = maximum(length(s.scenario) for s in scenarios) + 2
    println(io, rpad("scenario", w), rpad("arm", 16),
            join(lpad.(["pr", "re", "cv", "F1", "solved"], 9)))
    n = 0
    for s in scenarios
        ok = solved(s.scenario; tag = tag, needs_no_insurance = needs_no_insurance)
        n += ok
        println(io, rpad(s.scenario, w), rpad(s.arm, 16),
                join(lpad.([s.pr, s.re, s.cv, round(f1_score(s.pr, s.re), digits = 3),
                            ok ? "yes" : "NO"], 9)))
    end
    println(io, "\n", n, " of ", length(scenarios), " ready")
    return n
end


# Raster resolution for every PNG.  Plots treats `size` as px at a nominal
# 100 dpi, so size/100 is the figure's PHYSICAL size in inches and `dpi` only
# changes how finely that is rasterised.  600 rather than 300 because these
# figures are mostly thin grey lines and open markers, which is exactly what
# alias badly at 300.
#
# The PDF written alongside is vector and has no resolution at all - prefer it
# for submission, and use the PNG for drafts and screens.
include(joinpath(@__DIR__, "png_dpi.jl"))

const FIGURE_DPI = 600

# Plot styling shared by every panel, so the figures read as one set.
# size = (450, 300) is 4.5 x 3.0 inches, a sensible single panel print size, so
# the 10pt fonts below are 10pt on the printed page rather than shrinking when
# the figure is scaled into a column.
const PANEL = (size = (450, 300), legendfontsize = 10, tickfontsize = 10,
               dpi = FIGURE_DPI)
const PRECISION_LABEL = "Index accuracy (Precision & recall)"

# Every composite panel is drawn against F1, so the utility and bankruptcy
# panels share one accuracy axis and the two arms are comparable on it.  The
# label is a constant for the same reason f1_score is: two panels disagreeing
# about what the x axis means would be invisible in the finished figure.
const F1_LABEL = "Index accuracy (F1 score)"
const REDUCTION_LABEL = "percentage reduction of \n the bankruptcy rate"
