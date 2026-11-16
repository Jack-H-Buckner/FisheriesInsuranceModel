#####################################################################
#####################################################################
### The three panel composite figure.
###
### Replaces index_insurnace_composite.pptx, which was stitched
### together by hand in PowerPoint from two separately generated
### panels and so could not be regenerated when the numbers changed.
###
###     julia --project=. index_model/analysis/run_composite_figure.jl
###
###   A   increase in utility          RA base     iso arm, 0% vs 50% loading
###   B   reduction in bankruptcy      no_RA base  iso arm, 0% vs 50% loading
###   C   reduction in bankruptcy      no_RA base  50% loading, iso vs full coverage
###
### B and C read as a pair: B varies the LOADING FACTOR at a fixed arm,
### C varies the ARM at a fixed loading factor.
###
### All three panels share one x axis, F1 = 2*pr*re/(pr+re).  On the
### iso premium arm that reduces to pr; on full coverage it is
### 2pr/(pr+1).  A useful property: full coverage at pr = 1/3 gives
### F1 = 0.5 exactly, the same as iso at pr = 0.5, so panel C's two
### curves span an identical [0.5, 1.0] range instead of being compared
### over different domains.
###
### READ THE PANELS CAREFULLY.  A is the RISK AVERSE base and B and C
### are RISK NEUTRAL.  They are different agents: a certainty
### equivalent gain requires risk aversion, and no_RA has gamma = 0, so
### there is no risk premium to recover there; the bankruptcy channel
### works through the borrowing constraint and needs no risk aversion
### at all.  Coherent as a figure - two distinct channels by which
### basis risk bites - but NOT one agent's three outcomes.
###
### Reads the two CSVs; it recomputes nothing.  The bankruptcy table is
### hours of simulation, so re-cutting a panel has to be a seconds long
### operation, not a re-run.
###
### Options:
###   --utility=PATH      default index_model/figures/utility_by_f1.csv
###   --bankruptcy=PATH   default index_model/figures/bankruptcy_rates.csv
###   --out=PATH          default index_model/figures/index_insurance_composite.png
###   --layout=1x3|3x1    panel arrangement (default 1x3)
###   --tag=_cheaper      read the preliminary tables and suffix the output
###   --dpi=600           raster resolution of the PNG.  The PDF is
###                       vector and unaffected
###   --no-pdf            skip the PDF written alongside the PNG
###
### Jack H. Buckner, Oregon State University, 2026
### Generated with Claude Code.
#####################################################################
#####################################################################

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

ENV["GKSwstype"] = "100"
using Plots, LaTeXStrings

# Standalone: csv_io.jl and figure_support.jl load no model and no solver, so
# this script starts in about a second.  figure_support.jl calls output_dir in
# `solved`, which is not defined here - that function is simply not used.
include(joinpath(@__DIR__, "csv_io.jl"))
include(joinpath(@__DIR__, "png_dpi.jl"))

f1_score(pr, re) = 2 * pr * re / (pr + re)

# Raster resolution.  Plots treats `size` as px at a nominal 100 dpi, so
# size/100 is the figure's PHYSICAL size in inches and `dpi` only changes how
# finely it is rasterised.  600 rather than 300 because this figure is thin grey
# lines and open markers, which alias badly at 300.
#
# The panel sizes below are therefore chosen as PRINT inches, not screen pixels.
# Getting that wrong is the other half of "blurry": a figure laid out at 15
# inches and then scaled into a 7.5 inch column renders its 9pt labels at 4.5pt.
const DEFAULT_DPI = 600

const PANEL = (legendfontsize = 9, tickfontsize = 9, guidefontsize = 10)
const F1_LABEL = "Index accuracy (F1 score)"
const UTILITY_LABEL = "Increase in utility with insurance \n (certainty equivalents)"
const REDUCTION_LABEL = "percentage reduction of \n the bankruptcy rate"

# Marker convention, carried over from the hand made composite:
#   filled = 0% loading,  open white = 50% loading.
# Panel C is entirely 50% loading, so there the fill distinguishes the ARM
# instead - filled grey = full coverage, open white = iso premium - matching
# change_in_bankruptcy_rate_f1.png.
const FILLED = "grey"
const OPEN   = "white"


#####################################################################
### Arguments
#####################################################################

const FIGURE_DIR = joinpath(normpath(joinpath(@__DIR__, "..")), "figures")

function parse_args(args)
    opts = Dict{Symbol,Any}(:utility => nothing, :bankruptcy => nothing,
                            :out => nothing, :layout => "1x3", :tag => "",
                            :pdf => true, :dpi => DEFAULT_DPI)
    for a in args
        if a == "--no-pdf";                    opts[:pdf] = false
        elseif startswith(a, "--utility=");    opts[:utility] = String(a[11:end])
        elseif startswith(a, "--bankruptcy="); opts[:bankruptcy] = String(a[14:end])
        elseif startswith(a, "--out=");        opts[:out] = String(a[7:end])
        elseif startswith(a, "--layout=");     opts[:layout] = String(a[10:end])
        elseif startswith(a, "--tag=");        opts[:tag] = String(a[7:end])
        elseif startswith(a, "--dpi=");        opts[:dpi] = parse(Int, a[7:end])
        else
            throw(ArgumentError("unknown option $a; see the header of " *
                                "index_model/analysis/run_composite_figure.jl"))
        end
    end
    opts[:layout] in ("1x3", "3x1") || throw(ArgumentError(
        "--layout must be 1x3 or 3x1, got $(opts[:layout])"))

    t = opts[:tag]
    opts[:utility]    = something(opts[:utility],
                                  joinpath(FIGURE_DIR, "utility_by_f1$t.csv"))
    opts[:bankruptcy] = something(opts[:bankruptcy],
                                  joinpath(FIGURE_DIR, "bankruptcy_rates$t.csv"))
    opts[:out]        = something(opts[:out],
                                  joinpath(FIGURE_DIR, "index_insurance_composite$t.png"))
    return opts
end


#####################################################################
### Series selection
#####################################################################

# Rows are selected by their design columns rather than by the `key` string, so
# the composite does not depend on how either producing script happens to name
# its rows.
function series(rows; arm, cv, sort_by = :f1)
    rs = [r for r in rows if r.arm == arm && isapprox(r.cv, cv; atol = 1e-8)]
    return sort(rs, by = r -> getproperty(r, sort_by))
end


# A panel with a missing curve is the failure this guards against: it looks like
# a finished figure and is simply wrong.  Name what is absent instead.
function require_series(rs, what, path)
    isempty(rs) && throw(ErrorException(
        "no rows for $what in $(relpath(path)).\n" *
        "Run the producing script first, or check --tag."))
    return rs
end

# The iso arm's fair-price ladder often runs further down than its loaded one,
# because the low precision cv33 solves are the last to land.  That is fine and
# expected - the panels plot everything available - but it should be visible in
# the log rather than discovered in the finished figure.
function report_range(label, rs)
    isempty(rs) && return
    println("    ", rpad(label, 34), length(rs), " points, F1 ",
            round(minimum(r.f1 for r in rs), digits = 3), " - ",
            round(maximum(r.f1 for r in rs), digits = 3))
end


#####################################################################
### Panels
#####################################################################

# Two series on one panel: a connecting line plus markers, styled so the two are
# told apart by fill rather than by colour (the manuscript figures are greyscale).
function two_series_panel(a, b; field, ylab, labels, fills, title, ylim = nothing)
    plt = nothing
    for (rs, label, fill) in ((a, labels[1], fills[1]), (b, labels[2], fills[2]))
        isempty(rs) && continue
        x = [r.f1 for r in rs]
        y = [getproperty(r, field) for r in rs]
        plt = plt === nothing ?
            Plots.plot(x, y, color = "grey", label = "") :
            Plots.plot!(x, y, color = "grey", label = "")
        Plots.scatter!(x, y, color = fill, markerstrokecolor = "black", label = label)
    end
    Plots.plot!(xlab = F1_LABEL, ylab = ylab, title = title,
                titlelocation = :left, titlefontsize = 14; PANEL...)
    ylim === nothing || Plots.plot!(ylim = ylim)
    return plt
end


function panel_a(utility, path)
    fair = require_series(series(utility; arm = "iso_premium", cv = 0.0),
                          "iso premium, 0% loading", path)
    cv33 = series(utility; arm = "iso_premium", cv = 0.33)
    println("  panel A  utility, risk averse base")
    report_range("iso premium, 0% loading", fair)
    report_range("iso premium, 50% loading", cv33)
    return two_series_panel(fair, cv33; field = :pct_increase_CE,
                            ylab = UTILITY_LABEL,
                            labels = ("loading factor 0%", "loading factor 50%"),
                            fills = (FILLED, OPEN), title = "A")
end


function panel_b(bankruptcy, path)
    fair = require_series(series(bankruptcy; arm = "iso_premium", cv = 0.0),
                          "iso premium, 0% loading", path)
    cv33 = series(bankruptcy; arm = "iso_premium", cv = 0.33)
    println("  panel B  bankruptcy, risk neutral base")
    report_range("iso premium, 0% loading", fair)
    report_range("iso premium, 50% loading", cv33)
    return two_series_panel(fair, cv33; field = :pct_reduction,
                            ylab = REDUCTION_LABEL,
                            labels = ("loading factor 0%", "loading factor 50%"),
                            fills = (FILLED, OPEN), title = "B", ylim = (0, 100))
end


function panel_c(bankruptcy, path)
    # Both series are 50% loading; the perfect index belongs to both arms and is
    # solved once, under the iso label, so it is prepended to the full coverage
    # series to give that curve its endpoint.
    iso  = require_series(series(bankruptcy; arm = "iso_premium", cv = 0.33),
                          "iso premium, 50% loading", path)
    full = series(bankruptcy; arm = "full_coverage", cv = 0.33)
    perfect = [r for r in iso if isapprox(r.f1, 1.0; atol = 1e-8)]
    full = sort(vcat(perfect, full), by = r -> r.f1)
    require_series(full, "full coverage, 50% loading", path)

    println("  panel C  bankruptcy, both arms at 50% loading")
    report_range("full coverage (recall = 100%)", full)
    report_range("iso premium (recall = precision)", iso)
    return two_series_panel(full, iso; field = :pct_reduction,
                            ylab = REDUCTION_LABEL,
                            labels = ("recall = 100%", "recall = precision"),
                            fills = (FILLED, OPEN), title = "C", ylim = (0, 100))
end


#####################################################################
### Main
#####################################################################

function main(args)
    opts = parse_args(args)

    println("=" ^ 70)
    println("three panel composite figure")
    println("utility     ", opts[:utility])
    println("bankruptcy  ", opts[:bankruptcy])
    println("output      ", opts[:out])
    println("layout      ", opts[:layout])
    println("dpi         ", opts[:dpi])
    println("=" ^ 70)

    for (what, path) in (("utility", opts[:utility]), ("bankruptcy", opts[:bankruptcy]))
        isfile(path) || begin
            println("\nmissing the $what table: ", path)
            println("run index_model/analysis/run_",
                    what == "utility" ? "utility" : "bankruptcy", "_figures.jl first")
            return 1
        end
    end

    utility    = read_csv(opts[:utility])
    bankruptcy = read_csv(opts[:bankruptcy])
    println("\nread ", length(utility), " utility rows, ",
            length(bankruptcy), " bankruptcy rows\n")

    pa = panel_a(utility,    opts[:utility])
    pb = panel_b(bankruptcy, opts[:bankruptcy])
    pc = panel_c(bankruptcy, opts[:bankruptcy])

    # Sizes are PRINT inches x 100: 1x3 is 9.0 x 2.8 in, a full page width row;
    # 3x1 is 3.6 x 8.4 in, a single column.  The y labels are two lines and the
    # x label sits below the tick labels, so the margins have to be generous or
    # Plots crops them at the panel edge - silently, and only visible in the
    # finished file.
    layout, size_ = opts[:layout] == "1x3" ? ((1, 3), (900, 280)) : ((3, 1), (360, 840))
    plt = Plots.plot(pa, pb, pc; layout = layout, size = size_, dpi = opts[:dpi],
                     left_margin = 7Plots.mm, bottom_margin = 6Plots.mm,
                     top_margin = 3Plots.mm, right_margin = 3Plots.mm)

    mkpath(dirname(abspath(opts[:out])))
    savefig_dpi(plt, opts[:out]; dpi = opts[:dpi])
    println("\nwrote ", relpath(opts[:out]), "  (", opts[:dpi], " dpi, ",
            join(size_ .* (opts[:dpi] ÷ 100), " x "), " px)")
    if opts[:pdf]
        pdf = replace(opts[:out], r"\.png$" => ".pdf")
        pdf == opts[:out] || (Plots.savefig(plt, pdf);
                              println("wrote ", relpath(pdf), "  (vector)"))
    end

    println("\nNOTE: panel A is the RISK AVERSE base; panels B and C are RISK NEUTRAL.")
    println("      Two channels by which basis risk bites, not one agent's three outcomes.")
    return 0
end


if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
