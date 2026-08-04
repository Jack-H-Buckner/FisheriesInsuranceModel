using DataFrames, CSV, Plots
using Plots.Measures: mm

# Plot 10-year bankruptcy probabilities (P10) from table 4 as a grouped bar
# chart: one group per life history, two bars per group (with / without
# insurance), with whiskers showing the 95% bootstrap CI.

df = CSV.read("manuscript/table_4.csv", DataFrame)

# Order and human-readable labels for the life-history cases (match the
# naming used in figure 3).
cases  = ["no_RA_short", "no_RA_intermidiate", "no_RA_long"]
labels = ["Short", "Intermediate", "Long"]
n      = length(cases)

row_of(case) = findfirst(==(case), df.case)
col_for(case, col) = df[row_of(case), col]

# Convert proportions to percentages for the axis.
pct(x) = 100 * x

P10_with    = [pct(col_for(c, :P10_with_insurance)) for c in cases]
P10_with_lo = [pct(col_for(c, :P10_with_lower))     for c in cases]
P10_with_hi = [pct(col_for(c, :P10_with_upper))     for c in cases]
P10_no      = [pct(col_for(c, :P10_no_insurance))   for c in cases]
P10_no_lo   = [pct(col_for(c, :P10_no_lower))       for c in cases]
P10_no_hi   = [pct(col_for(c, :P10_no_upper))       for c in cases]

# Plots.jl yerror expects (distance_below, distance_above) from each bar top.
err_with = (P10_with .- P10_with_lo, P10_with_hi .- P10_with)
err_no   = (P10_no   .- P10_no_lo,   P10_no_hi   .- P10_no)

# Side-by-side bar layout: groups at integer x, arms offset by ±0.20.
group_x = collect(1:n)
offset  = 0.20
bar_w   = 0.36

ymax = 1.1 * maximum(vcat(P10_no_hi, P10_with_hi))

fig = Plots.bar(group_x .- offset, P10_no;
    yerror     = err_no,
    bar_width  = bar_w,
    color      = :grey,
    linecolor  = :black,
    label      = "No insurance",
    ylabel     = "10-year bankruptcy probability (%)",
    xlabel     = "Life history",
    xticks     = (group_x, labels),
    ylims      = (0, ymax),
    legend     = :topleft,
    framestyle = :box,
    size       = (700, 450),
    left_margin   = 8mm,
    bottom_margin = 6mm,
    right_margin  = 4mm,
    top_margin    = 3mm,
    tickfontsize  = 10,
    guidefontsize = 12,
    legendfontsize = 10,
)
Plots.bar!(fig, group_x .+ offset, P10_with;
    yerror    = err_with,
    bar_width = bar_w,
    color     = :white,
    linecolor = :black,
    label     = "With insurance",
)

Plots.savefig(fig, "manuscript/figure_6.png")
fig
