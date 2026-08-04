using Pkg
using Plots, CSV, DataFrames

# Helper: read a value_diff.csv and tag it with a parameter value.
function load_diff(path, col, val)
    df = CSV.read(path, DataFrame)
    df[!, col] .= val
    return df
end

# Panel A: Non-risk-averse, value increase vs fixed costs (y0 variants).
df_A = vcat(
    load_diff("results/no_RA_intermidiate_lower_y0/value_diff.csv",  :fixed_costs, 0.4),
    load_diff("results/no_RA_intermidiate/value_diff.csv",           :fixed_costs, 0.25),
    load_diff("results/no_RA_intermidiate_higher_y0/value_diff.csv", :fixed_costs, 0.1),
)

# Panel B: Non-risk-averse, value increase vs interest rate (k variants).
df_B = vcat(
    load_diff("results/no_RA_intermidiate/value_diff.csv",            :interest, 0.0),
    load_diff("results/no_RA_intermidiate_higher_k/value_diff.csv",   :interest, 0.02),
    load_diff("results/no_RA_intermidiate_higher_k_2/value_diff.csv", :interest, 0.04),
    load_diff("results/no_RA_intermidiate_higher_k_3/value_diff.csv", :interest, 0.06),
)

# Panel C: Risk-averse, value increase vs interest rate (k variants).
df_C = vcat(
    load_diff("results/RA_intermidiate/value_diff.csv",            :interest, 0.0),
    load_diff("results/RA_intermidiate_higher_k/value_diff.csv",   :interest, 0.02),
    load_diff("results/RA_intermidiate_higher_k_2/value_diff.csv", :interest, 0.04),
    load_diff("results/RA_intermidiate_higher_k_3/value_diff.csv", :interest, 0.06),
)

# Plot three lines per panel for savings = 1.5, 0.0, -1.5.
function add_savings_lines!(plt, df, xcol; show_legend = false)
    s_high, s_mid, s_low = 1.5, 0.0, -1.5
    sub(s) = sort!(df[df.savings .== s, :], xcol)

    d = sub(s_high)
    Plots.plot!(plt, d[!, xcol], d.increase,
        label = "", color = "black", linewidth = 1, linestyle = :dash)
    Plots.scatter!(plt, d[!, xcol], d.increase,
        label = show_legend ? "savings = 1.5"  : "",
        color = "black", markersize = 5)

    d = sub(s_mid)
    Plots.plot!(plt, d[!, xcol], d.increase,
        label = "", color = "black", linewidth = 1)
    Plots.scatter!(plt, d[!, xcol], d.increase,
        label = show_legend ? "savings = 0.0"  : "",
        color = "grey", markersize = 5)

    d = sub(s_low)
    Plots.plot!(plt, d[!, xcol], d.increase,
        label = "", color = "black", linewidth = 1, linestyle = :dot)
    Plots.scatter!(plt, d[!, xcol], d.increase,
        label = show_legend ? "savings = -1.5" : "",
        color = "white", markersize = 5)
    return plt
end

ylims_A = (0.0, 4.5)
ylims_B = (-0.1, 6.5)
ylims_C = (0.0, 6.5)

pltA = Plots.plot(; title = "A", titlelocation = :left,
    xlabel = "Fixed costs", ylabel = "Percent increase in value",
    ylims = ylims_A, legend = :bottomright,
    legend_foreground_color = :transparent)
add_savings_lines!(pltA, df_A, :fixed_costs; show_legend = true)

pltB = Plots.plot(; title = "B", titlelocation = :left,
    xlabel = "Interest rate", ylabel = "Percent increase in value",
    ylims = ylims_B, legend = :none)
add_savings_lines!(pltB, df_B, :interest; show_legend = false)

pltC = Plots.plot(; title = "C", titlelocation = :left,
    xlabel = "Interest rate", ylabel = "Percent increase in value",
    ylims = ylims_C, legend = :none)
add_savings_lines!(pltC, df_C, :interest; show_legend = false)

fig = Plots.plot(pltA, pltB, pltC, layout = (1, 3), size = (1500, 400),
    titlelocation = :left, titlefontsize = 14,
    left_margin = 10Plots.mm, bottom_margin = 10Plots.mm,
    top_margin = 5Plots.mm, right_margin = 5Plots.mm)
savefig(fig, "manuscript/figure_S8.png")
fig
