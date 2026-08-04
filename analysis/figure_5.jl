using Pkg

using Plots, CairoMakie, LaTeXStrings

Survival = exp.(-1*[0.4, 0.2, 0.11])
Value_Bankruptcy = [2.09, 3.46, 2.50]
Value_risk_averse = [2.56, 4.31, 9.10]

pltA = Plots.plot(Survival,Value_Bankruptcy, label = "Bankruptcy", color = "black", linewidth = 2,
    title = "A", titlelocation = :left)
Plots.scatter!(Survival,Value_Bankruptcy, label = "", color = "black", markersize = 5)
Plots.plot!(Survival,Value_risk_averse, label = "Risk Averse", color = "black", linestyle = :dash, linewidth = 2)
Plots.scatter!(Survival,Value_risk_averse, label = "", color = "white", markersize = 5)
Plots.plot!(size = (450,275), xlabel = "Survival", ylabel = "Percent increase in value")
Plots.annotate!([0.7, 0.875], [4, 4], [string("Fast LH"),string("Slow LH")], ylims = (0,9.5))
#savefig("~/Documents/value_by_LH.png")
pltA


using CSV, DataFrames
df = CSV.read("results/no_RA_intermidiate_lower_frequency/value_diff.csv",DataFrame)
df.frequency .= 0.01

df2 = CSV.read("results/no_RA_intermidiate/value_diff.csv",DataFrame)
df2.frequency .= 0.025
df = vcat(df,df2)

df3 = CSV.read("results/no_RA_intermidiate_higher_frequency/value_diff.csv",DataFrame)
df3.frequency .= 0.05
df = vcat(df,df3)

df3 = CSV.read("results/no_RA_intermidiate_higher_frequency_2/value_diff.csv",DataFrame)
df3.frequency .= 0.075
df = vcat(df,df3)

df4 = CSV.read("results/no_RA_intermidiate_higher_frequency_3/value_diff.csv",DataFrame)
df4.frequency .= 0.1
df = vcat(df,df4)


df_RA = CSV.read("results/RA_intermidiate_lower_frequency/value_diff.csv",DataFrame)
df_RA.frequency .= 0.01

df2 = CSV.read("results/RA_intermidiate/value_diff.csv",DataFrame)
df2.frequency .= 0.025
df_RA = vcat(df_RA,df2)

df3 = CSV.read("results/RA_intermidiate_higher_frequency/value_diff.csv",DataFrame)
df3.frequency .= 0.05
df_RA = vcat(df_RA,df3)

df3 = CSV.read("results/RA_intermidiate_higher_frequency_2/value_diff.csv",DataFrame)
df3.frequency .= 0.075
df_RA = vcat(df_RA,df3)

df4 = CSV.read("results/RA_intermidiate_higher_frequency_3/value_diff.csv",DataFrame)
df4.frequency .= 0.1
df_RA = vcat(df_RA,df4)

pltB = Plots.plot(df.frequency[df.savings .== -1.5],df.increase[df.savings .== -1.5],
            label = "", color = "black", linewidth = 1,
            title = "B", titlelocation = :left)
Plots.scatter!(df.frequency[df.savings .== -1.5],df.increase[df.savings .== -1.5],
            label = "", color = "black", markersize= 5, ylims = (0.0,4.0))
Plots.plot!(df_RA.frequency[df_RA.savings .== -1.5],df_RA.increase[df_RA.savings .== -1.5],
            label = "", color = "black", linewidth = 1, linestyle = :dash)
Plots.scatter!(df_RA.frequency[df_RA.savings .== -1.5],df_RA.increase[df_RA.savings .== -1.5],
            label = "", color = "white", markersize= 5, ylims = (0.0,11.0))
Plots.plot!(size = (450,275), xlabel = "Frequency", ylabel = "Percent increase in value")
#savefig("~/Documents/value_by_risk.png")
pltB

# Panel C: value of insurance under different harvest control rules
# (intermediate life history, savings = -1.5 to match panel B's slice).
# Two HCRs are available: baseline (b_limit = 0.20) and a more permissive
# rule (b_limit = 0.01, the "lower_blimit" parameter set). For each HCR
# we read the percent value increase at savings = -1.5 for both the
# bankruptcy and risk-averse utility specifications.
df_HCR_base_bank   = CSV.read("results/no_RA_intermidiate/value_diff.csv",              DataFrame)
df_HCR_lower_bank  = CSV.read("results/no_RA_intermidiate_lower_blimit/value_diff.csv", DataFrame)
df_HCR_base_RA     = CSV.read("results/RA_intermidiate/value_diff.csv",                 DataFrame)
df_HCR_lower_RA    = CSV.read("results/RA_intermidiate_lower_blimit/value_diff.csv",    DataFrame)

value_at(df, s) = df.increase[df.savings .== s][1]

# Order points by b_limit so the x-axis is monotone increasing (more
# permissive HCR on the left, more precautionary HCR on the right).
b_limits  = [0.01, 0.20]
bank_vals = [value_at(df_HCR_lower_bank, -1.5), value_at(df_HCR_base_bank, -1.5)]
RA_vals   = [value_at(df_HCR_lower_RA,   -1.5), value_at(df_HCR_base_RA,   -1.5)]

pltC = Plots.plot(b_limits, bank_vals, label = "", color = "black", linewidth = 1,
    title = "C", titlelocation = :left)
Plots.scatter!(b_limits, bank_vals, label = "", color = "black", markersize = 5)
Plots.plot!(b_limits, RA_vals, label = "", color = "black", linestyle = :dash, linewidth = 1)
Plots.scatter!(b_limits, RA_vals, label = "", color = "white", markersize = 5)
Plots.plot!(size = (450,275), xlabel = "Limit reference point (B/B0)", ylabel = "Percent increase in value")
pltC

# Shared y-axis range across all three panels so they're visually
# comparable. Take the max across every series plotted in any panel and
# add 10% headroom; pin the floor at 0 since "percent increase in value"
# is non-negative across all scenarios shown.
all_values = vcat(
    Value_Bankruptcy, Value_risk_averse,
    df.increase[df.savings .== -1.5], df_RA.increase[df_RA.savings .== -1.5],
    bank_vals, RA_vals,
)
ymax_global = 1.1 * maximum(all_values)

Plots.plot!(pltA, ylims = (0, ymax_global))
Plots.plot!(pltB, ylims = (0, ymax_global))
Plots.plot!(pltC, ylims = (0, ymax_global))

# Combined panel A + panel B + panel C figure for the manuscript
fig = Plots.plot(pltA, pltB, pltC, layout = (1, 3), size = (1400, 360),
    titlelocation = :left, titlefontsize = 14,
    left_margin = 10Plots.mm, bottom_margin = 8Plots.mm,
    right_margin = 5Plots.mm, top_margin = 4Plots.mm)
savefig(fig, "manuscript/figure_5.png")
fig

