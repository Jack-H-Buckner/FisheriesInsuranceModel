using Pkg 
#Pkg.activate(".")
include("../src/model.jl")
include("../src/analysis.jl")
include("../src/plotting.jl")
print("here")
coverage_1 = coverage_levels("RA_intermidiate", smax = 2.5)
coverage_2 = coverage_levels("RA_intermidiate_higher_k", smax = 2.5)
coverage_3 = coverage_levels("RA_intermidiate_higher_k_2", smax = 2.5)

fig = plot_policy(coverage_1, coverage_2, coverage_3, "Optimal policies";
    titles = ("low interst (k = 0%)", "intermedite interst (k = 2%)", "high inteest (k = 4%)"))
CairoMakie.save("manuscript/figure_S6.png", fig)
fig