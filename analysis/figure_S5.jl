using Pkg 
#Pkg.activate(".")
include("../src/model.jl")
include("../src/analysis.jl")
include("../src/plotting.jl")
print("here")
coverage_1 = coverage_levels("RA_short",smax = 2.5)
coverage_2 = coverage_levels("RA_intermidiate",smax = 2.5)
coverage_3 = coverage_levels("RA_long",smax = 2.5)

fig = plot_policy(coverage_1,coverage_2,coverage_3,"Optimal policies")
CairoMakie.save("manuscript/figure_S5.png", fig)
fig