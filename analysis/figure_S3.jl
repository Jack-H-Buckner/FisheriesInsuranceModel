using Pkg
#Pkg.activate(".")
include("../src/model.jl")
include("../src/analysis.jl")
include("../src/plotting.jl")
print("here")
coverage_1 = coverage_levels("no_RA_short")
coverage_2 = coverage_levels("no_RA_intermidiate")
coverage_3 = coverage_levels("no_RA_long")

fig = plot_policy(coverage_1,coverage_2,coverage_3,"Optimal policies")
CairoMakie.save("manuscript/figure_S3.png", fig)
fig
