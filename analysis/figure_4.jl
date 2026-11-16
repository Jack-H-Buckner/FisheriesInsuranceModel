using Pkg 
#Pkg.activate(".")
include("../src/model.jl")
include("../src/analysis.jl")
include("../src/plotting.jl")
# set_png_dpi!: stamp the PNG pHYs chunk so the file reports the resolution it
# was rendered at instead of defaulting to 72.  Guarded so that running several
# figure scripts in one session does not redefine its consts.
isdefined(Main, :savefig_dpi) || include("../index_model/analysis/png_dpi.jl")

coverage_1 = coverage_levels("RA_short")
coverage_2 = coverage_levels("RA_intermidiate")
coverage_3 = coverage_levels("RA_long")

# plot_policy_gr rather than plot_policy: same figure, drawn with Plots/GR
# instead of CairoMakie, which is not a dependency of this project and costs
# more in compile time than the whole rest of this script does to run.  The
# resolution is set inside plot_policy_gr (9 x 5 in at 300 dpi, 2700 x 1500 px).
fig = plot_policy_gr(coverage_1, coverage_2, coverage_3, "Optimal policies")

savefig_dpi(fig, "manuscript/figure_4.png"; dpi = 300)
Plots.savefig(fig, "manuscript/figure_4.pdf")   # vector companion
fig