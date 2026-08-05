#####################################################################
#####################################################################
### Solves one index insurance scenario WITH insurance and saves the
### solution to results_index/<scenario>/.
###
### Run from anywhere (paths are derived from @__DIR__):
###     julia --project=<repo> --threads=N index_model/run.jl <scenario> \
###           [config.toml]
###
### Grid sizes come from index_model/config.toml, not the command line,
### so every job in a sweep solves on the same grid.  The optional
### second argument points at a different config file.
### Run `make setup_index` once before submitting jobs: Pkg.instantiate
### is deliberately NOT called here, because dozens of workers racing on
### the depot to fetch the unregistered ValueFunctionIterations.jl is a
### real failure mode.
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

ENV["GKSwstype"] = "100" # headless GR, so plotting works on a compute node
using Plots, LaTeXStrings

include(joinpath(@__DIR__, "model.jl"))
include(joinpath(@__DIR__, "plotting.jl"))
include(joinpath(@__DIR__, "grids.jl"))
include(joinpath(@__DIR__, "solve.jl"))
include(joinpath(@__DIR__, "parameters", "generate_scenarios.jl"))

args = parse_runner_args(ARGS)   # also reads and validates the grid config
scenario_row(args.params)        # fail fast on a typo, before loading anything
println("grid sizes from ", args.config)

# Included at top level, not inside solve_scenario: a parameter file loaded
# during a function call would add its methods in a newer world age than the
# running caller.  This defines `p`.
include(joinpath(@__DIR__, "parameters", args.params * ".jl"))

solve_scenario(p, args.params; insurance = true, args.sizes...)
