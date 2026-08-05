#####################################################################
#####################################################################
### Solves one index insurance scenario WITHOUT insurance (the premium
### grid collapses to [0.0]) and saves the solution to
### results_index/no_insurance_<base>/.
###
###     julia --project=<repo> --threads=N index_model/run_no_insurance.jl \
###           <scenario> [config.toml]
###
### Takes a scenario name and the same optional config path as run.jl,
### so the two cannot be invoked inconsistently (the base model's
### Makefile passes 1 argument to src/run_no_insurance.jl, which wants
### 7, and errors).  Nη is read and ignored: with no insurance the
### premium grid is a single point.
###
### Many scenarios map to ONE output directory, because the solution is
### invariant to both (pr, re) and (cf, cv) - see no_insurance_name in
### parameters/generate_scenarios.jl.  Submit this once per base
### parameter set, not once per scenario: concurrent jobs would race to
### write the same sol.jld2.
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

# See the note in run.jl: top level, not inside solve_scenario.  Defines `p`.
include(joinpath(@__DIR__, "parameters", args.params * ".jl"))

solve_scenario(p, args.params; insurance = false, args.sizes...)
