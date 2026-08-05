#####################################################################
#####################################################################
### Shared solve driver for the index insurance model.
###
### run.jl and run_no_insurance.jl are thin wrappers over
### solve_scenario, so the two cannot drift apart.  In the base model
### src/run.jl, src/run_no_insurance.jl and src/run_variable_strike.jl
### are three near identical copies, which is how their argument
### contracts ended up inconsistent.
###
### The parameter set is passed in rather than included here: loading a
### parameter file inside a function would add its methods in a newer
### world age than the running caller.  The run scripts include it at
### top level and hand `p` over.
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

const RUNNER_USAGE = "usage: <params> [config.toml]\n" *
    "grid sizes are read from index_model/config.toml, not the command line"

# Positional CLI contract, shared by both runners:
#     <params> [config.toml]
# Unlike src/run.jl the sizes are not arguments; the optional second argument is
# a path to a config file to use instead of index_model/config.toml, which is
# how a test or a resolution check runs without editing the committed config.
function parse_runner_args(args)
    length(args) in (1, 2) || throw(ArgumentError(
        "expected 1 or 2 arguments, got $(length(args)): $(join(args, " "))\n" *
        RUNNER_USAGE))
    params = String(args[1])
    config = length(args) == 2 ? String(args[2]) : GRID_CONFIG_PATH
    return (params = params, config = config, sizes = read_grid_config(config))
end


results_root() = joinpath(normpath(joinpath(@__DIR__, "..")), "results_index")

# The insurance solve is keyed on the scenario; the no insurance solve collapses
# to one directory per base parameter set (see no_insurance_name).
output_dir(params; insurance) =
    joinpath(results_root(), insurance ? params : no_insurance_name(params))


"""
    solve_scenario(p, params; insurance = true, plots = true, sizes...)

Solve one index insurance scenario and write the solution, the grid metadata
and the diagnostic plots to `results_index/`.

`p` is the parameter set defined by `index_model/parameters/\$params.jl`, loaded
by the caller at top level.
"""
function solve_scenario(p, params::AbstractString;
                        insurance::Bool = true,
                        plots::Bool = true,
                        maxiter = nothing,
                        outdir = nothing,
                        sizes...)

    g = build_grids(params; insurance = insurance, sizes...)
    outdir = something(outdir, output_dir(params; insurance = insurance))
    mkpath(outdir)   # mkpath, not mkdir: race safe when jobs start together

    println("=" ^ 70)
    println("scenario   ", params)
    println("mode       ", insurance ? "with insurance" : "no insurance")
    println("output     ", outdir)
    println("threads    ", Threads.nthreads())
    println("=" ^ 70)
    print_grid(g)
    println()
    print_index_diagnostics(p.T)
    println()

    X, Xmc = init_ranom_variables(p; Nquad = g.sizes.Nquad)
    prob = init(g.action_grid, g.state_grid, g.δ, p, X, index)

    maxiter = something(maxiter, round(Int, 5 / (1 - g.δ)))
    println("solving (maxiter = ", maxiter, ") ...")
    t0 = time()
    solve!(prob, maxiter = maxiter)
    elapsed = time() - t0
    println("solved in ", round(elapsed / 60, digits = 2), " min")

    # Save before plotting: a plotting failure must not throw away a solve that
    # may have taken hours.
    save_solution(prob, joinpath(outdir, "sol.jld2"))
    write_grid_meta(joinpath(outdir, "grid.toml"), g, p, params)
    println("wrote ", joinpath(outdir, "sol.jld2"), " and grid.toml")

    if plots
        try
            diagnostic_plots(prob, p, g, X, Xmc, outdir)
            println("wrote diagnostic plots")
        catch err
            @warn "diagnostic plots failed; the solution is still saved" exception = err
        end
    end

    return prob, g, outdir
end


# Grid coverage checks and solution summaries, mirroring src/run.jl:35-102.
# check_b_grid is the one that matters: it confirms the state grid actually
# covers the ergodic set of the biological model for this parameter set.
function diagnostic_plots(prob, p, g, X, Xmc, outdir)
    s_grid, b_grid, Δb_grid = g.state_grid
    κ_grid, η_grid = g.action_grid

    states, m, f, h, y = simulate_trajectory(p, X, 5000)
    Plots.scatter(states[1, :], states[2, :], label = "simulated",
                  xlabel = "log biomass", ylabel = "Δ log biomass", markersize = 2)
    Plots.scatter!(b_grid, zeros(length(b_grid)), label = "b grid")
    plt = Plots.scatter!(zeros(length(Δb_grid)), Δb_grid, label = "Δb grid")
    Plots.savefig(plt, joinpath(outdir, "check_b_grid.png"))

    plt = Plots.scatter(exp.(s_grid), s_grid, label = "Savings",
                        xlabel = "Savings", ylabel = "log savings")
    Plots.savefig(plt, joinpath(outdir, "check_s_grid.png"))

    u = [(κ, η) for κ in κ_grid, η in η_grid if κ + η <= 1]
    plt = Plots.scatter(first.(u), last.(u), xlabel = "Consumption",
                        ylabel = "Premiums", label = "")
    Plots.savefig(plt, joinpath(outdir, "check_u_grid.png"))

    plt = plot_value_function(prob)
    Plots.savefig(plt, joinpath(outdir, "value_function.png"))

    x = simulation(prob, El, Xmc, 76; with_policy = false,
                        s0 = [1.0, 0.0, log(0.3), 0.0, 1.0])
    plt = plot_simulation(x, p, max_USD = 6.0, max_Bt = 1.5, T = 75)
    Plots.savefig(plt, joinpath(outdir, "simulation.png"))
    return nothing
end
