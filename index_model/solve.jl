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
    # The output tag comes from the config too, never from its own argument, so
    # the sizes and the directory they land in cannot be chosen independently.
    return (params = params, config = config,
            sizes = read_grid_config(config), tag = read_output_tag(config))
end


results_root() = joinpath(normpath(joinpath(@__DIR__, "..")), "results_index")

# The insurance solve is keyed on the scenario; the no insurance solve collapses
# to one directory per base parameter set (see no_insurance_name).
#
# `tag` is the suffix declared by the config the solve ran on: production solves
# use "" and keep their bare names, while a preliminary solve on
# config_cheaper.toml lands in results_index/<scenario>_cheaper/ and so cannot
# overwrite the production solution.
output_dir(params; insurance, tag::AbstractString = "") =
    joinpath(results_root(), (insurance ? params : no_insurance_name(params)) * tag)


# A solve overwrites whatever is already in its output directory, and rewrites
# grid.toml alongside, so the result stays self consistent and assert_grid_match
# has nothing to catch.  That is fine when re-solving on the same grid, and is
# how a preliminary and a production solve are kept apart (different tags, so
# different directories).  It is NOT fine when it happens by accident, so say so
# loudly.  A warning rather than an error: re-solving after a deliberate config
# change is legitimate, and a job that dies here has wasted a queue slot.
function warn_on_regrid(outdir, g)
    meta_path = joinpath(outdir, "grid.toml")
    isfile(meta_path) || return nothing
    old = try
        read_grid_meta(meta_path)
    catch
        return nothing            # unreadable metadata is not this check's job
    end
    recorded = get(old, "sizes", nothing)
    recorded isa AbstractDict || return nothing

    changed = [string(k) for (k, v) in pairs(g.sizes) if get(recorded, string(k), nothing) != v]
    isempty(changed) && return nothing

    @warn """
          overwriting an existing solution that was solved on a DIFFERENT grid.
          If this was not intended, the usual cause is running a preliminary
          config without its [output] tag, which puts cheap results in a
          production directory.""" dir = outdir differing = changed
    return nothing
end


"""
    solve_scenario(p, params; insurance = true, plots = true, tag = "", sizes...)

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
                        tag::AbstractString = "",
                        sizes...)

    g = build_grids(params; insurance = insurance, sizes...)
    outdir = something(outdir, output_dir(params; insurance = insurance, tag = tag))
    warn_on_regrid(outdir, g)
    mkpath(outdir)   # mkpath, not mkdir: race safe when jobs start together

    println("=" ^ 70)
    println("scenario   ", params)
    println("mode       ", insurance ? "with insurance" : "no insurance")
    println("output     ", outdir)
    isempty(tag) || println("tag        ", tag, "  (preliminary; not a production result)")
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
    write_grid_meta(joinpath(outdir, "grid.toml"), g, p, params; tag = tag)
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
