#####################################################################
### Reload a solved model saved by run.jl and re-make the solution plots.
###
### Usage (use the same grid arguments that were used for run.jl, minus maxiter):
###   julia -t 6 src/plot_solution.jl RA_intermidiate 10 10 12 12 6 7
###
### run.jl only saves the value and policy functions (sol.jld2), so the dynamic
### program is rebuilt with `init` on the same grids and the saved solution is
### loaded into it. The grids must match those of the original run. Rebuilding
### also recomputes the premium interpolation, which is the slow part.
#####################################################################
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Plots
params, Nκ, Nη, Ns_grid, Nb_grid, NΔb_grid, Nquad = ARGS
include("../src/model.jl"); include("../src/plotting.jl")
include("../parameters/V2/$params.jl"); include("../parameters/V2/insurer_capital_model.jl")

if !isfile("results/V2/$params/sol.jld2")
    error("no saved solution at results/V2/$params/sol.jld2, run src/run.jl first")
end

function main(params, Nκ, Nη, Ns_grid, Nb_grid, NΔb_grid, Nquad)
    X, Xmc = init_ranom_variables(p; Nquad = parse(Int,Nquad))

    # state grids, keep in sync with run.jl
    bmax = 0.5; bmin = -3.5; Nb_grid = parse(Int,Nb_grid)
    Δbmax = 1.25; Δbmin = -1.0; NΔb_grid = parse(Int,NΔb_grid)
    if occursin("short", params)
        bmax = 1.0; bmin = -3.0
        Δbmax = 2; Δbmin = -1.25
    elseif occursin("long", params)
        bmax = 0.0; bmin = -4.5
        Δbmax = 0.75; Δbmin = -0.75
    end
    b_grid = bmin:((bmax-bmin)/(Nb_grid-1)):bmax
    Δb_grid = Δbmin:((Δbmax-Δbmin)/(NΔb_grid-1)):Δbmax

    smax = 2.5; smin = -4.0; Ns_grid = parse(Int,Ns_grid)
    s_grid = smin:((smax-smin)/(Ns_grid-1)):smax
    state_grid = [s_grid,b_grid,Δb_grid]

    # action grid, keep in sync with run.jl
    κmin = 0.01; κmax = 0.99; Nκ = parse(Int,Nκ)
    κ_grid = κmin:((κmax-κmin)/(Nκ-1)):κmax
    ηmin = 0.0; ηmax=0.05; Nη = parse(Int,Nη)
    η_grid = ηmin:((ηmax-ηmin)/(Nη-1)):ηmax
    action_grid = [κ_grid,η_grid]

    function strike(s,p)
        return 0.5*(p.b_target + p.b_limit)
    end

    function index(s,X)
        return exp(s[1])+X[1]
    end

    # rebuild the program (without solving) and load the saved solution
    prob, El = init(action_grid,state_grid,p.δ,p,X,index,strike)
    load_solution(prob,"results/V2/$params/sol.jld2")

    ##################################
    ####  Plot solutions.         ####
    ##################################
    x=simulation(prob,X,Xmc,100;with_policy = false, s0 = [1.0,log(-p.s̄),log(0.4),0.0,1.0])
    plt = plot_simulation(x,max_USD=6.0,max_Bt=1.5,T = 95)
    savefig(plt,"results/V2/$params/simulation.png")

    plt = plot_premiums(prob)
    savefig(plt,"results/V2/$params/premiums.png")

    plt = plot_loading(El,b_grid,Δb_grid,X,prob.p,index,strike)
    savefig(plt,"results/V2/$params/loading.png")

    plt = plot_value_function(prob)
    savefig(plt,"results/V2/$params/value_function.png")
end
main(params,Nκ, Nη, Ns_grid, Nb_grid, NΔb_grid, Nquad)
