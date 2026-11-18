# Mirrors run.jl for the model with disaster relief (model_with_relief.jl): the discrete
# state is (Ib, ID, IM) with 8 levels, and p gets the declaration and relief probabilities
# pR from parameters/V2/disaster_relief.jl (p.Relief comes from the parameter file).
# A disaster is always declared after a high mortality event, so there is no pD.
# Usage: julia src/run_with_relief.jl <params> Nκ Nη Ns_grid Nb_grid NΔb_grid Nquad maxiter
using Pkg 
Pkg.activate(".")
Pkg.instantiate()
using Plots
params, Nκ, Nη, Ns_grid, Nb_grid, NΔb_grid, Nquad, maxiter = ARGS
include("../src/model_with_relief.jl"); include("../src/plotting.jl")
include("../parameters/V2/$params.jl"); include("../parameters/V2/insurer_capital_model.jl")
include("../parameters/V2/disaster_relief.jl") # pR
p = ComponentArray(merge(NamedTuple(p), (pR = pR,)))

# plotting.jl calls the two argument encode_discrete(Ib, IM); with relief the plots show
# the slice with no disaster declared (ID = 0).
encode_discrete(Ib,IM) = encode_discrete(Ib,0,IM)

outdir = "results/V2/relief_$params"
if !isdir(outdir)
    mkpath(outdir)
end 

function main(params,Nκ, Nη, Ns_grid, Nb_grid, NΔb_grid, Nquad,maxiter)
    ###################################
    ####  Define random variables #####
    ###################################
    X, Xmc = init_ranom_variables_relief(p; Nquad = parse(Int,Nquad))

    # define state variables and check the grid 
    ######################################
    ####  Define state varible grids #####
    ######################################
    # check paramters to get min and max values
    bmax = 0.5; bmin = -3.5; Nb_grid = parse(Int,Nb_grid)
    Δbmax = 1.25; Δbmin = -1.0; NΔb_grid = parse(Int,NΔb_grid)
    if occursin("short", params)
        bmax = 1.0; bmin = -3.0
        Δbmax = 2; Δbmin = -1.25
    elseif occursin("long", params)
        bmax = 0.0; bmin = -4.5
        Δbmax = 0.75; Δbmin = -0.75
    end 
    # stock 

    states, m, f, h, y = simulate_trajectory(p,X,5000)
    Plots.scatter(states[1,:], states[2,:])
    b_grid = bmin:((bmax-bmin)/(Nb_grid-1)):bmax
    Δb_grid = Δbmin:((Δbmax-Δbmin)/(NΔb_grid-1)):Δbmax
    plt = Plots.scatter!(b_grid,[0.0])
    Plots.scatter!(zeros(length(Δb_grid)),Δb_grid)
    savefig(plt,"$outdir/check_b_grid.png")

    # savings 
    smax = 2.5; smin = -4.0; Ns_grid = parse(Int,Ns_grid)
    s_grid = smin:((smax-smin)/(Ns_grid-1)):smax
    plt = Plots.scatter(exp.(s_grid),s_grid, label = "Savings", xlabel = "Savings", ylabel = "log savings")
    savefig(plt,"$outdir/check_s_grid.png")
    state_grid = [s_grid,b_grid,Δb_grid]

    ###############################
    ####  Define action grid  #####
    ###############################
    κmin = 0.01; κmax = 0.99; Nκ = parse(Int,Nκ)
    κ_grid = κmin:((κmax-κmin)/(Nκ-1)):κmax
    ηmin = 0.0; ηmax = η_max; Nη = parse(Int,Nη)
    η_grid = ηmin:((ηmax-ηmin)/(Nη-1)):ηmax
    action_grid = [κ_grid,η_grid]

    # plot grid
    k=0; u_grid = zeros(2,length(κ_grid)*length(η_grid))
    for (i,j) in Iterators.product(1:length(κ_grid),1:length(η_grid))
        k +=1
        u_grid[:,k] = vcat(κ_grid[i],η_grid[j]) 
    end
    u_grid=u_grid[:,reshape(mapslices(sum,u_grid,dims = 1).<=1,length(κ_grid)*length(η_grid))]
    plt = Plots.scatter(u_grid[1,:],u_grid[2,:], xlabel = "Consumption", ylabel= "Premiums", label = "")
    savefig(plt,"$outdir/check_u_grid.png")


    #######################################
    ####  Strike and index functions  #####
    #######################################

    function strike(s,p)
        return 0.5*(p.b_target + p.b_limit)
    end

    function index(s,X)
        return exp(s[1])+X[1] 
    end 

    ##################################
    ####  Set up dynamic program  ####
    ##################################
    δ = p.δ
    prob, El = init(action_grid,state_grid,δ,p,X,index,strike)
    ValueFunctionIterations.solve!(prob,maxiter = parse(Int,maxiter)) #round(Int,5/(1-δ)))
    save_solution(prob,"$outdir/sol.jld2")


    ##################################
    ####  Plot solutions.         ####
    ##################################
    # start solvent, no disaster declared, typical mortality
    x=simulation(prob,X,Xmc,250;with_policy = false, s0 = [encode_discrete(0,0,0),0.0,log(0.3),0.0])
    plt = plot_simulation(x,max_USD=6.0,max_Bt=1.5,T = 200)
    savefig(plt,"$outdir/simulation.png")

    plt = plot_premiums(prob)
    savefig(plt,"$outdir/premiums.png")

    plt = plot_loading(El,b_grid,Δb_grid,X,prob.p,index,strike)
    savefig(plt,"$outdir/loading.png")

    plt = plot_value_function(prob)
    savefig(plt,"$outdir/value_function.png")
    
end
main(params,Nκ, Nη, Ns_grid, Nb_grid, NΔb_grid, Nquad,maxiter)
