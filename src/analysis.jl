using Pkg, Random
using Plots, DataFrames, CSV
include("../src/model.jl"); include("../src/plotting.jl")
include("../src/insurnace_contracts.jl")
function coverage_levels(params; smax = log(2.0))
    #params = "no_RA_short"
    Nκ = "20"; Nη = "20"; Ns_grid = "20"; Nb_grid = "20"; NΔb_grid = "3"; Nquad = "10"
    include("../parameters/$params.jl")

    ###################################
    ####  Define random variables #####
    ###################################
    X, Xmc = init_ranom_variables(p; Nquad = parse(Int,Nquad))

    # define state variables and check the grid 
    ######################################
    ####  Define state varible grids #####
    ######################################
    # check paramters to get min and max values
    bmax = 0.0; bmin = -2.5; Nb_grid = parse(Int,Nb_grid)
    Δbmax = 0.5; Δbmin = -0.5; NΔb_grid = parse(Int,NΔb_grid)

    # Stock 
    Random.seed!(1234)
    b_grid = bmin:((bmax-bmin)/(Nb_grid-1)):bmax
    Δb_grid = Δbmin:((Δbmax-Δbmin)/(NΔb_grid-1)):Δbmax


    # Savings 
    smin = -2.5; Ns_grid = parse(Int,Ns_grid)
    s_grid = smin:((smax-smin)/(Ns_grid-1)):smax
    state_grid = [s_grid,b_grid,Δb_grid]

    ###############################
    ####  Define action grid  #####
    ###############################
    κmin = 0.01; κmax = 0.99; Nκ = parse(Int,Nκ)
    κ_grid = κmin:((κmax-κmin)/(Nκ-1)):κmax
    ηmin = 0.0; ηmax=0.2; Nη = parse(Int,Nη)
    η_grid = ηmin:((ηmax-ηmin)/(Nη-1)):ηmax
    action_grid = [κ_grid,η_grid]

    # plot grid
    k=0; u_grid = zeros(2,length(κ_grid)*length(η_grid))
    for (i,j) in Iterators.product(1:length(κ_grid),1:length(η_grid))
        k +=1
        u_grid[:,k] = vcat(κ_grid[i],η_grid[j]) 
    end
    u_grid=u_grid[:,reshape(mapslices(sum,u_grid,dims = 1).<=1,length(κ_grid)*length(η_grid))]


    #######################################
    ####  Strike and index functions  #####
    #######################################
    function strike(s,p)
        return 0.5*(p.b_target +p.b_limit) 
    end

    function index(s,X)
        return exp(s[1])+X[1]
    end 

    ##################################
    ####  Set up dynamic program  ####
    ##################################
    δ = 0.95
    prob, El = init(action_grid,state_grid,δ,p,X,index,strike)
    states = prob.V.states[:,prob.V.states[1,:].==1]

    dims = prob.V.Bslines[1].dims

    ##################################
    ####  load in solutions       ####
    ##################################
    load_solution(prob,"results/$params/sol.jld2")

    ##################################
    ####  calcaulte policy        ####
    ##################################

    u_opt = ValueFunctionIterations.policy(states,prob.u,prob.X,prob.R,prob.F,prob.p,prob.δ,prob.V)
    η_s = reshape(u_opt[2,:],dims...)
    η̄_s =  mapslices(s -> price_per_exposure(El(s[3],s[4],s[5]),prob.p.cf,prob.p.cv),states,dims = 1)
    η̄_s = reshape(η̄_s,dims...)
    u = [1e-10,0.0]
    function insurance_claims(s)
        nodes = prob.X(s,prob.p).nodes
        weights = prob.X(s,prob.p).weights
        weights = weights[nodes[5,:].==1.0] 
        weights = weights ./ sum(weights)
        nodes = nodes[:,nodes[5,:].==1.0] 
        claim = mapslices(x -> index_insurance(index(prob.F(s,u,x,prob.p)[3:5],[0.0]),strike(s,prob.p),1,1), nodes, dims = 1)
        sum(claim' .* weights )
    end 
    claims =  mapslices(s -> insurance_claims(s),states,dims = 1)
    claims = reshape(claims,dims...)

    return η_s, η̄_s, claims, state_grid
end 




function coverage_levels_variable_strike(params; smax = log(2.0))
    #params = "no_RA_short"
    Nκ = "20"; Nη = "20"; Ns_grid = "20"; Nb_grid = "20"; NΔb_grid = "3"; Nquad = "10"
    include("../parameters/$params.jl")

    ###################################
    ####  Define random variables #####
    ###################################
    X, Xmc = init_ranom_variables(p; Nquad = parse(Int,Nquad))

    # define state variables and check the grid 
    ######################################
    ####  Define state varible grids #####
    ######################################
    # check paramters to get min and max values
    bmax = 0.0; bmin = -2.5; Nb_grid = parse(Int,Nb_grid)
    Δbmax = 0.5; Δbmin = -0.5; NΔb_grid = parse(Int,NΔb_grid)

    # Stock 
    Random.seed!(1234)
    b_grid = bmin:((bmax-bmin)/(Nb_grid-1)):bmax
    Δb_grid = Δbmin:((Δbmax-Δbmin)/(NΔb_grid-1)):Δbmax

    # Savings 
    smin = -2.5; Ns_grid = parse(Int,Ns_grid)
    s_grid = smin:((smax-smin)/(Ns_grid-1)):smax
    state_grid = [s_grid,b_grid,Δb_grid]

    ###############################
    ####  Define action grid  #####
    ###############################
    κmin = 0.01; κmax = 0.99; Nκ = parse(Int,Nκ)
    κ_grid = κmin:((κmax-κmin)/(Nκ-1)):κmax
    ηmin = 0.0; ηmax=0.2; Nη = parse(Int,Nη)
    η_grid = ηmin:((ηmax-ηmin)/(Nη-1)):ηmax
    action_grid = [κ_grid,η_grid]

    # plot grid
    k=0; u_grid = zeros(2,length(κ_grid)*length(η_grid))
    for (i,j) in Iterators.product(1:length(κ_grid),1:length(η_grid))
        k +=1
        u_grid[:,k] = vcat(κ_grid[i],η_grid[j]) 
    end
    u_grid=u_grid[:,reshape(mapslices(sum,u_grid,dims = 1).<=1,length(κ_grid)*length(η_grid))]


    #######################################
    ####  Strike and index functions  #####
    #######################################
    function strike(s,p)
        s,h,f = update_stock(s,[0.0,1.0],p)
        return 1.1*exp(s[1])
    end

    function index(s,X)
        return exp(s[1])+X[1]
    end 

    ##################################
    ####  Set up dynamic program  ####
    ##################################
    δ = 0.95
    prob, El = init(action_grid,state_grid,δ,p,X,index,strike)
    states = prob.V.states[:,prob.V.states[1,:].==1]

    dims = prob.V.Bslines[1].dims

    ##################################
    ####  load in solutions       ####
    ##################################
    load_solution(prob,"results/variable_strike_$params/sol.jld2")

    ##################################
    ####  calcaulte policy        ####
    ##################################

    u_opt = ValueFunctionIterations.policy(states,prob.u,prob.X,prob.R,prob.F,prob.p,prob.δ,prob.V)
    η_s = reshape(u_opt[2,:],dims...)
    η̄_s =  mapslices(s -> price_per_exposure(El(s[3],s[4],s[5]),prob.p.cf,prob.p.cv),states,dims = 1)
    η̄_s = reshape(η̄_s,dims...)
    u = [1e-10,0.0]
    function insurance_claims(s)
        nodes = prob.X(s,prob.p).nodes
        weights = prob.X(s,prob.p).weights
        weights = weights[nodes[5,:].==1.0] 
        weights = weights ./ sum(weights)
        nodes = nodes[:,nodes[5,:].==1.0] 
        claim = mapslices(x -> index_insurance(index(prob.F(s,u,x,prob.p)[3:5],[0.0]),strike(s[3:5],prob.p),1,1), nodes, dims = 1)
        sum(claim' .* weights )
    end 
    claims =  mapslices(s -> insurance_claims(s),states,dims = 1)
    claims = reshape(claims,dims...)

    return η_s, η̄_s, claims, state_grid
end 


function certianty_equivelent(prob,state)
    inverse_utility(prob.V(state)*(1-prob.δ),prob.p.γ)
end 

function value_diffs(params)
    #params = "no_RA_short"
    Nκ = "20"; Nη = "2"; Ns_grid = "24"; Nb_grid = "24"; NΔb_grid = "18"; Nquad = "50"
    include("../parameters/$params.jl")

    ###################################
    ####  Define random variables #####
    ###################################
    X, Xmc = init_ranom_variables(p; Nquad = parse(Int,Nquad))

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


    # savings 
    smax = 2.5; smin = -4.0; Ns_grid = parse(Int,Ns_grid)
    s_grid = smin:((smax-smin)/(Ns_grid-1)):smax
    plt = Plots.scatter(exp.(s_grid),s_grid, label = "Savings", xlabel = "Savings", ylabel = "log savings")

    state_grid = [s_grid,b_grid,Δb_grid]
    ###############################
    ####  Define action grid  #####
    ###############################
    κmin = 0.01; κmax = 0.99; Nκ = parse(Int,Nκ)
    κ_grid = κmin:((κmax-κmin)/(Nκ-1)):κmax
    ηmin = 0.0; ηmax=1e-10; Nη = parse(Int,Nη)
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


    #######################################
    ####  Strike and index functions  #####
    #######################################
    function strike(s,p)
        return 0.5*(p.b_target +p.b_limit) 
    end

    function index(s,X)
        return exp(s[1])+X[1]
    end 

    ##################################
    ####  Set up dynamic program  ####
    ##################################
    δ = 0.95
    prob, El = init(action_grid,state_grid,δ,p,X,index,strike)
    prob_insurance, El = init(action_grid,state_grid,δ,p,X,index,strike)

    ##################################
    ####  load in solutions       ####
    ##################################
    load_solution(prob,"results/no_insurance_$params/sol.jld2")
    load_solution(prob_insurance,"results/$params/sol.jld2")
    Random.seed!(123)
    sim_states, m, f, h, y = simulate_trajectory(p,Xmc,2000)
    s0 = vcat([1.0,log(-1*p.s̄)],sim_states[:,rand(500:2000)])

    function average_diff(s,states)
        st = log(s - p.s̄)
        V = mean(mapslices(s -> certianty_equivelent(prob_insurance,vcat([1.0,st],s)), states, dims = 1))
        V_no_I = mean(mapslices(s -> certianty_equivelent(prob,vcat([1.0,st],s)), states, dims = 1))
        return (V - V_no_I)/V_no_I
    end

    st_vals = [0.25, 0.5, 1.0, 2.0, 3.5] .+ p.s̄
    values = broadcast(s -> 100*average_diff(s,sim_states),st_vals )
    diff = DataFrame(savings = st_vals, increase = values)
    CSV.write("results/$params/value_diff.csv", diff)
    
    return nothing
end 




function value_diffs(params1, params2)
    #params = "no_RA_short"
    Nκ = "20"; Nη = "2"; Ns_grid = "24"; Nb_grid = "24"; NΔb_grid = "18"; Nquad = "50"
    include("../parameters/$params1.jl")
    include("../parameters/$params2.jl")

    ###################################
    ####  Define random variables #####
    ###################################
    X, Xmc = init_ranom_variables(p; Nquad = parse(Int,Nquad))

    # define state variables and check the grid 
    ######################################
    ####  Define state varible grids #####
    ######################################
    # check paramters to get min and max values
    bmax = 0.5; bmin = -3.5; Nb_grid = parse(Int,Nb_grid)
    Δbmax = 1.25; Δbmin = -1.0; NΔb_grid = parse(Int,NΔb_grid)
    if occursin("short", params1)
        bmax = 1.0; bmin = -3.0
        Δbmax = 2; Δbmin = -1.25
    elseif occursin("long", params1)
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


    # savings 
    smax = 2.5; smin = -4.0; Ns_grid = parse(Int,Ns_grid)
    s_grid = smin:((smax-smin)/(Ns_grid-1)):smax
    plt = Plots.scatter(exp.(s_grid),s_grid, label = "Savings", xlabel = "Savings", ylabel = "log savings")

    state_grid = [s_grid,b_grid,Δb_grid]
    ###############################
    ####  Define action grid  #####
    ###############################
    κmin = 0.01; κmax = 0.99; Nκ = parse(Int,Nκ)
    κ_grid = κmin:((κmax-κmin)/(Nκ-1)):κmax
    ηmin = 0.0; ηmax=1e-10; Nη = parse(Int,Nη)
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


    #######################################
    ####  Strike and index functions  #####
    #######################################
    function strike(s,p)
        return 0.5*(p.b_target +p.b_limit) 
    end

    function index(s,X)
        return exp(s[1])+X[1]
    end 

    ##################################
    ####  Set up dynamic program  ####
    ##################################
    δ = 0.95
    prob, El = init(action_grid,state_grid,δ,p,X,index,strike)
    prob_insurance, El = init(action_grid,state_grid,δ,p,X,index,strike)

    ##################################
    ####  load in solutions       ####
    ##################################
    load_solution(prob,"results/no_insurance_$params2/sol.jld2")
    load_solution(prob_insurance,"results/$params1/sol.jld2")
    Random.seed!(123)
    sim_states, m, f, h, y = simulate_trajectory(p,Xmc,2000)
    s0 = vcat([1.0,log(-1*p.s̄)],sim_states[:,rand(500:2000)])

    function average_diff(s,states)
        st = log(s - p.s̄)
        V = mean(mapslices(s -> certianty_equivelent(prob_insurance,vcat([1.0,st],s)), states, dims = 1))
        V_no_I = mean(mapslices(s -> certianty_equivelent(prob,vcat([1.0,st],s)), states, dims = 1))
        return (V - V_no_I)/V_no_I
    end

    st_vals = [0.25, 0.5, 1.0, 2.0, 3.5] .+ p.s̄
    values = broadcast(s -> 100*average_diff(s,sim_states),st_vals )
    diff = DataFrame(savings = st_vals, increase = values)
    CSV.write("results/$params1/value_diff+$params2.csv", diff)
    
    return nothing
end 