#####################################################################
#####################################################################
### This file combines the population model, income model, and 
### insurance contracts to define the fisheries insurance demand model.
### and uses the ValueFunctionIterations.jl package to solve the model.
###
### The full model is defined in model.jl. I have chosen to split the 
### code in this way to keep individual files from being too long.  
### Jack H. Buckner, Oregon State University, 2025
#####################################################################
#####################################################################

using ValueFunctionIterations, Interpolations, StatsBase, Distributions

include("population_model_order_2.jl")
include("income_model.jl")
include("derived_parameters.jl")
include("random_variables.jl")
include("insurnace_contracts.jl")

# Simulate a trajectory from the biological model to help
# deterine a good grid for a given parameter set.
# p: paramters
# X: random variables
# T: length of trajectory
function simulate_trajectory(p,X::ValueFunctionIterations.RandomVariableFunction,T::Int)
    states = zeros(3,T+1)
    states[:,1] = [log(p.b_target),0.0,1.0]
    y = zeros(T); h = zeros(T); f = zeros(T); m = zeros(T)
    for t in 1:T
        Xt = X(vcat([1.0,1.0],states[:,t]),p)
        ind = sample(1:length(Xt.weights),Weights(Xt.weights))
        states[:,t+1] = update_stock(states[:,t],Xt.nodes[4:5,ind],p)[1]
        m[t] = environment_mortaltiy(states[3,t+1],p.m̄,p.m̲)
        f[t] = control_rule(exp(states[1,t+1]),p.b_target,p.f_target,p.b_limit)
        h[t] = baranov(exp(states[1,t+1]),m[t],f[t])
        y[t] = income(h[t], f[t], p.μ_P, p.μ_Ce)
    end 
    return states, m, f, h, y
end 


function simulate_trajectory(p,X::ValueFunctionIterations.MCRandomVariable,T::Int)
    states = zeros(3,T+1)
    states[:,1] = [log(p.b_target),0.0,1.0]
    y = zeros(T); h = zeros(T); f = zeros(T); m = zeros(T)
    for t in 1:T
        X()
        ind = sample(1:length(X.weights),Weights(X.weights))
        states[:,t+1] = update_stock(states[:,t],X.nodes[4:5,ind],p)[1]
        m[t] = environment_mortaltiy(states[3,t+1],p.m̄,p.m̲)
        f[t] = control_rule(exp(states[1,t+1]),p.b_target,p.f_target,p.b_limit)
        h[t] = baranov(exp(states[1,t+1]),m[t],f[t])
        y[t] = income(h[t], f[t], p.μ_P, p.μ_Ce)
    end 
    return states, m, f, h, y
end 

function simulate_trajectory(p,X,T)
    states = zeros(3,T+2)
    states[:,1] = [p.b_target,p.b_target,1.0]
    for t in 1:(T+1)
        states[:,t+1] = update_stock(states[:,t],X[1:2],p)[1]
    end 
    return states
end 


function init(action_grid,state_grid,δ,p,X,index,strike)

    # unpack grid and set grid for mortlatity and bankruptcy 
    κ_grid,η_grid = action_grid
    s_grid,b_grid,Δb_grid = state_grid 
    m_grid = 1.0:2.0; Ib_grid = 1.0:2.0

    # calculate expected losses
    bmax = b_grid[end]+1; bmin = b_grid[1]-1; Nb_grid = 150
    b_grid_EL = bmin:((bmax-bmin)/(Nb_grid-1)):bmax
    Δbmax = Δb_grid[end]+0.5; Δbmin = Δb_grid[1]-0.5; NΔb_grid = 150
    Δb_grid_EL = Δbmin:((Δbmax-Δbmin)/(NΔb_grid-1)):Δbmax
    El = iterpolate_losses(b_grid_EL,Δb_grid_EL,m_grid,X,p,index,strike)

    # set up the action grid 
    k=0; action_grid = zeros(2,length(κ_grid)*length(η_grid))
    for (i,j) in Iterators.product(1:length(κ_grid),1:length(η_grid))
        k +=1
        action_grid[:,k] = vcat(κ_grid[i],η_grid[j]) 
    end
    action_grid=action_grid[:,reshape(mapslices(sum,action_grid,dims = 1).<=1,length(κ_grid)*length(η_grid))]

    # Scale action grid so actions dd up to th budget constraint
    # and if bankrupts don't both searching over actions 
    function u(s,p) 
        Ib,st,bt,Δbt,mt=s
        
        budget_constraint = exp(st)
        if p.y0 > 0.0
            budget_constraint += p.y0
        end
        if (Ib > 1.0 )| (budget_constraint <= 0.0)
            return [0.0; 0.0;;]
        end
        return budget_constraint .* action_grid
    end

    # Define the full state update 
    function F(s,u,X,p)

        Ib,st,Xt,ΔXt,Mt_1=s # Unpack states
        st = exp(st)+p.s̄; bt = exp(Xt); bt_1 = exp(Xt-ΔXt)
        if Ib > 1.0 # If bankrupt stay bankrupt skip remaining computations
            return [2.0, 0.0, 0.0, 0.0, 1.0]
        end

        if !(Mt_1 in [1.0,2.0]) # Must be either 1 or 2
            throw("Mt_1 must be one or two")
        end 
        # Unpack actions
        budget_constraint = st - p.s̄
        if p.y0 > 0.0
            budget_constraint += p.y0
        end
        κt,ηt=u 
        if (κt <=0) | (ηt < 0) # Check consumption is positive
            throw("consmption and premiums must be positive")
        elseif (κt+ηt) > budget_constraint
            throw("κt+ηt must be less than st")
        end 
        Zeps,pt,c_et,Rt,mrng =  X
        # calcualte premius 
        zc = strike(s[3:5],p)
        η̄ = price_per_exposure(El(Xt,ΔXt,Mt_1),p.cf,p.cv)
        # update population
        s_next,ht,ft = update_stock(s[3:5],X[4:5],p)
        yt = income(ht, ft, pt, c_et)

        # calcualte insurance 
        wt = index_insurance(index(s_next,X),zc,ηt,η̄)

        # update savings
        st1 = update_savings(st,κt, ηt, yt, wt, p.y0, p.k)

        # if budget constaint is violated go bankrupt 
        if bankrupt(st1,p.s̄,Ib)
            return vcat([2.0,-6.0],s_next)
        end 
        return vcat([1.0,log(st1-p.s̄)],s_next)  
    end 

    function R(s,u,X,p)
        Ib=s[1] # get bankruptcy
        κt=u[1] # get consumption
        utility(κt,Ib,p.κ_bankrupt,p.γ)
    end 

    V = ValueFunction(2,[[s_grid, b_grid, Δb_grid, 1.0:2.0],[0:1,0:1,0:1,1:2]])
    prob = DynamicProgram(V, R, F, p, u,  X, δ; solve = false)

    return prob, El

end 


