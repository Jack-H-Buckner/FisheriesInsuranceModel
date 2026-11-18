#####################################################################
#####################################################################
### Defines insurance claims payments and premium rates for 
### the fisheries insurance demand model.
### The full model is defined in model.jl. I have chosen to split the 
### code in this way to keep individual files from being too long.  
### Jack H. Buckner, Oregon State University, 2025
#####################################################################
#####################################################################

# The popualtion model (population_model.jl) is needed to update the index of the stock 
# import this file after popualtion_model.jl

# define the index used as the basis of the contract
# depend on the state of the system and a random variable
# that can be used to include basis risk 
function index(s,X)
    return exp(s[1]) + X[1]
end 

# claims paments are proportional to the ammount the
# index of the stock (zt) falls below a threshold (zc)
# The coverage level ϕ determines the proportionality constant
# and is scaled so that the expected indemnity is one when ϕ = 1
# the input El = max(0,zc-zt) enforces this scaling
function index_insurance(zt,zc,ϕ)
    ϕ*max(0,zc-zt)
end 

# The price per unit of exposure equals one plus variable costs like taxes
# and commisions (c_v), and surplus capital costs (c_k). The capital costs are a 
# complex funciton of the stock's abundance and other factors that determine a
# market equilibrim for the insurance preiums
function unit_premiums(s0, strike_, pars, K, n, Q_base, CV, α, τ,c_v; Nmc = 20_000, seed = 1)
    c_k = capital_costs(s0, strike_, pars, K, n, Q_base, CV, α, τ; Nmc = Nmc, seed = seed)
    return 1 + c_v + c_k
end 

# Coverage ϕ bought with premium spending η at price η̄ per unit of coverage.
# The nugget keeps the division finite where the expected payout is ~0.
function coverage(η̄,η)
    nugget = 1e-12
    η/max(η̄,nugget)
end 

# set strike on the slope of the HCR 
# payout begin when the stoke is on the slope of the harvest 
# control rule but before the fishery is fully closed. 
function strike(s,p)
    return 0.5*(p.b_limit+p.b_target)
end 


# Premium per unit of coverage (ϕ = 1) in a single state s = [log stock, Δ log stock, mortality state]:
# the unit premium (1 + c_v + c_k) from the insurer capital model times the
# expected payout El = E[max(0, zc - zt)] computed with the quadrature random variables X.
# Returns (η̄, El); the loading is η̄/El = 1 + c_v + c_k.
function premium_per_coverage(s,X,pars,index,strike; Nmc = 20_000, seed = 1)
    Xt = X(vcat(ones(2),s),pars)
    zt = mapslices(x -> index(update_stock_control_rule(s,x[4:5],pars)[1],x)[1], Xt.nodes, dims = 1)
    zt = reshape(zt,length(Xt.weights))
    zc = strike(s,pars)
    η1 = unit_premiums(s, zc, pars, pars.K, round(Int,pars.n), pars.Q_base, pars.CV,
                       pars.α, pars.τ, pars.c_v; Nmc = Nmc, seed = seed)
    El = sum(index_insurance.(zt,zc,1).*Xt.weights)
    return η1 * El, El
end

# Calculate the premium per unit of coverage over a grid of stock states and
# fit an interpolation function to the grid.
# The premium for coverage ϕ = 1 (payout max(0, zc - zt)) is the unit premium
# (1 + c_v + c_k) times the expected payout El = E[max(0, zc - zt)].
# The insurer parameters (K, n, Q_base, CV, α, τ, c_v) are read from `pars`
# (see `add_insurer_parameters`). Returns an interpolation object called as
# f(log stock, Δ log stock, mortality state).
function iterpolate_premiums(bt_grid,Δb_grid,mt_grid,X,pars,index,strike; Nmc = 20_000, seed = 1)
    η   = zeros(length(bt_grid),length(Δb_grid),length(mt_grid))
    Nbt = length(bt_grid); NΔb = length(Δb_grid); Nmt = length(mt_grid)
    for (i,j,k) in Iterators.product(1:Nbt,1:NΔb,1:Nmt)
        s = [bt_grid[i],Δb_grid[j],mt_grid[k]]
        η[i,j,k] = premium_per_coverage(s,X,pars,index,strike; Nmc = Nmc, seed = seed)[1]
    end
    interp_η  = extrapolate(Interpolations.scale(interpolate(η, BSpline(Linear())),
                                                 bt_grid,Δb_grid,mt_grid), Flat())
    return interp_η
end

# Append the insurer capital-model parameters to the parameter set so they travel
# with `p` (and `prob.p`). `insurer` is a NamedTuple with fields
# K, n, Q_base, CV, α, τ, c_v (see parameters/insurer_capital_model.jl).
function add_insurer_parameters(p, insurer)
    return ComponentArray(p; pairs(insurer)...)
end
