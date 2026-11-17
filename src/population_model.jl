#####################################################################
#####################################################################
### Defines the population dynamics subcomponent of the fisheires
### insurance demand model. 
### The full model is defined in model.jl. I have chosen to split the 
### code in this way to keep individual files from being too long.  
### Jack H. Buckner, Oregon State University, 2025
#####################################################################
#####################################################################

# Harvest and fishing mortality (the harvest control rule `control_rule`, the 
# Baranov catch equation `baranov`, and profit maximizing effort `fishing_mortality`)
# are defined in harvest_effort_model.jl so they are shared with that module.
include("harvest_effort_model.jl")
include("income_model.jl") # `income`, used by simulate_trajectory
using Random, ValueFunctionIterations # `sample_markov_chain` is used by update_stock

# The Deriso-Schnute model simulated an age structured fish popuatlion by tracking
# the stock biomass bt and lagged biomass bt_1. The model also tracks the fishing 
# mortality ft, lagged fishing mortality ft_1, as well as current and lagged natural 
# mortaltiy mt and mt_1
function deriso_schnute(bt,bt_1,mt,mt_1,ft,ft_1,g,FR)
    return (1+g)*bt*exp(-mt-ft) - g*bt_1*exp(-mt-ft-mt_1-ft_1) + FR
end 

# The beverton hold model describes density dependet recruitment. 
# the paramter z describes the effect of reducing the stock on future recruitmnet. 
# the paramter bt_τ is the stock biomass lagged by τ periods. b0 and r0 are the 
# unfied biomass and recruitment and rt is the random effect on recruitment. 
function beverton_holt(bt_τ,b0,r0,z,rt)
    return 4*z*r0*bt_τ*exp(rt-0.5*p.σ_R^2)/((1-z)*b0+(5*z-1)*bt_τ)
end 

# defines the natrual mortlatiy rate given a draw from a Markov chain 
# with values 1 or 2. 
function environment_mortaltiy(xt,m̄,m̲)
    if xt == 1
        return m̄
    end
    return m̲
end

# Fishing mortality rules used by `update_stock`. Each takes the biomass bt, natural
# mortality mt and the parameter set p, and returns fishing mortality.
# - profit_max_effort: the fleet maximizes profit at the mean price p.μ_P and effort
#   cost p.μ_Ce, but cannot exceed the harvest allowed by the control rule.
# - control_rule_effort: the fleet fishes exactly at the control rule.
profit_max_effort(bt,mt,p) = fishing_mortality(bt,mt,p,p.μ_P,p.μ_Ce)
control_rule_effort(bt,mt,p) = control_rule(bt,p.b_target,p.f_target,p.b_limit)

# one year update of the stock given the effort rule, the current (new) natural
# mortality mt with its state Mt, and lagged mortality state Mt_1. Shared by the 
# `update_stock` methods below.
function stock_step(s,X,p,effort_rule,mt,Mt)
    Xt,ΔXt,Mt_1=s # unpack states
    bt = exp(Xt); bt_1 = exp(Xt-ΔXt)
    rt = X[1] # recruitment shock
    mt_1 = environment_mortaltiy(Mt_1,p.m̄,p.m̲)#,p.m̲)
    ft = effort_rule(bt,mt,p) # fishing mortlatiy 
    ft_1 = effort_rule(bt_1,mt_1,p) # lagged fishing mortlatiy     
    FR = beverton_holt(bt_1,p.b0,p.r0,p.z,rt) # recruitmemt
    bt1 = deriso_schnute(bt,bt_1,mt,mt_1,ft,ft_1,p.g,FR) # pop model 
    ht = baranov(bt,mt,ft)
    if bt1 < 0 
        # The popualtion cannot go below zero, but it is possible
        # for the model to predit this if the current biomass is inconsistent with 
        # lagged abundance and mortlatiy rates. 
        return [0.0, bt1-bt, Mt], ht, ft 
    end 
    return [log(bt1), log(bt1)-log(bt), Mt], ht, ft
end

# describes changes in the abundance of the stock over one time period
# given the stock's current state defined by its biomass bt, lagged biomass bt_1,
#  and lagged natural mortaltiy mt_1. The new stock also depeds on the relization
# of a random variable X which describes both recruitment shocks and natural mortality. 
# `effort_rule` sets fishing mortality; the default is profit maximization subject to
# the harvest control rule. See `update_stock_control_rule` for the control rule only.
function update_stock(s,X,p; effort_rule = profit_max_effort)
    Mt = sample_markov_chain(round(Int,s[3]),p.T,X[2])
    mt = environment_mortaltiy(Mt,p.m̄,p.m̲)
    return stock_step(s,X,p,effort_rule,mt,Mt)
end

# update the stock with a known mortaltiy rate 
function update_stock(s,X,p,I; effort_rule = profit_max_effort)
    mt = [p.m̄,p.m̲][I+1] # select moraltiayt rate base don event I
    return stock_step(s,X,p,effort_rule,mt,I+1)
end

# the same updates with fishing mortality set by the harvest control rule only
update_stock_control_rule(s,X,p) = update_stock(s,X,p; effort_rule = control_rule_effort)
update_stock_control_rule(s,X,p,I) = update_stock(s,X,p,I; effort_rule = control_rule_effort)


"""
    simulate_trajectory(p, T; seed = nothing) -> (states, m, f, h, y)

Simulate `T` years of the stock under the parameter set `p`, starting from the
target biomass, with fishing mortality set by `profit_max_effort`. Mirrors the
`simulate_trajectory(p, X, T)` methods in model.jl but draws its own random numbers
(a recruitment shock `ε_R ~ Normal(0, p.σ_R)` and a uniform draw for the mortality
Markov chain), so it needs no random variable object and can be used from files
that are loaded before model.jl, e.g. for calibrating parameters.

`seed` fixes the random numbers (common random numbers across calls and across
values of the parameters being tuned); `nothing` uses the global RNG.

# Returns
- `states`: 3 x (T+1) matrix of (log stock, change in log stock, mortality state)
- `m`, `f`, `h`, `y`: natural mortality, fishing mortality, harvest and fishing income
  (price x harvest - cost x effort), for years 1..T, evaluated at the new stock.
"""
function simulate_trajectory(p,T::Int; seed = nothing)
    rng = seed === nothing ? Random.default_rng() : Random.Xoshiro(seed)
    states = zeros(3,T+1)
    states[:,1] = [log(p.b_target),0.0,1.0]
    y = zeros(T); h = zeros(T); f = zeros(T); m = zeros(T)
    for t in 1:T
        X = [p.σ_R*randn(rng), rand(rng)]
        states[:,t+1] = update_stock(states[:,t],X,p)[1]
        m[t] = environment_mortaltiy(states[3,t+1],p.m̄,p.m̲)
        f[t] = profit_max_effort(exp(states[1,t+1]),m[t],p)
        h[t] = baranov(exp(states[1,t+1]),m[t],f[t])
        y[t] = income(h[t], f[t], p.μ_P, p.μ_Ce)
    end 
    return states, m, f, h, y
end

"""
    mean_income(p; T = 2000, burnin = 50, seed = 1) -> Float64

Average fishing income `y` over `T` years after a burn in of `burnin` years, from
`simulate_trajectory` with a fixed seed.
"""
function mean_income(p; T = 2000, burnin = 50, seed = 1)
    y = simulate_trajectory(p,burnin+T; seed = seed)[5]
    return sum(y[(burnin+1):end])/T
end
