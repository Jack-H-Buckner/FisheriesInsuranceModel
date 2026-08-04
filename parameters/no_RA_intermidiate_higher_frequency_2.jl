#####################################################################
#####################################################################
### This file defines a parameter set for a fish population with and
### intermidate somatic growth and natrual mortaltiy rate.  
### Jack H. Buckner, Oregon State University, 2025
#####################################################################
#####################################################################

using ComponentArrays
include("../src/derived_parameters.jl")
include("universal_params.jl")

g = exp(-0.162) # somatic growth rate 
m̄ = 0.2 # base line mortality rate  
r0 = R0(b0,m̄,g) # recruitment without fishing
f_target = f_b_target(z,r0,b0,b_target,m̄,g)
μ_P = mean_price(b_target,f_target,y_target,m̄,by0 )
T = [0.925 0.5; 0.075 0.5]
# Put paramters in a component array
p = ComponentArray(
    (   
        b0 = b0, g = g, m̄ = m̄, # growth paramters 
        m̲ = m̄ - log(low_survival), # higherr mortality rate
        r0 = r0, # recruitment without fishing
        z = z, # stock recruitment steepness
        b_target = b_target, # target reference point 
        b_limit = b_limit, # limit reference point 
        f_target = f_target, # target fishign mortality rate
        T = T, # mortality rate transitions
        σ_R = σ_R, # variance of recruitment deviations
        μZ = 0.0, # mean price for catch 
        μ_P = μ_P , # mean price for catch 
        μ_Ce = cost_By0(μ_P,m̄,by0), # Mean cost of effort 
        s̄ = s̄, # borrowing budget_constraint
        k = k, # inflation adjusted interest rate
        y0 = y0_no_RA, # fixed cost / alterantive income 
        cf = cf, # fixed insuance costs
        cv = cv, # variable insurance costs,
        κ_bankrupt = κ_bankrupt_no_RA, # consumption when bankrupt
        γ = γ_no_RA # Risk aversion 
    )
)



