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
# The payout is decided by the joint (mortality, payout) draw for THIS period,
# not by anything carried in the state.  update_stock draws the same joint
# outcome from the same state and the same uniform M_rng, so the payout and the
# mortality realisation are automatically consistent - that coupling is the
# basis risk the contract is built around.
function index(s,X,p)
    _,_,_,_,Mt_1=s
    _,_,_,_,M_rng = X
    return joint_payout(draw_joint(Mt_1,p,M_rng)) ? 1.0 : 0.0
end

# calims paments are proportional to the ammount the
# index of the stock (zt) falls below a threshold (zc)
# given the  
function index_insurance(zt,ηt,η̄t)
    ηt/η̄t*max(0,zt)
end 

# calculate the price per exposure given the 
# expected losses (El), fixed costs (cf), and variable costs / subsidies (cv)
function price_per_exposure(EL,cf,cv)
    nugget = 1e-12
    (max(nugget,EL)+cf)/(1-cv)
end 


# Expected losses: the probability that next period's joint draw lands in a
# payout state, i.e. the mass on rows 2 and 4 (x = 2) of the column the carried
# mortality state draws from.  This is the actuarially fair premium rate.
#
# Under joint_transition it is the same for every state and equals
# p_index = re*p_low/pr, so the fair price is flat and the variation in coverage
# is demand rather than price.  Written state by state anyway, so a contract
# with a state dependent payout rate would still be priced correctly.
function El(s,p)
    _,_,_,_,Mt_1=s
    sum(p.T[[2,4],joint_column(Mt_1)])
end 