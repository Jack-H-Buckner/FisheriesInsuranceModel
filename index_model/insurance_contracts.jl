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
function index(s,X,p)
    _,_,_,_,Mt_1=s 
    _,_,_,_,M_rng = X 
    M_t = sample_markov_chain(round(Int,Mt_1),p.T,M_rng)
    if (M_t == 2) | (M_t == 4)
        return 1.0
    end
    return 0
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


#Calculate the expected losses over a grid of stock states
# and fits an interpolation function to the grid.
function El(s,p)
    #print(s)
    _,_,_,_,Mt_1=s 
    sum(p.T[[2,4],Int(Mt_1)])
end 