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

# calims paments are proportional to the ammount the
# index of the stock (zt) falls below a threshold (zc)
# given the  
function index_insurance(zt,zc,ηt,η̄t)
    ηt/η̄t*max(0,zc-zt)
end 

# calculate the price per exposure given the 
# expected losses (El), fixed costs (cf), and variable costs / subsidies (cv)
function price_per_exposure(EL,cf,cv)
    nugget = 1e-12
    (max(nugget,EL)+cf)/(1-cv)
end 

# set strike on the slope of the HCR 
# payout begin when the stoke is on the slope of the harvest 
# control rule but before the fishery is fully closed. 
function strike(s,p)
    return 0.5*(p.b_limit+p.b_target)
end 


# Calculate the expected losses over a grid of stock states
# and fits an interpolation function to the grid.

function iterpolate_losses(bt_grid,Δb_grid,mt_grid,X,p,index,strike)
    El = zeros(length(bt_grid),length(Δb_grid),length(mt_grid))
    Nbt = length(bt_grid); NΔb = length(Δb_grid); Nmt = length(mt_grid)
    for (i,j,k) in Iterators.product(1:Nbt,1:NΔb,1:Nmt)
        bt = bt_grid[i]; Δb = Δb_grid[j]; mt = mt_grid[k]
        s = [bt,Δb,mt]
        Xt = X(vcat(ones(2),s),p) 
        zt = mapslices(x -> index(update_stock(s,x[4:5],p)[1],x)[1], Xt.nodes, dims = 1) 
        zt = reshape(zt,length(Xt.weights))
        zc = strike(s,p)
        El[i,j,k] = sum(index_insurance.(zt,zc,1,1).*Xt.weights)
    end
    interp=interpolate(El,BSpline(Cubic(Line(OnGrid()))))
    interp=Interpolations.scale(interp,bt_grid,Δb_grid,mt_grid)
    interp= extrapolate(interp, Flat())
    return interp
end 