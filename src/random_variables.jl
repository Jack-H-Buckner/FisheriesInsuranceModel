#####################################################################
#####################################################################
### Define random variables used to quantify underitnaty in the 
### fisheries insurance demand model.
### The full model is defined in model.jl. I have chosen to split the 
### code in this way to keep individual files from being too long.  
### Jack H. Buckner, Oregon State University, 2025
#####################################################################
#####################################################################

using ValueFunctionIterations

# Returns a random variable object with the firts three rows have the values of the
# basis risk added to the index, prices for catch, and cost of effort.  The next two rows have
# the recruitment deviations and the mortlatiy rates. 
function init_ranom_variables(p;Nquad = 10)
    # recruitment deviation are normal random variables
    X1 = GaussHermiteRandomVariable(Nquad,[0.0],[p.σ_R^2;;])
    X2 = MarkovChain(p.T) # mortlatiy rates and the index
    X3 = RandomVariable([p.μZ; p.μ_P; p.μ_Ce;;], [1.0,1.0,1.0]) # basis risk, prices and costs are fixed at their means
    Xmat = product(product(X3,X1),X2)

    # Define a function so that the full set of quadrature nodes
    # is only used when the firm is not bankrupt to save on compute. 
    function Xf(s,p)
        disc = s[1] # discrete state: bankrupt states have an even index
        if iseven(round(Int,disc)) # If bankrupt don't both calculating over all nodes
            return RandomVariable([0.0;;], [1.0])
        end 
        return Xmat 
    end 

    function sampler()
            
        return [p.μZ, p.μ_P, p.μ_Ce, rand(Distributions.Normal(0,p.σ_R))]
    end 
    Xmc = MCRandomVariable(sampler, X2, 1)
    return ValueFunctionIterations.RandomVariableFunction(Xf), Xmc
end

# Random variables for model_with_relief.jl. Same as `init_ranom_variables`, plus one
# Bernoulli draw used by `relief` in disaster_relief.jl:
#   row 6, IRrng: relief payment; paid when `p.pR > IRrng`, so the node 0.0
#                 (probability p.pR) means paid and 1.0 means not paid.
# Disaster declaration after a high mortality event is certain, so it has no draw.
# The two-point distribution replaces a uniform draw without approximation.
# Requires p.pR.
function init_ranom_variables_relief(p;Nquad = 10)
    X1 = GaussHermiteRandomVariable(Nquad,[0.0],[p.σ_R^2;;])
    X2 = MarkovChain(p.T) # mortlatiy rates and the index
    X3 = RandomVariable([p.μZ; p.μ_P; p.μ_Ce;;], [1.0,1.0,1.0])
    XD = RandomVariable([0.0 1.0], [p.pR, 1-p.pR])
    Xmat = product(product(product(X3,X1),X2),XD)

    # only use the full set of quadrature nodes when the firm is not bankrupt
    # (bankrupt states have an even discrete index)
    function Xf(s,p)
        if iseven(round(Int,s[1]))
            return RandomVariable([0.0;;], [1.0])
        end 
        return Xmat 
    end 

    # the sampler supplies rows 1-4; the discrete product (Markov chain and the
    # relief draw) supplies rows 5-6
    function sampler()
        return [p.μZ, p.μ_P, p.μ_Ce, rand(Distributions.Normal(0,p.σ_R))]
    end 
    Xmc = MCRandomVariable(sampler, product(X2,XD), 1)
    return ValueFunctionIterations.RandomVariableFunction(Xf), Xmc
end
