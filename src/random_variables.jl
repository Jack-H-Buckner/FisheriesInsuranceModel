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
        Ib,st,bt,Δbt,mt=s
        if Ib >1.0 # If bankrupt don't both calculating over all nodes
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