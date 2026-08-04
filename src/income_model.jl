#####################################################################
#####################################################################
### Defines the income and savings model for the fisheries insurance
### demand model.
### The full model is defined in model.jl. I have chosen to split the 
### code in this way to keep individual files from being too long.  
### Jack H. Buckner, Oregon State University, 2025
#####################################################################
#####################################################################

# Income from fishing is the revenue from harvest pt x ht - cost of effort ft x c_et
function income(ht, ft, pt, c_et)
    return pt * ht - c_et * ft 
end 

# borrowing cannot exceed s̄
function budge_constraint(st,s̄)
    return st-s̄
end 

# savings equal current saving minus premiums and consumption
function savings(st,κt, ηt)
    return st - κt - ηt
end 

# update savings - subtract premiums ηt and consumption κt
# add income from fishing yt, insurance calims wt, and alt income / costs y0
# savings are increased by a factor k to account for interest.
function update_savings(st,κt, ηt, yt, wt, y0, k)
    s = savings(st,κt, ηt)
    return (1+k)*s + yt + wt + y0
end 

# if saving sfall below the borrowing constaint at the end of the period the
# the firm become bankrupt. This can occur if saving from the prior period,
# income from fishign and insurance claims cannot cover interest paymetns on 
# debt (ks_t when s_t < 0) and fixed costs if y0 < 0 .
function bankrupt(st,s̄, Ib)
    return (st-s̄  <= 0) | (Ib > 1.0)
end 

# utility from consumption is given by a CRRA utility
# bankruptcy is characterized by a constant, but lower utility.
function utility(κt,I_bankrupt,κ_bankrupt,γ)
    if I_bankrupt > 1.0 
        return (κ_bankrupt^(1-γ)-1)/(1-γ)+1
    end
    return (κt^(1-γ)-1)/(1-γ)+1
end 


function inverse_utility(Ut,γ)
    return (Ut*(1-γ)+1)^(1/(1-γ))
end 