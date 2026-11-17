using Roots

function R0(b0,m̄,g)
    return b0*(1-exp(-m̄,))*(1-g*exp(-m̄))
end 

function f_spr(spr,m̄,g)
    a = spr*g*exp(-2*m̄)
    b = -1*spr*(1+g)*exp(-m̄)
    c = spr - 1 + exp(-m̄) + g*exp(-m̄) - g*exp(-2*m̄)
    -log((-b-sqrt(b^2-4*a*c))/(2*a))
end 

function cost_By0(Pt,m̄,by0)
    (Pt / m̄)*(1-exp(-m̄))*by0
end 

# requires population_model.jl
function mean_price(b_target,f_target,y_target,m̄,by0)
    H_target = baranov(b_target,m̄,f_target)
    y_target/(H_target - f_target*(1-exp(-m̄))*by0/ m̄)
end 


function f_b_target(z,r0,b0,b_target,m̄,g)
    a = g*exp(-2*m̄)
    b = -1*(1+g)*exp(-m̄)
    c = 1 - 4*z*r0/((1-z)*b0+(5*z-1)*b_target) 
    -log((-b-sqrt(b^2-4*a*c))/(2*a))
end


# Calibrate the price μ_P so that average fishing income equals `target`, holding the
# cost of effort μ_Ce (and everything else in p) fixed. Average income comes from a
# long simulation (`mean_income`: `burnin` years discarded, then `T` years) that uses
# the same random numbers (`seed`) for every price, so it is a smooth function of
# the price and can be solved with a bracketing root finder.
# requires population_model.jl; returns a copy of p with the calibrated μ_P.
function calibrate_price(p; target = 1.0, T = 2000, burnin = 50, seed = 1, tol = 1e-8)
    resid(P) = begin
        q = copy(p); q.μ_P = P
        mean_income(q; T = T, burnin = burnin, seed = seed) - target
    end
    hi = p.μ_P > 0 ? p.μ_P : 1.0
    while resid(hi) < 0                      # income increases with the price
        hi *= 2
        hi < 1e8 || error("no price reaches average income $target")
    end
    lo = hi/2
    while resid(lo) > 0
        lo /= 2
        lo > 1e-12 || error("average income exceeds $target even at a price near zero")
    end
    P = find_zero(resid, (lo, hi), Bisection(); xatol = tol)
    q = copy(p); q.μ_P = P
    return q
end
