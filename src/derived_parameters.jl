
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