#####################################################################
#####################################################################
### Defines the population dynamics subcomponent of the fisheires
### insurance demand model. 
### The full model is defined in model.jl. I have chosen to split the 
### code in this way to keep individual files from being too long.  
### Jack H. Buckner, Oregon State University, 2025
#####################################################################
#####################################################################

# Fishing effort/mortaltiy is determined by a piece wise function of the 
#stock bt called a harvest control rule
function control_rule(bt,b_target,f_target,b_limit)
    if bt > b_target
        return f_target
    else
        return max(0,f_target*(bt - b_limit) / (b_target - b_limit))
    end
end 

# The Baranov equaion relates fishing mortaltiy to catch give 
# the stock biomass bt and natrual mortaltiy mt 
function baranov(bt,mt,ft)
    return ft/(ft+mt)*(1-exp(-mt-ft))*bt
end

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

# describes chagnes in the abundance of the stock over one time period
# given the stock's current state defined by its biomass bt, lagged biomass bt_1,
#  and lagged natural mortaltiy mt_1. The new stock also depeds on the relization
# of a random variable X which describes both recruitment shocks and natural mortality. 
function update_stock(s,X,p)
    Xt,ΔXt,Mt_1=s # unpack states
    bt = exp(Xt); bt_1 = exp(Xt-ΔXt)
    rt,M_rng = X # unpack random variable
    Mt = sample_markov_chain(round(Int,Mt_1),p.T,M_rng)
    mt = environment_mortaltiy(Mt,p.m̄,p.m̲)#,p.m̲)
    mt_1 = environment_mortaltiy(Mt_1,p.m̄,p.m̲)#,p.m̲)
    ft = control_rule(bt,p.b_target,p.f_target,p.b_limit) # lagged fishing mortlatiy 
    ft_1 = control_rule(bt_1,p.b_target,p.f_target,p.b_limit) # lagged fishing mortlatiy     
    FR = beverton_holt(bt_1,p.b0,p.r0,p.z,rt) # recruitmemt
    bt1 = deriso_schnute(bt,bt_1,mt,mt_1,ft,ft_1,p.g,FR) # pop model 
    ht = baranov(bt,mt,ft)
    if bt1 < 0 
        # The popualtion cannot go below zero, but it is possible
        # for the model to predit this if the current biomass is inconsistent with 
        # lagged abundnce and mortlatiy rates. 
        return [0.0, bt1-bt, Mt], ht, ft 
    end 
    return [log(bt1), log(bt1)-log(bt), Mt], ht, ft
end


