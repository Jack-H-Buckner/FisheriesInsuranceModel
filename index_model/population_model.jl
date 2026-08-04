###################################################
###################################################
## Modified population model for index insurance 
## with basis risk. 
##
## This pulls from the existing model as much as 
## possible updating functions where needed
##
## 
###################################################
###################################################
include("../src/population_model.jl")


# defines the natrual mortlatiy rate given a draw from a Markov chain 
# with values 1, 2, 3 and 4. 
# states | 1 | 2 | 3 | 4 |
# -------|---|---|---|---|
# M      | 1 | 1 | 2 | 2 |
# x      | 1 | 2 | 1 | 2 |
#
# M = 1 typical mortality 
# M = 2 high mortality
# x = 1 no payout 
# x = 2 high payout
function environment_mortaltiy_index(Mt,m̄,m̲)
    
    if (Mt == 1) | (Mt == 2)
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
    mt = environment_mortaltiy_index(Mt,p.m̄,p.m̲)#,p.m̲)
    mt_1 = environment_mortaltiy_index(Mt_1,p.m̄,p.m̲)#,p.m̲)
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

