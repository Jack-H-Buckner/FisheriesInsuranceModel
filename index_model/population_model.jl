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


# The contract's transition matrix T is over the JOINT state
#
#   j      | 1 | 2 | 3 | 4 |
#   -------|---|---|---|---|
#   M      | 1 | 1 | 2 | 2 |
#   x      | 1 | 2 | 1 | 2 |
#
#   M = 1 typical mortality      x = 1 no payout
#   M = 2 high mortality         x = 2 payout
#
# but the model only ever CARRIES M as a state variable.  The payout indicator
# x is drawn fresh each period jointly with M, is paid into savings
# immediately, and (because joint_transition gives T_x identical columns) says
# nothing about the future, so carrying it would double the state space for no
# information.  See draw_joint / joint_mortality / joint_payout below.
#
# If a serially correlated payout signal is ever wanted - a T_x with distinct
# columns - x becomes informative and would have to be carried again.

# Column of the joint T to draw from, given the carried mortality state
# M in {1,2}.  Columns 2M-1 and 2M are identical, so either represents M.
joint_column(Mt_1) = 2*round(Int,Mt_1) - 1

# Split a joint draw j in 1:4 back into its mortality state and payout flag.
joint_mortality(j) = (j + 1) ÷ 2    # 1,2 -> 1 ;  3,4 -> 2
joint_payout(j)    = iseven(j)      # x = 2

# Draw the joint (mortality, payout) outcome for one period.  update_stock and
# index both call this with the same state and the same uniform draw, so the
# payout and the mortality realisation stay coupled - which is exactly the
# basis risk the contract is built around.
draw_joint(Mt_1,p,M_rng) = sample_markov_chain(joint_column(Mt_1),p.T,M_rng)

# Natural mortality rate for a carried state M in {1,2}.  This is the two state
# environment_mortaltiy from src/population_model.jl; kept as a named wrapper so
# the index model reads explicitly rather than relying on the base model's name.
environment_mortaltiy_index(Mt,m̄,m̲) = environment_mortaltiy(Mt,m̄,m̲)

# describes chagnes in the abundance of the stock over one time period
# given the stock's current state defined by its biomass bt, lagged biomass bt_1,
#  and lagged natural mortaltiy mt_1. The new stock also depeds on the relization
# of a random variable X which describes both recruitment shocks and natural mortality. 
function update_stock(s,X,p)
    Xt,ΔXt,Mt_1=s # unpack states
    bt = exp(Xt); bt_1 = exp(Xt-ΔXt)
    rt,M_rng = X # unpack random variable
    # Draw the joint (mortality, payout) outcome, then carry only the mortality
    # component forward as the state.
    Mt = joint_mortality(draw_joint(Mt_1,p,M_rng))
    mt = environment_mortaltiy_index(Mt,p.m̄,p.m̲)
    mt_1 = environment_mortaltiy_index(Mt_1,p.m̄,p.m̲)
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

