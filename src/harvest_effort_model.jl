#####################################################################
#####################################################################
### Defines harvest effort (fishing mortality) as the outcome of profit
### maximization subject to the harvest control rule, as an alternative to
### the control-rule-only policy in population_model.jl (`control_rule`).
### The fleet picks the effort that maximizes profit p*h(f) - c*f, but is
### not allowed to harvest more than the control rule permits.
### This file is self-contained and is included by population_model.jl, which
### uses `control_rule`, `baranov` (the catch equation) and `fishing_mortality`
### from here so the two modules cannot drift apart.
### Jack H. Buckner, Oregon State University, 2025
#####################################################################
#####################################################################

using Roots

# Fishing effort/mortaltiy under the harvest control rule: a piece wise function
# of the stock bt. Zero below b_limit, rising linearly to f_target at b_target.
function control_rule(bt,b_target,f_target,b_limit)
    if bt > b_target
        return f_target
    else
        return max(0,f_target*(bt - b_limit) / (b_target - b_limit))
    end
end 

# Harvest from effort f (fishing mortality) at natural mortality m and biomass b.
# This is the Baranov catch equation, written with z = f + m.
function harvest(f,m,b)
    z = f+m
    return f*b*(1-exp(-z))/z
end

# The Baranov equaion relates fishing mortaltiy to catch given the stock
# biomass bt and natrual mortaltiy mt (same as `harvest`, different argument order).
baranov(bt,mt,ft) = harvest(ft,mt,bt)

# dH/df for H = f*b*(1-exp(-z))/z, z = f+m
function dharvest_df(f,m,b)
    z = f+m
    n = b*(m + exp(-z)*(f^2+f*m-m))
    return n/z^2
end

# profit from price p, unit effort cost c, effort f, natural mortality m, biomass b
function profit(p,c,f,m,b)
    p*harvest(f,m,b)-c*f
end

# first-order condition of profit in f: zero at the unconstrained optimum,
# positive when marginal profit is negative
function profit_foc(p,c,f,m,b)
    c-p*dharvest_df(f,m,b)
end

# harvest allowed by the control rule at biomass b
function harvest_constaint(b,b_limit,b_target,f_target,m)
    return harvest(control_rule(b,b_target,f_target,b_limit),m,b)
end 

"""
    fisher_response(p, c, m, b; b_limit, b_target, f_target, f_max = 100.0)

Profit-maximizing effort from the first-order condition `profit_foc = 0`,
capped by the harvest control rule `harvest_constaint`.
Returns `(harvest, effort, profit)`.

# Arguments
- `p`, `c` : price per unit harvest and cost per unit effort.
- `m`      : natural mortality rate.
- `b`      : stock biomass.
# Keywords
- `b_limit`, `b_target`, `f_target` : harvest control rule parameters.
- `f_max`  : upper bound on effort for the root search.
"""
function fisher_response(p, c, m, b; b_limit=0.0, b_target=1.0, f_target=1.0, f_max=100.0)
    # 1. unconstrained optimum: root of the FOC in f
    foc(f) = profit_foc(p,c,f,m,b)
    if foc(1e-9) >= 0          # marginal profit negative even at zero effort
        f_opt = 0.0
    elseif foc(f_max) <= 0     # no interior root within bounds
        f_opt = f_max
    else
        f_opt = find_zero(foc, (1e-9, f_max), Bisection())
    end
    h_opt = harvest(f_opt,m,b)

    # 2. check against the harvest constraint
    h_con = harvest_constaint(b,b_limit,b_target,f_target,m)
    if h_opt <= h_con
        h, f = h_opt, f_opt
    else
        h = h_con
        # effort that yields the constrained harvest (harvest increasing in f)
        f = h_con <= 0 ? 0.0 : find_zero(x -> harvest(x,m,b) - h_con, (1e-9, f_opt), Bisection())
    end
    return (harvest=h, effort=f, profit=p*h - c*f)
end

"""
    fishing_mortality(bt, mt, p, price, cost) -> Float64

Fishing mortality chosen by the fleet at biomass `bt` and natural mortality `mt`:
the profit-maximizing effort, capped by the control rule in the parameter set `p`
(fields `b_limit`, `b_target`, `f_target`). A drop-in alternative to
`control_rule(bt, p.b_target, p.f_target, p.b_limit)`.
"""
function fishing_mortality(bt, mt, p, price, cost)
    return fisher_response(price, cost, mt, bt; b_limit = p.b_limit,
                           b_target = p.b_target, f_target = p.f_target).effort
end
