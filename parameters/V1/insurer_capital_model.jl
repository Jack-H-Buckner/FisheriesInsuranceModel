#####################################################################
### Parameters for the insurer capital-cost model
### (see src/insurer_capital_model.jl). Stock parameters are defined
### separately in the RA_*.jl / no_RA_*.jl files.
#####################################################################

# Base portfolio (all lines excluding the fisheries line)
n      = 20       # number of independent lines of business
Q_base = 1.0      # expected loss of each line
CV     = 0.5      # coefficient of variation of each line's loss

# Capital requirement and cost of capital
α = 0.99          # SST Expected Shortfall level
τ = 0.06          # cost of holding surplus capital (SST market value margin)

# Market structure
K = 0.025         # fixed participation cost per insurer

# Variable insurance costs (taxes, commissions), per unit of exposure.
# 1/3 gives a 50% loading factor (see premiums() in src/insurnace_contracts.jl)
c_v = 0.05

# Bundle for model.jl's `init`, which appends these to the parameter set `p`
insurer_pars = (K = K, n = n, Q_base = Q_base, CV = CV, α = α, τ = τ, c_v = c_v)
