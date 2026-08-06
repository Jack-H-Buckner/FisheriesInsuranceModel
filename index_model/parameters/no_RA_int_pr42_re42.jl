#####################################################################
#####################################################################
### GENERATED FILE - do not edit by hand.
### Regenerate with:
###     julia --project=. index_model/parameters/generate_scenarios.jl
###
### Index insurance scenario: no_RA_int_pr42_re42
###   base       no_RA_intermidiate
###   arm        iso_premium
###            re = pr, so the payout rate is pinned to the event rate and the premium is held at the perfect index level.  Isolates hedge quality from price.
###   precision  0.42
###   recall     0.42
###   cost       inherited from the base parameter set (actuarially fair)
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

# parameters/*.jl call mean_price -> baranov, so the population model has to
# be loaded first.  Use the index_model copy, never src/population_model.jl:
# the latter would overwrite update_stock with the two state version.
include(joinpath(@__DIR__, "..", "population_model.jl"))
include(joinpath(@__DIR__, "..", "..", "parameters", "no_RA_intermidiate.jl"))
include(joinpath(@__DIR__, "..", "index_transitions.jl"))

# p_low and p_stay are inherited from the base parameter set's 2x2 T, and
# cf / cv from parameters/universal_params.jl unless overridden below.
p = set_index_T(p; pr = 0.42, re = 0.42)
