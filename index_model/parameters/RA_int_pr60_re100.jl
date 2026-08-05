#####################################################################
#####################################################################
### GENERATED FILE - do not edit by hand.
### Regenerate with:
###     julia --project=. index_model/parameters/generate_scenarios.jl
###
### Index insurance scenario: RA_int_pr60_re100
###   base       RA_intermidiate
###   arm        full_coverage
###            re = 1, so every event is covered and declining precision buys that coverage with false positives.  Isolates price from hedge quality.
###   precision  0.6
###   recall     1.0
###   cost       inherited from the base parameter set (actuarially fair)
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

# parameters/*.jl call mean_price -> baranov, so the population model has to
# be loaded first.  Use the index_model copy, never src/population_model.jl:
# the latter would overwrite update_stock with the two state version.
include(joinpath(@__DIR__, "..", "population_model.jl"))
include(joinpath(@__DIR__, "..", "..", "parameters", "RA_intermidiate.jl"))
include(joinpath(@__DIR__, "..", "index_transitions.jl"))

# p_low and p_stay are inherited from the base parameter set's 2x2 T, and
# cf / cv from parameters/universal_params.jl unless overridden below.
p = set_index_T(p; pr = 0.6, re = 1.0)
