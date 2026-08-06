#####################################################################
#####################################################################
### GENERATED FILE - do not edit by hand.
### Regenerate with:
###     julia --project=. index_model/parameters/generate_scenarios.jl
###
### Index insurance scenario: no_RA_int_pr33_re100_cv50
###   base       no_RA_intermidiate
###   arm        full_coverage
###            re = 1, so every event is covered and declining precision buys that coverage with false positives.  Isolates price from hedge quality.
###   precision  0.333
###   recall     1.0
###   cost       cv = 0.5 (override)
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
p = set_index_T(p; pr = 0.333, re = 1.0)
p.cv = 0.5
