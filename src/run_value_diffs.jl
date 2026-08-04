#####################################################################
### Runs `value_diffs(params)` for every parameter set called (un-
### commented) in notebooks/tests.ipynb. Writes
### results/<params>/value_diff.csv for each.
###
### Run from src/:
###     julia run_value_diffs.jl
#####################################################################

cd(@__DIR__)

include("analysis.jl")

param_sets = [
    "no_RA_short",
    "no_RA_intermidiate",
    "no_RA_long",
    "RA_short",
    "RA_intermidiate",
    "RA_long",
    "RA_intermidiate_higher_severity",
    "RA_long_higher_frequency",
    "RA_long_higher_severity",
    "RA_long_longer_duration",
    "RA_short_higher_severity",
    "no_RA_long_higher_frequency",
    "no_RA_long_higher_severity",
    "no_RA_long_longer_duration",
    "RA_intermidiate_lower_frequency",
    "RA_intermidiate_higher_frequency_2",
    "RA_intermidiate_higher_frequency_3",
    "no_RA_intermidiate_higher_y0",
    "RA_intermidiate_lower_blimit",
    "no_RA_intermidiate_higher_k",
    "no_RA_intermidiate_lower_y0",
    "RA_intermidiate_higher_k",
    "RA_intermidiate_higher_k_2",
    "RA_intermidiate_higher_k_3",
]

for params in param_sets
    println("=== value_diffs($params) ===")
    try
        value_diffs(params)
    catch err
        @warn "value_diffs failed" params exception=(err, catch_backtrace())
    end
end
