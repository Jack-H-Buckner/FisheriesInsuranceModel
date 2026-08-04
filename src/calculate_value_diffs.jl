println(VERSION)
println("Current environment: ", Base.active_project())
using ValueFunctionIterations
println(ValueFunctionIterations.RandomVariableFunction)
include("../src/model.jl")
include("../src/analysis.jl")

println(ValueFunctionIterations.RandomVariableFunction)
println("running!")
### Change in probabilities (table S3) ###
# no risk aversion 
value_diffs("no_RA_short")
println("set 0")
value_diffs("no_RA_short_higher_frequency")
value_diffs("no_RA_short_higher_severity")
value_diffs("no_RA_short_longer_duration")
println("set 1")

value_diffs("no_RA_intermidiate")
value_diffs("no_RA_intermidiate_lower_frequency")
value_diffs("no_RA_intermidiate_higher_frequency")
value_diffs("no_RA_intermidiate_higher_severity")
value_diffs("no_RA_intermidiate_longer_duration")

value_diffs("no_RA_long")
value_diffs("no_RA_long_higher_frequency")
value_diffs("no_RA_long_higher_severity")
value_diffs("no_RA_long_longer_duration")

# risk aversion 
value_diffs("RA_short")
value_diffs("RA_short_higher_frequency")
value_diffs("RA_short_higher_severity")
value_diffs("RA_short_longer_duration")

value_diffs("RA_intermidiate")
value_diffs("RA_intermidiate_lower_frequency")
value_diffs("RA_intermidiate_higher_frequency")

value_diffs("RA_intermidiate_higher_severity")
value_diffs("RA_intermidiate_longer_duration")

value_diffs("RA_long")
value_diffs("RA_long_higher_frequency")
value_diffs("RA_long_higher_severity")
value_diffs("RA_long_longer_duration")

### change in costs and management (table S4, figure S8.A) ###
value_diffs("no_RA_intermidiate_higher_y0")
value_diffs("no_RA_intermidiate_lower_y0")
value_diffs("no_RA_intermidiate_higher_k")
value_diffs("no_RA_intermidiate_lower_blimit")

value_diffs("RA_intermidiate_lower_y0")
value_diffs("RA_intermidiate_higher_RA")
value_diffs("RA_intermidiate_higher_k")
value_diffs("RA_intermidiate_lower_blimit")



### Frequency of shocks figure 4.B ###
value_diffs("no_RA_intermidiate_higher_frequency_2")
value_diffs("no_RA_intermidiate_higher_frequency_3")

value_diffs("RA_intermidiate_higher_frequency_2")
value_diffs("RA_intermidiate_higher_frequency_3")


### Interest rates changd  (Figure S8.B,C) ###
value_diffs("no_RA_intermidiate_higher_k_2")
value_diffs("no_RA_intermidiate_higher_k_3")

value_diffs("RA_intermidiate_higher_k_2")
value_diffs("RA_intermidiate_higher_k_3")