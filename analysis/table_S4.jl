using DataFrames, CSV

# Table S4: percent change in expected present utility (with vs. without
# insurance) across economic-parameter variants on the intermediate life
# history. Reads the same `results/<params>/value_diff.csv` outputs produced
# by `value_diffs` in src/analysis.jl. Table columns are xt budgets
# {0.25, 0.5, 1.0, 2}; the 3.5 row is dropped.

xt_vals = [0.25, 0.5, 1.0, 2.0]

aversion_blocks = [
    (
        label = "Bankruptcy",
        scenarios = [
            (label = "Base",                            params = "no_RA_intermidiate"),
            (label = "Lower fixed costs (higher y0)",   params = "no_RA_intermidiate_higher_y0"),
            (label = "Higher fixed costs (lower y0)",   params = "no_RA_intermidiate_lower_y0"),
            (label = "Higher interest rate k",          params = "no_RA_intermidiate_higher_k"),
        ],
    ),
    (
        label = "Risk-averse",
        scenarios = [
            (label = "Base",                            params = "RA_intermidiate"),
            (label = "Higher fixed costs (lower y0)",   params = "RA_intermidiate_lower_y0"),
            (label = "Higher risk-aversion a",          params = "RA_intermidiate_higher_RA"),
            (label = "Higher interest rate k",          params = "RA_intermidiate_higher_k"),
            (label = "Lower B_limit",                   params = "RA_intermidiate_lower_blimit"),
        ],
    ),
]

table_S4 = DataFrame(
    aversion     = String[],
    scenario     = String[],
    params       = String[],
    xt           = Float64[],
    pct_increase = Float64[],
)

for av in aversion_blocks, sc in av.scenarios
    path = "results/$(sc.params)/value_diff.csv"
    if !isfile(path)
        @warn "Missing value_diff.csv" sc.params path
        for xt in xt_vals
            push!(table_S4, (av.label, sc.label, sc.params, xt, NaN))
        end
        continue
    end
    df = CSV.read(path, DataFrame)
    # value_diffs writes savings rows in fixed order [0.25,0.5,1.0,2.0,3.5] .+ s̄;
    # take the first four positionally.
    for (i, xt) in enumerate(xt_vals)
        push!(table_S4, (av.label, sc.label, sc.params, xt, df.increase[i]))
    end
end

CSV.write("manuscript/table_S4_long.csv", table_S4)

wide = unstack(table_S4, [:aversion, :scenario], :xt, :pct_increase)
rename!(wide, Dict(
    Symbol("0.25") => "xt_0.25",
    Symbol("0.5")  => "xt_0.5",
    Symbol("1.0")  => "xt_1.0",
    Symbol("2.0")  => "xt_2.0",
))
CSV.write("manuscript/table_S4.csv", wide)

println(wide)
wide
