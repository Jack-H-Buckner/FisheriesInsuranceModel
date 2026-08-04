using DataFrames, CSV

# Table S3 reads percent-change-in-utility values that were computed by
# `value_diffs(params)` in src/analysis.jl and saved to
# `results/<params>/value_diff.csv`. Each csv has 5 rows of (savings, increase)
# pairs at savings st_vals = [0.25, 0.5, 1.0, 2.0, 3.5] .+ p.s̄. Table S3 uses
# the first four (xt ∈ {0.25, 0.5, 1.0, 2}); the 3.5 row is dropped.

xt_vals = [0.25, 0.5, 1.0, 2.0]

life_histories = [
    (label = "Intermediate life history", suffix = "intermidiate"),
    (label = "Fast life history",         suffix = "short"),
    (label = "Slow life history",         suffix = "long"),
]

aversion_blocks = [
    (label = "Bankruptcy",  prefix = "no_RA"),
    (label = "Risk-averse", prefix = "RA"),
]

scenarios = [
    (label = "Base risk",       tag = ""),
    (label = "High frequency",  tag = "_higher_frequency"),
    (label = "High mortality",  tag = "_higher_severity"),
    (label = "High recurrence", tag = "_longer_duration"),
]

# Long-format dataframe: one row per (life_history, aversion, scenario, xt).
table_S3 = DataFrame(
    life_history = String[],
    aversion     = String[],
    scenario     = String[],
    params       = String[],
    xt           = Float64[],
    pct_increase = Float64[],
)

for lh in life_histories, av in aversion_blocks, sc in scenarios
    params = string(av.prefix, "_", lh.suffix, sc.tag)
    path   = "results/$(params)/value_diff.csv"
    if !isfile(path)
        @warn "Missing value_diff.csv" params path
        for xt in xt_vals
            push!(table_S3, (lh.label, av.label, sc.label, params, xt, NaN))
        end
        continue
    end
    df = CSV.read(path, DataFrame)
    # value_diffs writes savings rows in fixed order [0.25,0.5,1.0,2.0,3.5] .+ s̄.
    # Take the first four rows positionally; the absolute `savings` column
    # depends on s̄ and is not what we want to key on.
    for (i, xt) in enumerate(xt_vals)
        push!(table_S3, (lh.label, av.label, sc.label, params, xt,
                         df.increase[i]))
    end
end

CSV.write("manuscript/table_S3_long.csv", table_S3)

# Wide format that mirrors the printed table layout:
# rows = (life_history, aversion, scenario), columns = xt budgets.
wide = unstack(table_S3, [:life_history, :aversion, :scenario], :xt,
               :pct_increase)
rename!(wide, Dict(
    Symbol("0.25") => "xt_0.25",
    Symbol("0.5")  => "xt_0.5",
    Symbol("1.0")  => "xt_1.0",
    Symbol("2.0")  => "xt_2.0",
))
CSV.write("manuscript/table_S3.csv", wide)

println(wide)
wide
