#####################################################################
#####################################################################
### Defines the index insurance scenario sweep and writes one
### parameter file per scenario.
###
### This file is the single source of truth for the sweep.  Analysis
### scripts `include` it to get `scenario_table()` without writing
### anything; running it as a script regenerates the parameter files.
###
###     julia --project=. index_model/parameters/generate_scenarios.jl
###
### Jack H. Buckner, Oregon State University, 2026
### generated with Claude Code
#####################################################################
#####################################################################

# Loaded at top level, not inside a function: a package first loaded during a
# function call adds its methods to a newer world age than the caller, so
# ComponentArrays' getproperty would not be visible to summarize_scenarios.
include(joinpath(@__DIR__, "..", "population_model.jl"))
include(joinpath(@__DIR__, "..", "insurance_contracts.jl")) # price_per_exposure
include(joinpath(@__DIR__, "..", "index_transitions.jl"))

# Precision ladder shared by both arms.  Spaced so that the arm B premium
# multiplier 1/pr lands on an even 1.00 / 1.33 / 1.67 / 2.00 sequence.
const PR_LADDER = [1.00, 0.75, 0.6, 0.50, 0.42, 0.333]

# An uninformative falsification scenario (pr = re = p_low, which gives
# fpr = b0 = p_low, i.e. a payout statistically independent of the mortality
# event, at the same premium as the iso_premium arm) is deliberately not part
# of the sweep.  To add it back, push a row with kind = :uninformative and
# restore the corresponding branch in scenario_source.

const BASES = [
    (file = "no_RA_intermidiate", tag = "no_RA_int"),
    (file = "RA_intermidiate",    tag = "RA_int"),
]

# Insurance cost structure.  `nothing` means inherit from the base parameter
# set, which gets cf = 1e-6 and cv = 0.0 (actuarially fair) from
# parameters/universal_params.jl.  Inheriting is the default so those constants
# have exactly one home; an override is written into the scenario file only
# when it is actually different.
#
# Each variant is crossed with every arm and rung, and its `tag` is appended to
# the scenario name so a loaded premium is visible in the filename, the results
# directory and the manifest instead of silently changing every scenario at
# once.  To add a 10% load alongside the fair-price sweep, uncomment the second
# row: that doubles the job count and produces e.g. no_RA_int_pr60_re100_cv10.
#
# The premium actually charged is η̄ = (p_index + cf)/(1 - cv).
const COST_VARIANTS = [
    (tag = "",     cf = nothing, cv = nothing),  # inherit from the base set
    # (tag = "cv10", cf = nothing, cv = 0.10),   # 10% variable load
    (tag = "cv50", cf = nothing, cv = 0.50),   # 25% variable load
]

const ARM_DOC = Dict(
    "iso_premium" =>
        "re = pr, so the payout rate is pinned to the event rate and the premium " *
        "is held at the perfect index level.  Isolates hedge quality from price.",
    "full_coverage" =>
        "re = 1, so every event is covered and declining precision buys that " *
        "coverage with false positives.  Isolates price from hedge quality.",
)

tag100(x) = string(round(Int, 100 * x))

function scenario_name(base_tag, pr, re, cost_tag = "")
    name = "$(base_tag)_pr$(tag100(pr))_re$(tag100(re))"
    return isempty(cost_tag) ? name : "$(name)_$(cost_tag)"
end


# Full design of the sweep, one row per (base, arm, rung).  The perfect index
# (pr = re = 1) belongs to both arms and so appears twice, under each arm
# label; it is solved once.  Derived quantities are deliberately absent -
# they depend on p_low and p_stay, which live in the base parameter file, and
# are recovered with index_moments at load time.
function scenario_table()
    rows = NamedTuple[]
    for b in BASES, c in COST_VARIANTS
        for (arm, res) in (("iso_premium", PR_LADDER),        # re = pr
                           ("full_coverage", fill(1.0, length(PR_LADDER))))
            for (pr, re) in zip(PR_LADDER, res)
                push!(rows, (name = scenario_name(b.tag, pr, re, c.tag),
                             base = b.file, base_tag = b.tag, arm = arm,
                             pr = pr, re = re,
                             cost_tag = c.tag, cf = c.cf, cv = c.cv))
            end
        end
    end
    return rows
end


# One row per parameter file to write and per job to submit.
function unique_scenarios(rows = scenario_table())
    seen = Set{String}()
    out = NamedTuple[]
    for r in rows
        r.name in seen && continue
        push!(seen, r.name)
        push!(out, r)
    end
    return out
end

arms_for(name, rows = scenario_table()) = unique([r.arm for r in rows if r.name == name])

# Look up the design row for a scenario name.  Runners and analysis scripts go
# through this so a typo in a job submission fails immediately rather than
# producing an empty or wrongly named results directory.
function scenario_row(name, rows = scenario_table())
    i = findfirst(r -> r.name == name, rows)
    i === nothing && throw(ArgumentError(
        "unknown scenario \"$name\"; known scenarios: " *
        join(sort(unique([r.name for r in rows])), ", ")))
    return rows[i]
end

scenario_base(name, rows = scenario_table()) = scenario_row(name, rows).base

# The no insurance solution depends on neither the contract design nor its cost
# structure, so one solve is shared by every scenario on a given base:
#
#   (pr, re): with η_grid = [0.0] the index never enters F, and because T_x has
#             identical columns the x dimension is non persistent, so M' | M is
#             the marginal chain [1-p_low 1-p_stay; p_low p_stay] regardless.
#   (cf, cv): they enter only through η̄ = price_per_exposure(El, cf, cv), which
#             appears solely as ηt/η̄ in index_insurance.  With ηt = 0 the claim
#             is 0 for any η̄.
#
# What is left is p_low, p_stay and the economics of the base parameter set, all
# of which are fixed by the base file.  test_runners.jl asserts this rather than
# taking it on trust.
no_insurance_name(name, rows = scenario_table()) =
    "no_insurance_" * scenario_base(name, rows)


function scenario_source(r, arms)
    armline = join(arms, " and ")
    doc = join(["###            " * ARM_DOC[a] for a in arms], "\n")
    prline = "p = set_index_T(p; pr = $(r.pr), re = $(r.re))"

    # cf and cv are written only when the scenario overrides them, so that the
    # base parameter set stays the single home for the fair-price defaults.
    costlines = String[]
    r.cf === nothing || push!(costlines, "p.cf = $(r.cf)")
    r.cv === nothing || push!(costlines, "p.cv = $(r.cv)")
    costblock = isempty(costlines) ? "" : "\n" * join(costlines, "\n")
    costdesc  = isempty(costlines) ?
        "inherited from the base parameter set (actuarially fair)" :
        join([r.cf === nothing ? "" : "cf = $(r.cf)",
              r.cv === nothing ? "" : "cv = $(r.cv)"] |> x -> filter(!isempty, x), ", ") *
        " (override)"

    return """
    #####################################################################
    #####################################################################
    ### GENERATED FILE - do not edit by hand.
    ### Regenerate with:
    ###     julia --project=. index_model/parameters/generate_scenarios.jl
    ###
    ### Index insurance scenario: $(r.name)
    ###   base       $(r.base)
    ###   arm        $(armline)
    $(doc)
    ###   precision  $(r.pr)
    ###   recall     $(r.re)
    ###   cost       $(costdesc)
    ### Jack H. Buckner, Oregon State University, 2026
    #####################################################################
    #####################################################################

    # parameters/*.jl call mean_price -> baranov, so the population model has to
    # be loaded first.  Use the index_model copy, never src/population_model.jl:
    # the latter would overwrite update_stock with the two state version.
    include(joinpath(@__DIR__, "..", "population_model.jl"))
    include(joinpath(@__DIR__, "..", "..", "parameters", "$(r.base).jl"))
    include(joinpath(@__DIR__, "..", "index_transitions.jl"))

    # p_low and p_stay are inherited from the base parameter set's 2x2 T, and
    # cf / cv from parameters/universal_params.jl unless overridden below.
    $(prline)$(costblock)
    """
end


function write_scenarios(; dir = @__DIR__, verbose = true)
    rows = scenario_table()
    written = String[]
    for r in unique_scenarios(rows)
        path = joinpath(dir, r.name * ".jl")
        open(io -> write(io, scenario_source(r, arms_for(r.name, rows))), path, "w")
        push!(written, r.name)
        verbose && println("wrote ", relpath(path))
    end
    verbose && println("\n", length(written), " scenario files (",
                       length(rows), " design rows; the perfect index is shared by both arms)")
    return written
end


# Load each base parameter set and print the realised contract properties.
# Nothing here is hardcoded: p_low and p_stay come from the base files.
function summarize_scenarios(; io = stdout)
    rows = scenario_table()
    w = maximum(length.([r.name for r in rows])) + 2
    println(io, rpad("scenario", w), rpad("arm", 16),
            join(lpad.(["pr", "re", "p_index", "ret_yr", "fpr", "b0", "cv", "eta_bar"], 10)))
    for b in BASES
        include(joinpath(@__DIR__, "..", "..", "parameters", b.file * ".jl"))
        p_low, p_stay = p.T[2, 1], p.T[2, 2]
        for r in rows
            r.base == b.file || continue
            m = index_moments(r.pr, r.re, p_low, p_stay)
            # effective cost structure: the base set unless the scenario overrides
            cf = r.cf === nothing ? p.cf : r.cf
            cv = r.cv === nothing ? p.cv : r.cv
            η̄ = price_per_exposure(m.p_index, cf, cv)
            println(io, rpad(r.name, w), rpad(r.arm, 16),
                    join(lpad.(string.(round.([r.pr, r.re, m.p_index, m.return_period,
                                               m.fpr, m.b0, cv, η̄], digits = 4)), 10)))
        end
    end
    return nothing
end


if abspath(PROGRAM_FILE) == @__FILE__
    write_scenarios()
    println()
    summarize_scenarios()
end
