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
#####################################################################
#####################################################################

# Precision ladder shared by both arms.  Spaced so that the arm B premium
# multiplier 1/pr lands on an even 1.00 / 1.33 / 1.67 / 2.00 sequence.
const PR_LADDER = [1.00, 0.75, 0.60, 0.50]

# Recall for the uninformative falsification scenario.  Its precision is not
# a free parameter: it is pinned to p_low, which is what makes the payout
# uninformative, so the generated file reads it out of the base parameter set.
const UNINFORMATIVE_RE = 0.50

const BASES = [
    (file = "no_RA_intermidiate", tag = "no_RA_int"),
    (file = "RA_intermidiate",    tag = "RA_int"),
]

const ARM_DOC = Dict(
    "iso_premium" =>
        "re = pr, so the payout rate is pinned to the event rate and the premium " *
        "is held at the perfect index level.  Isolates hedge quality from price.",
    "full_coverage" =>
        "re = 1, so every event is covered and declining precision buys that " *
        "coverage with false positives.  Isolates price from hedge quality.",
    "uninformative" =>
        "pr = p_low, so the payout carries no information about next period's " *
        "mortality.  Demand must be zero at cv = 0; falsification test, not a design point.",
)

tag100(x) = string(round(Int, 100 * x))

function scenario_name(base_tag, arm, pr, re)
    arm == "uninformative" && return "$(base_tag)_uninformative"
    return "$(base_tag)_pr$(tag100(pr))_re$(tag100(re))"
end


# Full design of the sweep, one row per (base, arm, rung).  The perfect index
# (pr = re = 1) belongs to both arms and so appears twice, under each arm
# label; it is solved once.  Derived quantities are deliberately absent -
# they depend on p_low and p_stay, which live in the base parameter file, and
# are recovered with index_moments at load time.
function scenario_table()
    rows = NamedTuple[]
    for b in BASES
        for pr in PR_LADDER
            push!(rows, (name = scenario_name(b.tag, "iso_premium", pr, pr),
                         base = b.file, base_tag = b.tag, arm = "iso_premium",
                         pr = pr, re = pr, kind = :design))
        end
        for pr in PR_LADDER
            push!(rows, (name = scenario_name(b.tag, "full_coverage", pr, 1.0),
                         base = b.file, base_tag = b.tag, arm = "full_coverage",
                         pr = pr, re = 1.0, kind = :design))
        end
        push!(rows, (name = scenario_name(b.tag, "uninformative", NaN, UNINFORMATIVE_RE),
                     base = b.file, base_tag = b.tag, arm = "uninformative",
                     pr = NaN, re = UNINFORMATIVE_RE, kind = :uninformative))
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


function scenario_source(r, arms)
    armline = join(arms, " and ")
    doc = join(["###            " * ARM_DOC[a] for a in arms], "\n")
    prline = r.kind == :uninformative ?
        "p = set_index_T(p; pr = p.T[2, 1], re = $(r.re), max_p_index = 0.6)" :
        "p = set_index_T(p; pr = $(r.pr), re = $(r.re))"
    prdesc = r.kind == :uninformative ? "p_low (uninformative)" : string(r.pr)

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
    ###   precision  $(prdesc)
    ###   recall     $(r.re)
    ### Jack H. Buckner, Oregon State University, 2026
    #####################################################################
    #####################################################################

    # parameters/*.jl call mean_price -> baranov, so the population model has to
    # be loaded first.  Use the index_model copy, never src/population_model.jl:
    # the latter would overwrite update_stock with the two state version.
    include(joinpath(@__DIR__, "..", "population_model.jl"))
    include(joinpath(@__DIR__, "..", "..", "parameters", "$(r.base).jl"))
    include(joinpath(@__DIR__, "index_universal_params.jl"))
    include(joinpath(@__DIR__, "..", "index_transitions.jl"))

    # p_low and p_stay are inherited from the base parameter set's 2x2 T.
    $(prline)
    p.cf = index_cf
    p.cv = index_cv
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
    include(joinpath(@__DIR__, "..", "index_transitions.jl"))
    include(joinpath(@__DIR__, "..", "population_model.jl"))
    rows = scenario_table()
    println(io, rpad("scenario", 26), rpad("arm", 16),
            join(lpad.(["pr", "re", "p_index", "ret_yr", "prem_x", "fpr", "b0"], 10)))
    for b in BASES
        include(joinpath(@__DIR__, "..", "..", "parameters", b.file * ".jl"))
        p_low, p_stay = p.T[2, 1], p.T[2, 2]
        for r in rows
            r.base == b.file || continue
            pr = r.kind == :uninformative ? p_low : r.pr
            m = index_moments(pr, r.re, p_low, p_stay)
            println(io, rpad(r.name, 26), rpad(r.arm, 16),
                    join(lpad.(string.(round.([pr, r.re, m.p_index, m.return_period,
                                               m.premium_ratio, m.fpr, m.b0], digits = 4)), 10)))
        end
    end
    return nothing
end


if abspath(PROGRAM_FILE) == @__FILE__
    write_scenarios()
    println()
    summarize_scenarios()
end
