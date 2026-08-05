#####################################################################
#####################################################################
### State and action grids for the index insurance model, plus the
### metadata written alongside each solution.
###
### This is the SINGLE source of truth for grid construction.  Only V
### and P are serialised by save_solution, so every consumer of a
### sol.jld2 has to rebuild the problem on a byte identical grid or it
### silently reads the wrong numbers.  In the base model src/run.jl,
### src/analysis.jl and analysis/table_4.jl each hardcode different
### sizes; here the runners and every analysis script call build_grids,
### the sizes come from index_model/config.toml rather than from each
### caller's command line, and the sizes actually used are recorded
### next to the solution so a mismatch raises instead of going
### unnoticed.
### Jack H. Buckner, Oregon State University, 2026
### generated with Claude Code
#####################################################################
#####################################################################

using TOML

const DELTA = 0.95

# Action bounds, from src/run.jl:53-56.  ηmax applies to the insurance solve
# only; the no insurance solve uses the degenerate grid [0.0], matching
# src/run_no_insurance.jl:53.
const ACTION_BOUNDS = (κmin = 0.01, κmax = 0.99, ηmin = 0.0, ηmax = 0.3)


#####################################################################
### Grid sizes, read from index_model/config.toml
#####################################################################

# Sizes live in a config file rather than in the runners' command lines: the
# Makefile submits one job per scenario, and a sweep is only comparable if every
# job used the same grid.  See the header of config.toml.
const GRID_CONFIG_PATH = joinpath(@__DIR__, "config.toml")

const SIZE_KEYS = (:Nκ, :Nη, :Ns, :Nb, :NΔb, :Nquad)

# TOML bare keys are ASCII, so `Nκ = 35` is a parse error and the key has to be
# written `"Nκ" = 35`.  ASCII spellings are accepted so an unquoted config still
# works rather than failing on a rule that has nothing to do with the model.
const SIZE_NAMES = merge(Dict(string(k) => k for k in SIZE_KEYS),
                         Dict("Nkappa" => :Nκ, "Neta" => :Nη,
                              "NDeltab" => :NΔb, "Ndeltab" => :NΔb))

# Nquad is a node count, the rest are grid point counts and gridrange needs 2.
min_size(k::Symbol) = k === :Nquad ? 1 : 2


"""
    read_grid_config(path = GRID_CONFIG_PATH)

Read the grid sizes from a config file and return them as a NamedTuple ready to
splat into `build_grids` or `solve_scenario`.

Every size must be present: a partially specified config would silently mix
configured and default sizes, which is the failure this file exists to prevent.
Unknown keys are an error too, so a typo is caught rather than ignored.
"""
function read_grid_config(path::AbstractString = GRID_CONFIG_PATH)
    isfile(path) || throw(ArgumentError(
        "no grid config at $path; the index model reads its grid sizes from " *
        "index_model/config.toml"))

    cfg = try
        TOML.parsefile(path)
    catch err
        err isa TOML.ParserError || rethrow()
        throw(ArgumentError(
            "could not parse $path: $(err.type)  (line $(err.line))\n" *
            "note: TOML bare keys are ASCII only, so Greek names must be quoted, " *
            "as in \"Nκ\" = 35"))
    end

    sizes = get(cfg, "sizes", nothing)
    sizes isa AbstractDict || throw(ArgumentError(
        "$path has no [sizes] table; expected keys " * join(string.(SIZE_KEYS), ", ")))

    out = Dict{Symbol, Int}()
    for (name, v) in sizes
        k = get(SIZE_NAMES, name, nothing)
        k === nothing && throw(ArgumentError(
            "unknown size `$name` in $path; expected one of " *
            join(string.(SIZE_KEYS), ", ")))
        haskey(out, k) && throw(ArgumentError("$k given twice in $path"))
        v isa Integer || throw(ArgumentError(
            "$k in $path must be an integer, got $(repr(v))"))
        v >= min_size(k) || throw(ArgumentError(
            "$k in $path must be at least $(min_size(k)), got $v"))
        out[k] = v
    end

    absent = [string(k) for k in SIZE_KEYS if !haskey(out, k)]
    isempty(absent) || throw(ArgumentError(
        "$path is missing " * join(absent, ", ") *
        "; all six sizes must be given, so a solve never mixes configured and " *
        "default grids"))

    return NamedTuple{SIZE_KEYS}(Tuple(out[k] for k in SIZE_KEYS))
end

# Sizes used when a caller does not pass its own.  Read once, at include time,
# so the runners, the tests and every analysis script agree by construction.
const GRID_DEFAULTS = read_grid_config()


# Life history is selected by substring match on the scenario name, exactly as
# src/run.jl:26-31 does.  Scenario names generated for the index sweep are of
# the form no_RA_int_pr60_re100 and so resolve to :intermediate; a future base
# named *_short or *_long picks up the matching bounds automatically.
function life_history(params::AbstractString)
    occursin("short", params) && return :short
    occursin("long", params)  && return :long
    return :intermediate
end

# State bounds, from src/run.jl:24-32 and 44.
function state_bounds(params::AbstractString)
    lh = life_history(params)
    lh === :short && return (bmin = -3.0, bmax = 1.0,  Δbmin = -1.25, Δbmax = 2.0,
                             smin = -4.0, smax = 2.5)
    lh === :long  && return (bmin = -4.5, bmax = 0.0,  Δbmin = -0.75, Δbmax = 0.75,
                             smin = -4.0, smax = 2.5)
    return                  (bmin = -3.5, bmax = 0.5,  Δbmin = -1.0,  Δbmax = 1.25,
                             smin = -4.0, smax = 2.5)
end


asint(x::Integer) = Int(x)
asint(x::AbstractString) = parse(Int, x)

# Matches the `lo:((hi-lo)/(n-1)):hi` idiom used throughout src/.
function gridrange(lo, hi, n)
    n >= 2 || throw(ArgumentError("grid needs at least 2 points, got $n"))
    return lo:((hi - lo) / (n - 1)):hi
end


"""
    build_grids(params; Nκ, Nη, Ns, Nb, NΔb, Nquad, insurance = true, δ = DELTA)

Build the action and state grids for scenario `params`.  Sizes default to the
ones in `index_model/config.toml` and accept either integers or strings, so a
size read back out of a `grid.toml` can be passed straight through.

`insurance = false` replaces the premium grid with the degenerate `[0.0]`, which
is what makes the no insurance solve a genuine no insurance solve rather than a
cheap-insurance one.

Returns a NamedTuple carrying the grids plus the sizes and bounds that produced
them, so the caller can hand the whole thing to `write_grid_meta`.
"""
function build_grids(params::AbstractString;
                     Nκ    = GRID_DEFAULTS.Nκ,
                     Nη    = GRID_DEFAULTS.Nη,
                     Ns    = GRID_DEFAULTS.Ns,
                     Nb    = GRID_DEFAULTS.Nb,
                     NΔb   = GRID_DEFAULTS.NΔb,
                     Nquad = GRID_DEFAULTS.Nquad,
                     insurance::Bool = true,
                     δ = DELTA)

    Nκ, Nη, Ns, Nb, NΔb, Nquad = asint.((Nκ, Nη, Ns, Nb, NΔb, Nquad))
    sb = state_bounds(params)
    ab = ACTION_BOUNDS

    s_grid  = gridrange(sb.smin,  sb.smax,  Ns)
    b_grid  = gridrange(sb.bmin,  sb.bmax,  Nb)
    Δb_grid = gridrange(sb.Δbmin, sb.Δbmax, NΔb)
    κ_grid  = gridrange(ab.κmin,  ab.κmax,  Nκ)
    η_grid  = insurance ? gridrange(ab.ηmin, ab.ηmax, Nη) : [0.0]

    # Effective sizes: the no insurance premium grid is a single point, so
    # record 1 rather than the Nη that was passed in and ignored.
    Nη_eff = length(η_grid)

    # Number of actions actually searched, after the budget simplex filter in
    # index_model/model.jl:90, and number of continuous states.  Reported in the
    # job log because together they set the solve cost.
    n_actions = count(κ + η <= 1 for κ in κ_grid, η in η_grid)
    # 2 mortality states, + the 16 point degenerate bankrupt branch.  The
    # contract's T is 4x4 over the joint (mortality, payout) outcome, but the
    # payout indicator is never carried as a state; see population_model.jl.
    n_states  = Ns * Nb * NΔb * 2 + 16

    return (action_grid = [κ_grid, η_grid],
            state_grid  = [s_grid, b_grid, Δb_grid],
            δ           = δ,
            insurance   = insurance,
            sizes       = (Nκ = Nκ, Nη = Nη_eff, Ns = Ns, Nb = Nb, NΔb = NΔb, Nquad = Nquad),
            bounds      = merge(sb, (κmin = ab.κmin, κmax = ab.κmax,
                                     ηmin = ab.ηmin,
                                     ηmax = insurance ? ab.ηmax : 0.0)),
            n_actions   = n_actions,
            n_states    = n_states,
            life_history = string(life_history(params)))
end


# Marginal mortality chain implied by a 4x4 joint T, so the metadata records
# the event frequency and duration the solve actually ran with.
function marginal_chain(T)
    size(T) == (2, 2) && return (p_low = T[2, 1], p_stay = T[2, 2])
    return (p_low = T[3, 1] + T[4, 1], p_stay = T[3, 3] + T[4, 3])
end


"""
    grid_meta(g, p, params)

Metadata describing a solve: the grid that produced it, and enough of the
parameter set to confirm an analysis script is reading the solution it thinks
it is.  `save_solution` persists only V and P, so this is the record of
everything else.
"""
function grid_meta(g, p, params::AbstractString)
    mc = marginal_chain(p.T)
    return Dict(
        "scenario"     => params,
        "insurance"    => g.insurance,
        "delta"        => g.δ,
        "life_history" => g.life_history,
        "n_actions"    => g.n_actions,
        "n_states"     => g.n_states,
        "sizes"  => Dict(string(k) => v for (k, v) in pairs(g.sizes)),
        "bounds" => Dict(string(k) => v for (k, v) in pairs(g.bounds)),
        "parameters" => Dict(
            "T"      => collect(vec(p.T)),   # column major
            "T_size" => collect(size(p.T)),
            "p_low"  => mc.p_low,
            "p_stay" => mc.p_stay,
            "cf"     => p.cf,
            "cv"     => p.cv,
            "gamma"  => p.γ,
            "y0"     => p.y0,
            "s_bar"  => p.s̄,
        ),
        "provenance" => Dict(
            "julia"    => string(VERSION),
            "nthreads" => Threads.nthreads(),
            "written"  => Libc.strftime("%Y-%m-%d %H:%M:%S", time()),
        ),
    )
end

write_grid_meta(path::AbstractString, g, p, params::AbstractString) =
    open(io -> TOML.print(io, grid_meta(g, p, params)), path, "w")

read_grid_meta(path::AbstractString) = TOML.parsefile(path)


"""
    grids_from_meta(meta; check_bounds = true)

Rebuild the grids a solution was produced with, straight from its grid.toml.
By default this also re-derives the bounds from the scenario name and checks
them against the recorded ones, so editing the bound constants in this file
after a solve is caught rather than silently changing what analysis assumes.
"""
function grids_from_meta(meta::AbstractDict; check_bounds = true)
    sz = meta["sizes"]
    g = build_grids(meta["scenario"];
                    Nκ = sz["Nκ"], Nη = sz["Nη"], Ns = sz["Ns"],
                    Nb = sz["Nb"], NΔb = sz["NΔb"], Nquad = sz["Nquad"],
                    insurance = meta["insurance"], δ = meta["delta"])
    check_bounds && assert_grid_match(meta, g)
    return g
end


"""
    assert_grid_match(meta, g)

Throw unless the grid `g` is the one recorded in `meta`.  Every analysis script
runs this before trusting a loaded value function.
"""
function assert_grid_match(meta::AbstractDict, g)
    bad = String[]

    meta["insurance"] == g.insurance ||
        push!(bad, "insurance: recorded $(meta["insurance"]), rebuilt $(g.insurance)")
    meta["delta"] ≈ g.δ ||
        push!(bad, "delta: recorded $(meta["delta"]), rebuilt $(g.δ)")

    for (k, v) in pairs(g.sizes)
        rec = get(meta["sizes"], string(k), nothing)
        rec == v || push!(bad, "sizes.$k: recorded $rec, rebuilt $v")
    end
    for (k, v) in pairs(g.bounds)
        rec = get(meta["bounds"], string(k), nothing)
        (rec !== nothing && rec ≈ v) ||
            push!(bad, "bounds.$k: recorded $rec, rebuilt $v")
    end

    # n_states is not implied by the sizes above: it also depends on how many
    # mortality states the model carries.  Checking it is what catches a
    # solution produced before the model dropped the redundant payout state,
    # which would otherwise load silently and shift M = 2 from high mortality to
    # typical-mortality-with-a-payout.
    rec_n = get(meta, "n_states", nothing)
    rec_n == g.n_states || push!(bad,
        "n_states: recorded $rec_n, rebuilt $(g.n_states)" *
        (rec_n == 2 * g.n_states - 16 ?
         " (the solution predates the 2 state model; re-solve it)" : ""))

    isempty(bad) || throw(ErrorException(
        "grid mismatch for scenario $(meta["scenario"]) - the loaded value function was " *
        "solved on a different grid:\n  " * join(bad, "\n  ")))
    return nothing
end


# Compact description for the job log.
function print_grid(g; io = stdout)
    s, b, Δb = g.state_grid
    κ, η = g.action_grid
    println(io, "grid (", g.life_history, ", ",
            g.insurance ? "with insurance" : "no insurance", ", δ = ", g.δ, ")")
    println(io, "  s   ", length(s),  " in [", first(s),  ", ", last(s),  "]")
    println(io, "  b   ", length(b),  " in [", first(b),  ", ", last(b),  "]")
    println(io, "  Δb  ", length(Δb), " in [", first(Δb), ", ", last(Δb), "]")
    println(io, "  κ   ", length(κ),  " in [", first(κ),  ", ", last(κ),  "]")
    println(io, "  η   ", length(η),  " in [", first(η),  ", ", last(η),  "]")
    println(io, "  Nquad ", g.sizes.Nquad,
            "  states ", g.n_states, "  actions ", g.n_actions)
    return nothing
end
