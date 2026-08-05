#####################################################################
#####################################################################
### Shared setup for the index insurance analysis scripts.
###
### Rebuilds a solved problem from its results directory and produces
### the stationary distribution of biological states that every demand
### metric is averaged over.
###
### save_solution persists only V and P, so a consumer has to rebuild
### the problem on the same grid.  rebuild_problem does that from the
### grid.toml written next to the solution and asserts the match, so a
### stale or mismatched solution raises instead of quietly producing
### the wrong numbers.
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

using Random, Statistics

include(joinpath(@__DIR__, "..", "model.jl"))
include(joinpath(@__DIR__, "..", "plotting.jl"))
include(joinpath(@__DIR__, "..", "grids.jl"))
include(joinpath(@__DIR__, "..", "solve.jl"))
include(joinpath(@__DIR__, "..", "analysis.jl"))
include(joinpath(@__DIR__, "..", "parameters", "generate_scenarios.jl"))

# Base parameter sets, loaded once at top level.  Scenario parameters are then
# rebuilt in memory with set_index_T rather than by including a scenario file
# per scenario: a parameter file loaded inside a function, or inside a loop
# body, redefines methods in a newer world age than the running caller.
# test_coverage.jl checks this reproduces the generated files exactly.
const BASE_PARAMETERS = Dict{String,Any}()
for b in BASES
    include(joinpath(@__DIR__, "..", "..", "parameters", b.file * ".jl"))
    BASE_PARAMETERS[b.file] = copy(p)
end


"""
    scenario_parameters(scenario)

Rebuild the parameter set for `scenario` in memory.  Equivalent to including
`index_model/parameters/\$scenario.jl`, without the include.
"""
function scenario_parameters(scenario::AbstractString)
    r = scenario_row(scenario)
    q = set_index_T(BASE_PARAMETERS[r.base]; pr = r.pr, re = r.re)
    r.cf === nothing || (q.cf = r.cf)
    r.cv === nothing || (q.cv = r.cv)
    return q
end


# The no insurance solution is shared across contracts, so its grid.toml records
# whichever scenario happened to produce it.  Only the quantities that solution
# genuinely depends on are checked in that case.
function assert_parameters_match(meta, p; insurance::Bool)
    rec = meta["parameters"]
    mc = marginal_chain(p.T)
    bad = String[]

    isapprox(rec["p_low"],  mc.p_low;  atol = 1e-12) ||
        push!(bad, "p_low: recorded $(rec["p_low"]), rebuilt $(mc.p_low)")
    isapprox(rec["p_stay"], mc.p_stay; atol = 1e-12) ||
        push!(bad, "p_stay: recorded $(rec["p_stay"]), rebuilt $(mc.p_stay)")
    isapprox(rec["gamma"], p.γ; atol = 1e-12) ||
        push!(bad, "gamma: recorded $(rec["gamma"]), rebuilt $(p.γ)")
    isapprox(rec["y0"], p.y0; atol = 1e-12) ||
        push!(bad, "y0: recorded $(rec["y0"]), rebuilt $(p.y0)")

    if insurance
        # the contract itself only matters when insurance can be bought
        isapprox(reshape(rec["T"], 4, 4), p.T; atol = 1e-12) ||
            push!(bad, "T does not match the recorded transition matrix")
        isapprox(rec["cf"], p.cf; atol = 1e-12) ||
            push!(bad, "cf: recorded $(rec["cf"]), rebuilt $(p.cf)")
        isapprox(rec["cv"], p.cv; atol = 1e-12) ||
            push!(bad, "cv: recorded $(rec["cv"]), rebuilt $(p.cv)")
    end

    isempty(bad) || throw(ErrorException(
        "parameter mismatch for $(meta["scenario"]): the solution on disk was " *
        "produced with different parameters:\n  " * join(bad, "\n  ")))
    return nothing
end


"""
    rebuild_problem(scenario; insurance = true)

Rebuild the dynamic program for `scenario` and load its solved value function.

Returns `(prob, p, g, X, Xmc, meta, dir)`.  Raises if the grid or the parameters
recorded next to the solution disagree with what is rebuilt here.
"""
function rebuild_problem(scenario::AbstractString; insurance::Bool = true, dir = nothing)
    dir = something(dir, output_dir(scenario; insurance = insurance))
    sol = joinpath(dir, "sol.jld2")
    isfile(sol) || throw(ErrorException(
        "no solution at $sol; run index_model/run" *
        (insurance ? "" : "_no_insurance") * ".jl for $scenario first"))

    meta = read_grid_meta(joinpath(dir, "grid.toml"))
    p = scenario_parameters(scenario)
    g = grids_from_meta(meta)                       # asserts grid + bounds
    assert_parameters_match(meta, p; insurance = insurance)

    X, Xmc = init_ranom_variables(p; Nquad = g.sizes.Nquad)
    prob = init(g.action_grid, g.state_grid, g.δ, p, X, index)
    load_solution(prob, sol)

    return (prob = prob, p = p, g = g, X = X, Xmc = Xmc, meta = meta, dir = dir)
end


"""
    stationary_states(p, Xmc; n = 2000, seed = 123, burnin = 0)

Simulate the biological model and return the resulting cloud of states as a
`3 x m` matrix of `(log b, Δ log b, M)`, the empirical stationary distribution
that the demand metrics are averaged over.

Defaults match `src/analysis.jl:288` (n = 2000, seed = 123, no burn in) so the
index numbers stay comparable with the base model.  The chain starts at
`(log b_target, 0, M = 1)`; raise `burnin` to discard the approach to
stationarity.
"""
function stationary_states(p, Xmc; n::Int = 2000, seed::Int = 123, burnin::Int = 0)
    Random.seed!(seed)
    states, m, f, h, y = simulate_trajectory(p, Xmc, n)
    burnin < size(states, 2) || throw(ArgumentError(
        "burnin = $burnin exceeds the $(size(states,2)) simulated states"))
    return states[:, (burnin + 1):end]
end


# Evenly thin a state cloud to at most `nmax` columns.  Averages over the
# stationary distribution converge quickly, and `policy` re-solves the one
# period problem at every state, so thinning keeps the sweep affordable.
function thin_states(sim_states, nmax::Int)
    n = size(sim_states, 2)
    n <= nmax && return sim_states
    return sim_states[:, round.(Int, range(1, n, length = nmax))]
end
