#####################################################################
#####################################################################
### Demand metric 1: coverage levels.
###
### Two views of the same policy:
###
###   coverage_grid    optimal coverage over the whole (savings, biomass)
###                    grid, for heat maps
###   coverage_by_budget  optimal coverage integrated over the stationary
###                    distribution of biological states, as a function of
###                    the budget available above the borrowing floor
###
### Coverage is η / η̄ throughout: premium spending divided by the price
### per unit exposure, i.e. the units of exposure bought.
### Jack H. Buckner, Oregon State University, 2026
### Code generated with Calude 
#####################################################################
#####################################################################

include(joinpath(@__DIR__, "setup.jl"))


# ------------------------------------------------------------------
# Heat maps over the state grid
# ------------------------------------------------------------------

"""
    coverage_grid(prob)

Coverage over the full non bankrupt state grid.  Thin wrapper on
`coverage_levels` (index_model/analysis.jl), which returns
`(η, η̄, coverage, grid)` with `coverage` of size `(Ns, Nb, NΔb, 2)` and `grid`
the `(s, b, Δb, M)` axes read off the value function.
"""
coverage_grid(prob) = coverage_levels(prob)


# The model carries mortality only: the payout indicator is drawn fresh each
# period and never enters the state (see population_model.jl).
const M_STATE_LABELS = Dict(1 => "typical mortality", 2 => "high mortality")

# Index of the Δb axis closest to a target, defaulting to a stock that was flat
# last period.
function Δb_index(grid, Δb = 0.0)
    ax = collect(grid[3])
    return argmin(abs.(ax .- Δb))
end


"""
    coverage_heatmap(cl; Δb = 0.0, M = 1, kwargs...)

Heat map of optimal coverage over biomass (x) and budget above the borrowing
floor (y), sliced at one Δ log b and one mortality state.

Both axes are plotted in levels rather than logs: `exp` of the `b` and `s` grids.
"""
function coverage_heatmap(cl; Δb = 0.0, M::Int = 1, title = nothing, kwargs...)
    s_ax, b_ax = collect(cl.grid[1]), collect(cl.grid[2])
    k = Δb_index(cl.grid, Δb)
    # coverage[:, :, k, M] is (Ns, Nb) = (length(y), length(x)), which is the
    # orientation Plots.heatmap expects
    return Plots.heatmap(exp.(b_ax), exp.(s_ax), cl.coverage[:, :, k, M];
        xlabel = "Biomass", ylabel = "Budget (savings above the floor)",
        title = something(title, M_STATE_LABELS[M]),
        colorbar_title = "Coverage (units of exposure)", kwargs...)
end


"""
    coverage_heatmap_panel(cl; Δb = 0.0, states = (1, 2))

Heat maps for the typical and high mortality states side by side.
"""
function coverage_heatmap_panel(cl; Δb = 0.0, states = (1, 2), kwargs...)
    plts = [coverage_heatmap(cl; Δb = Δb, M = M, kwargs...) for M in states]
    return Plots.plot(plts..., layout = (1, length(plts)),
                      size = (420 * length(plts), 340))
end


# ------------------------------------------------------------------
# Coverage integrated over the stationary distribution of the stock
# ------------------------------------------------------------------

# Budget levels used by the base model's welfare table (src/analysis.jl:298),
# expressed as savings above the borrowing floor, s - s̄.
const BASE_MODEL_BUDGETS = [0.25, 0.5, 1.0, 2.0, 3.5]

"""
    default_budgets(g; n = 15)

Log spaced budget levels spanning the solved savings grid, so the sweep never
extrapolates outside the state space.
"""
function default_budgets(g; n::Int = 15)
    s_ax = g.state_grid[1]
    return exp.(range(first(s_ax), last(s_ax), length = n))
end


"""
    coverage_by_budget(prob, sim_states, budgets; nmax = 500)

Optimal coverage as a function of the budget available above the borrowing
floor, integrated over the stationary distribution of biological states.

For each budget `s - s̄`, the savings state is held at `log(s - s̄)` while
`(log b, Δ log b, M)` ranges over `sim_states`, and the optimal policy is
re-solved at every one of those states.  Returns one NamedTuple per budget.

`nmax` thins `sim_states`: `policy` re-solves the one period problem at each
state, so the cost is `length(budgets) * nmax` Bellman evaluations.
"""
function coverage_by_budget(prob, sim_states, budgets = default_budgets(prob);
                            nmax::Int = 500)
    states = thin_states(sim_states, nmax)
    n = size(states, 2)
    p = prob.p

    # One policy call over every (budget, state) pair: the solver threads
    # internally, so a single wide matrix beats a loop of narrow ones.
    S = zeros(5, n * length(budgets))
    for (i, s) in enumerate(budgets)
        cols = ((i - 1) * n + 1):(i * n)
        S[1, cols] .= 1.0            # not bankrupt
        S[2, cols] .= log(s)         # savings state is log(s - s̄)
        S[3:5, cols] .= states
    end

    u_opt = ValueFunctionIterations.policy(S, prob.u, prob.X, prob.R, prob.F,
                                           p, prob.δ, prob.V)
    η_all = u_opt[2, :]
    η̄_all = vec(mapslices(s -> price_per_exposure(El(s, p), p.cf, p.cv), S, dims = 1))

    out = NamedTuple[]
    for (i, s) in enumerate(budgets)
        cols = ((i - 1) * n + 1):(i * n)
        η, η̄ = η_all[cols], η̄_all[cols]
        cov = η ./ η̄
        # The action grid caps premium spending at ηmax (0.2) of the budget, so
        # coverage above ηmax*s/η̄ is not representable.  Where the cap binds the
        # reported optimum is censored by the grid rather than chosen, so track
        # how often that happens instead of leaving it implicit.
        ηcap = maximum(prob.u([1.0, log(s), 0.0, 0.0, 1.0], p)[2, :])
        push!(out, (
            budget           = s,
            savings          = s + p.s̄,          # savings on the model's own scale
            mean_coverage    = mean(cov),
            median_coverage  = median(cov),
            q10_coverage     = quantile(cov, 0.1),
            q90_coverage     = quantile(cov, 0.9),
            mean_premium     = mean(η),
            premium_share    = mean(η) / s,      # fraction of the budget spent
            premium_rate     = mean(η̄),          # price per unit exposure
            frac_buying      = mean(cov .> 0),
            max_coverage     = ηcap / mean(η̄),   # highest coverage the grid allows
            frac_at_cap      = mean(η .>= ηcap * (1 - 1e-9)),
            n_states         = n,
        ))
    end
    return out
end

# coverage_by_budget(prob, sim_states) needs the savings axis; pull it from the
# value function so the sweep matches the solved grid without passing g around.
default_budgets(prob::ValueFunctionIterations.DynamicProgram; n::Int = 15) =
    exp.(range(extrema(prob.V.states[2, prob.V.states[1, :] .== 1])..., length = n))


"""
    plot_coverage_by_budget(rows; label = "")

Mean coverage against budget, with a 10-90% band over the stationary
distribution of biological states.
"""
function plot_coverage_by_budget(rows; label = "", plt = nothing, kwargs...)
    b = [r.budget for r in rows]
    m = [r.mean_coverage for r in rows]
    lo = [r.q10_coverage for r in rows]
    hi = [r.q90_coverage for r in rows]
    f = plt === nothing ? Plots.plot : Plots.plot!
    out = f(b, m; ribbon = (m .- lo, hi .- m), label = label,
            xlabel = "Budget (savings above the floor)",
            ylabel = "Coverage (units of exposure)",
            xscale = :log10, linewidth = 2, fillalpha = 0.15, kwargs...)
    return out
end
