#####################################################################
#####################################################################
### Tests for the coverage analysis: rebuilding a solved problem,
### the stationary state cloud, the heat map grid, and coverage
### integrated over the stock distribution by budget level.
###
### Run from the repository root:
###     julia --project=. index_model/test/test_coverage.jl
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

using Test

ENV["GKSwstype"] = "100"
using Plots, LaTeXStrings

include(joinpath(@__DIR__, "..", "analysis", "coverage.jl"))

const SCEN = "no_RA_int_pr60_re100"
const TINY = (Nκ = 8, Nη = 6, Ns = 5, Nb = 5, NΔb = 4, Nquad = 5)

# Solve once into a temp directory, then rebuild from it exactly as an analysis
# script would.
const WORKDIR = mktempdir()
let q = scenario_parameters(SCEN)
    solve_scenario(q, SCEN; insurance = true, plots = false,
                   outdir = joinpath(WORKDIR, "ins"), TINY...)
    solve_scenario(q, SCEN; insurance = false, plots = false,
                   outdir = joinpath(WORKDIR, "noins"), TINY...)
end


@testset "coverage analysis" begin

    @testset "scenario_parameters reproduces the generated file" begin
        # the analysis rebuilds p in memory instead of including a scenario file
        # per scenario; it must agree exactly with what the file produces
        for name in [r.name for r in unique_scenarios()]
            include(joinpath(@__DIR__, "..", "parameters", name * ".jl"))  # defines p
            q = scenario_parameters(name)
            @test q.T ≈ p.T
            @test q.cf == p.cf
            @test q.cv == p.cv
            @test q.γ  == p.γ
            @test q.y0 == p.y0
        end
    end

    @testset "rebuild_problem" begin
        r = rebuild_problem(SCEN; insurance = true, dir = joinpath(WORKDIR, "ins"))
        @test r.g.sizes == (Nκ = TINY.Nκ, Nη = TINY.Nη, Ns = TINY.Ns,
                            Nb = TINY.Nb, NΔb = TINY.NΔb, Nquad = TINY.Nquad)
        @test r.meta["scenario"] == SCEN
        @test r.p.T ≈ scenario_parameters(SCEN).T
        # the loaded value function is not the zero initialisation.  The solved
        # values live on the per level B splines; DiscreteAndContinuous.values
        # is a zero initialised wrapper that solve! never writes back to.
        @test any(abs.(r.prob.V.Bslines[1].values) .> 0)

        # a missing solution is a clear error, not a silent empty result
        @test_throws ErrorException rebuild_problem(SCEN; dir = mktempdir())
    end

    @testset "parameter mismatch is caught" begin
        meta = read_grid_meta(joinpath(WORKDIR, "ins", "grid.toml"))
        good = scenario_parameters(SCEN)
        @test assert_parameters_match(meta, good; insurance = true) === nothing

        # a different contract on the same base must not pass as this scenario
        other = scenario_parameters("no_RA_int_pr50_re50")
        @test_throws ErrorException assert_parameters_match(meta, other; insurance = true)
        # ... but it is legitimate for the shared no insurance solution, which
        # depends only on the marginal chain
        @test assert_parameters_match(meta, other; insurance = false) === nothing
    end

    @testset "stationary states" begin
        r = rebuild_problem(SCEN; dir = joinpath(WORKDIR, "ins"))
        a = stationary_states(r.p, r.Xmc; n = 300, seed = 7)
        b = stationary_states(r.p, r.Xmc; n = 300, seed = 7)
        @test a == b                                   # reproducible
        @test size(a) == (3, 301)
        @test all(m -> m in 1.0:2.0, a[3, :])          # mortality state only
        @test all(isfinite, a)
        @test a[:, 1] ≈ [log(r.p.b_target), 0.0, 1.0]  # documented initial state

        c = stationary_states(r.p, r.Xmc; n = 300, seed = 7, burnin = 100)
        @test size(c, 2) == 201
        @test c == a[:, 101:end]
        @test_throws ArgumentError stationary_states(r.p, r.Xmc; n = 10, burnin = 99)

        @test size(thin_states(a, 50), 2) == 50
        @test thin_states(a, 10_000) === a             # no copy when under the cap
    end

    @testset "coverage grid for heat maps" begin
        r = rebuild_problem(SCEN; dir = joinpath(WORKDIR, "ins"))
        cl = coverage_grid(r.prob)

        @test size(cl.coverage) == (TINY.Ns, TINY.Nb, TINY.NΔb, 2)
        @test cl.coverage ≈ cl.η ./ cl.η̄
        @test all(cl.η .>= 0)
        @test all(isfinite, cl.coverage)

        # El is state independent for this model, so price is flat and every bit
        # of variation in coverage is demand rather than price
        m = index_moments(0.6, 1.0, 0.025, 0.5)
        @test all(cl.η̄ .≈ price_per_exposure(m.p_index, r.p.cf, r.p.cv))

        # the two mortality states are genuinely distinct: high mortality is a
        # different problem, so the policies must NOT coincide
        @test maximum(abs.(cl.coverage[:, :, :, 1] .- cl.coverage[:, :, :, 2])) > 1e-8

        @test coverage_heatmap(cl; M = 1) !== nothing
        @test coverage_heatmap(cl; M = 2) !== nothing
        @test coverage_heatmap_panel(cl; states = 1:2) !== nothing
        # Δb slice picks the axis point nearest the request
        @test Δb_index(cl.grid, -1.0) == 1
        @test Δb_index(cl.grid, 1e6) == length(cl.grid[3])
    end

    @testset "coverage by budget" begin
        r = rebuild_problem(SCEN; dir = joinpath(WORKDIR, "ins"))
        sim = stationary_states(r.p, r.Xmc; n = 200)
        rows = coverage_by_budget(r.prob, sim, BASE_MODEL_BUDGETS; nmax = 60)

        @test length(rows) == length(BASE_MODEL_BUDGETS)
        @test [x.budget for x in rows] == BASE_MODEL_BUDGETS
        @test all(x -> x.n_states == 60, rows)
        # savings is reported on the model's own scale, budget is above the floor
        @test all(x -> x.savings ≈ x.budget + r.p.s̄, rows)

        for x in rows
            @test 0 <= x.frac_buying <= 1
            @test 0 <= x.frac_at_cap <= 1
            @test x.q10_coverage <= x.median_coverage <= x.q90_coverage
            @test isfinite(x.mean_coverage)
            # the action grid caps premium spending, so coverage cannot exceed
            # ηmax * budget / η̄
            @test x.mean_coverage <= x.max_coverage + 1e-8
            @test x.q90_coverage  <= x.max_coverage + 1e-8
            # reference the constant, not a literal, so raising the cap does not
            # leave a stale bound silently passing here
            @test x.premium_share <= ACTION_BOUNDS.ηmax + 1e-8
        end

        # the cap scales with the budget, so it never binds for the same reason
        # at every level
        caps = [x.max_coverage for x in rows]
        @test issorted(caps)

        # default budgets stay inside the solved savings grid, so the sweep never
        # extrapolates outside the state space
        b = default_budgets(r.g; n = 9)
        @test length(b) == 9
        s_ax = r.g.state_grid[1]
        @test log(first(b)) ≈ first(s_ax)
        @test log(last(b))  ≈ last(s_ax)
        @test default_budgets(r.prob; n = 9) ≈ b

        @test plot_coverage_by_budget(rows; label = "test") !== nothing
    end

    @testset "no insurance solve buys nothing" begin
        r = rebuild_problem(SCEN; insurance = false, dir = joinpath(WORKDIR, "noins"))
        sim = stationary_states(r.p, r.Xmc; n = 100)
        rows = coverage_by_budget(r.prob, sim, [1.0]; nmax = 40)
        @test rows[1].mean_coverage == 0.0
        @test rows[1].frac_buying == 0.0
        @test rows[1].max_coverage == 0.0
    end
end
