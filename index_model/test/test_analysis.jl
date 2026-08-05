#####################################################################
#####################################################################
### Tests for the routines that read a solved index insurance
### problem: the simulation wrapper, the plots, and coverage_levels.
###
### These cover the defects fixed in the index_model clean up, so a
### regression shows up here rather than as a silently wrong figure.
###
### Run from the repository root:
###     julia --project=. index_model/test/test_analysis.jl
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

using Test

ENV["GKSwstype"] = "100"
using Plots, LaTeXStrings

include(joinpath(@__DIR__, "..", "model.jl"))
include(joinpath(@__DIR__, "..", "plotting.jl"))
include(joinpath(@__DIR__, "..", "grids.jl"))
include(joinpath(@__DIR__, "..", "solve.jl"))
include(joinpath(@__DIR__, "..", "analysis.jl"))
include(joinpath(@__DIR__, "..", "parameters", "generate_scenarios.jl"))
include(joinpath(@__DIR__, "..", "parameters", "no_RA_int_pr60_re100.jl"))

const TINY = (Nκ = 5, Nη = 4, Ns = 4, Nb = 4, NΔb = 3, Nquad = 4)
const NSIM = 40

prob, g, = solve_scenario(p, "no_RA_int_pr60_re100";
                          insurance = true, plots = false,
                          outdir = mktempdir(), TINY...)


@testset "reading a solved problem" begin

    X, Xmc = init_ranom_variables(p; Nquad = g.sizes.Nquad)
    x = simulation(prob, El, Xmc, NSIM; with_policy = false,
                   s0 = [1.0, 0.0, log(0.3), 0.0, 1.0])

    @testset "simulation series are length consistent" begin
        # simulate returns states of size (5, T+1) but actions and randos of
        # size (., T).  Every per period series must be length T; wt used to be
        # T-1, which silently truncated plot_simulation by one period.
        per_period = (:st, :bt, :bt1, :κt, :ηt, :Elt, :η̄t, :Zt, :wt,
                      :ft, :mt, :ht, :pt, :c_et, :yt)
        for k in per_period
            @test length(getproperty(x, k)) == NSIM
        end
        @test length(x.wt) == length(x.ηt)
        # Ibt is deliberately T+1: the bankruptcy rate calculation reads the
        # indicator of the final state.
        @test length(x.Ibt) == NSIM + 1
    end

    @testset "mortality uses the 4 state mapping" begin
        # environment_mortaltiy_index, not the 2 state environment_mortaltiy,
        # which would return m̲ for state 2 (typical mortality with a payout).
        @test all(m -> m ≈ p.m̄ || m ≈ p.m̲, x.mt)
        @test environment_mortaltiy_index(1, p.m̄, p.m̲) ≈ p.m̄
        @test environment_mortaltiy_index(2, p.m̄, p.m̲) ≈ p.m̄   # payout, typical M
        @test environment_mortaltiy_index(3, p.m̄, p.m̲) ≈ p.m̲
        @test environment_mortaltiy_index(4, p.m̄, p.m̲) ≈ p.m̲
    end

    @testset "plots take p explicitly" begin
        # plot_simulation used to read a global `p`, picking up whatever
        # parameter set happened to be in scope.
        @test plot_simulation(x, p, T = NSIM - 1) !== nothing
        @test plot_value_function(prob) !== nothing
    end

    @testset "coverage_levels" begin
        cl = coverage_levels(prob)

        @test size(cl.η) == (TINY.Ns, TINY.Nb, TINY.NΔb, 4)
        @test size(cl.η) == size(cl.η̄) == size(cl.coverage)
        @test all(isfinite, cl.coverage)
        @test all(cl.η .>= 0)

        # the axes come from the value function, so they agree with the grid the
        # problem was solved on; this used to return an undefined global
        @test length(cl.grid) == 4
        for d in 1:3
            @test collect(cl.grid[d]) ≈ collect(g.state_grid[d])
        end
        @test collect(cl.grid[4]) == collect(1.0:4.0)

        # El is state independent for this model, so the fair price is flat and
        # every bit of variation in `coverage` is demand rather than price
        m = index_moments(0.6, 1.0, 0.025, 0.5)
        @test all(cl.η̄ .≈ price_per_exposure(m.p_index, p.cf, p.cv))
        @test cl.coverage ≈ cl.η ./ cl.η̄
    end

    @testset "deterministic simulate_trajectory" begin
        # starts in log space, like the two stochastic methods
        st = simulate_trajectory(p, [0.0, 0.5], 5)
        @test st[1, 1] ≈ log(p.b_target)
        @test st[2, 1] == 0.0
        @test size(st) == (3, 7)
        # typed, so a mis-typed call to the stochastic methods errors instead of
        # silently falling through to the deterministic one
        @test_throws MethodError simulate_trajectory(p, "not a shock vector", 5)
    end

    @testset "joint_transition_from_matrix agrees with the 4 argument form" begin
        @test joint_transition(0.05, 0.8, 0.8, 0.4) ≈ joint_transition_from_matrix(
            [0.94/0.95 0.6; 1-0.94/0.95 0.4], [0.2 0.6; 0.8 0.4], [0.95 0.95; 0.05 0.05])
    end
end
