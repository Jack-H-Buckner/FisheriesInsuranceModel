#####################################################################
#####################################################################
### Tests for the runner CLI contract, the output directory layout,
### and the invariance that lets one no insurance solve be shared by
### every scenario on a base parameter set.
###
### Run from the repository root:
###     julia --project=. index_model/test/test_runners.jl
###
### The invariance testset actually solves, on a deliberately tiny
### grid, so this is slower than the other test files.
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
include(joinpath(@__DIR__, "..", "parameters", "generate_scenarios.jl"))

# Small enough to solve in seconds; this checks wiring, not economics.
const TINY = (Nκ = 5, Nη = 4, Ns = 4, Nb = 4, NΔb = 3, Nquad = 4)


@testset "runners" begin

    @testset "CLI argument contract" begin
        # sizes come from index_model/config.toml, not the command line
        a = parse_runner_args(["no_RA_int_pr60_re100"])
        @test a.params == "no_RA_int_pr60_re100"
        @test a.config == GRID_CONFIG_PATH
        @test a.sizes == GRID_DEFAULTS
        @test build_grids(a.params; a.sizes...).sizes == GRID_DEFAULTS

        # an explicit config path overrides the committed one
        cfg = joinpath(mktempdir(), "config.toml")
        write(cfg, """
              [sizes]
              "Nκ" = 5
              "Nη" = 4
              "Ns" = 4
              "Nb" = 4
              "NΔb" = 3
              Nquad = 4
              """)
        b = parse_runner_args(["no_RA_int_pr60_re100", cfg])
        @test b.config == cfg
        @test b.sizes == TINY

        @test_throws ArgumentError parse_runner_args(String[])
        # the old contract took 7 positional sizes; passing them must fail
        # loudly rather than be silently ignored
        @test_throws ArgumentError parse_runner_args(["s", "35", "25", "24", "24", "18", "50"])
        # a config path that does not exist is caught before anything is solved
        @test_throws ArgumentError parse_runner_args(["s", joinpath(mktempdir(), "nope.toml")])
    end

    @testset "unknown scenario fails fast" begin
        @test_throws ArgumentError scenario_row("no_RA_int_pr99_re99")
        @test_throws ArgumentError scenario_base("typo")
        @test scenario_base("no_RA_int_pr60_re100") == "no_RA_intermidiate"
        @test scenario_base("RA_int_pr50_re50") == "RA_intermidiate"
    end

    @testset "output directory layout" begin
        # insurance solves are keyed on the scenario
        @test basename(output_dir("no_RA_int_pr60_re100"; insurance = true)) ==
              "no_RA_int_pr60_re100"
        # no insurance solves collapse to one directory per base
        @test basename(output_dir("no_RA_int_pr60_re100"; insurance = false)) ==
              "no_insurance_no_RA_intermidiate"
        @test basename(output_dir("RA_int_pr50_re50"; insurance = false)) ==
              "no_insurance_RA_intermidiate"

        # every no_RA scenario maps to the same no insurance directory, and the
        # two bases never collide
        no_RA = unique([no_insurance_name(r.name) for r in unique_scenarios()
                        if scenario_base(r.name) == "no_RA_intermidiate"])
        RA    = unique([no_insurance_name(r.name) for r in unique_scenarios()
                        if scenario_base(r.name) == "RA_intermidiate"])
        @test length(no_RA) == 1
        @test length(RA) == 1
        @test isempty(intersect(no_RA, RA))
    end

    @testset "no insurance solution is invariant to the contract" begin
        # This is what makes one no insurance solve per base legitimate rather
        # than a shortcut.  Two scenarios on the same base, with very different
        # contracts, must produce the same value function when η is pinned to 0.
        dir = mktempdir()

        include(joinpath(@__DIR__, "..", "parameters", "no_RA_int_pr100_re100.jl"))
        p_a = copy(p)
        include(joinpath(@__DIR__, "..", "parameters", "no_RA_int_pr50_re100.jl"))
        p_b = copy(p)

        # different contracts ...
        @test !isapprox(p_a.T, p_b.T; atol = 1e-8)
        # ... but the same marginal mortality chain
        @test marginal_chain(p_a.T).p_low  ≈ marginal_chain(p_b.T).p_low
        @test marginal_chain(p_a.T).p_stay ≈ marginal_chain(p_b.T).p_stay

        prob_a, = solve_scenario(p_a, "no_RA_int_pr100_re100"; insurance = false,
                                 plots = false, outdir = joinpath(dir, "a"), TINY...)
        prob_b, = solve_scenario(p_b, "no_RA_int_pr50_re100"; insurance = false,
                                 plots = false, outdir = joinpath(dir, "b"), TINY...)

        # value functions agree over the whole grid
        states = prob_a.V.states
        va = [prob_a.V(states[:, i]) for i in 1:size(states, 2)]
        vb = [prob_b.V(states[:, i]) for i in 1:size(states, 2)]
        @test maximum(abs.(va .- vb)) < 1e-8

        # and with insurance switched on they must NOT agree, otherwise the
        # test above would pass for the wrong reason
        prob_c, = solve_scenario(p_a, "no_RA_int_pr100_re100"; insurance = true,
                                 plots = false, outdir = joinpath(dir, "c"), TINY...)
        vc = [prob_c.V(states[:, i]) for i in 1:size(states, 2)]
        @test maximum(abs.(va .- vc)) > 1e-8
    end

    @testset "solve writes solution and metadata" begin
        dir = mktempdir()
        include(joinpath(@__DIR__, "..", "parameters", "no_RA_int_pr60_re100.jl"))
        prob, g, outdir = solve_scenario(p, "no_RA_int_pr60_re100"; insurance = true,
                                         plots = false, outdir = dir, TINY...)
        @test isfile(joinpath(dir, "sol.jld2"))
        @test isfile(joinpath(dir, "grid.toml"))

        # the metadata describes the grid that was actually used, and an
        # analysis script rebuilding from it lands on the same grid
        meta = read_grid_meta(joinpath(dir, "grid.toml"))
        @test meta["scenario"] == "no_RA_int_pr60_re100"
        @test meta["insurance"] == true
        @test meta["sizes"]["Ns"] == TINY.Ns
        @test assert_grid_match(meta, g) === nothing
        @test grids_from_meta(meta).state_grid == g.state_grid

        # the solution reloads onto a freshly built problem
        X, _ = init_ranom_variables(p; Nquad = g.sizes.Nquad)
        fresh = init(g.action_grid, g.state_grid, g.δ, p, X, index)
        load_solution(fresh, joinpath(dir, "sol.jld2"))
        s = prob.V.states[:, 1]
        @test fresh.V(s) ≈ prob.V(s)
    end

    @testset "budget constraint holds across the production grid" begin
        # F used to recompute the budget as (exp(st)+s̄)-s̄, which can round below
        # the exp(st) that u scaled the actions by and throw on a full budget
        # action.  Check F accepts every action u offers, at the production Ns.
        include(joinpath(@__DIR__, "..", "parameters", "no_RA_int_pr60_re100.jl"))
        g = build_grids("no_RA_int_pr60_re100")
        X, _ = init_ranom_variables(p; Nquad = 5)
        prob = init(g.action_grid, g.state_grid, g.δ, p, X, index)

        node = X([0.0, 0.0, 0.0, 0.0, 1.0], p).nodes[:, 3]
        n = 0
        for st in g.state_grid[1]
            s = [1.0, st, 0.0, 0.0, 1.0]
            actions = prob.u(s, p)
            for c in 1:size(actions, 2)
                @test prob.F(s, actions[:, c], node, p) isa Vector   # must not throw
                n += 1
            end
        end
        @test n == length(g.state_grid[1]) * g.n_actions
    end
end
