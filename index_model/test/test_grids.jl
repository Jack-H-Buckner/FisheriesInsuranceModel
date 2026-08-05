#####################################################################
#####################################################################
### Tests for the index insurance state/action grids and the grid
### metadata that guards against solving on one grid and analysing on
### another.
###
### Run from the repository root:
###     julia --project=. index_model/test/test_grids.jl
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

using Test, TOML

include(joinpath(@__DIR__, "..", "model.jl"))   # ValueFunction, init, population model
include(joinpath(@__DIR__, "..", "grids.jl"))
include(joinpath(@__DIR__, "..", "parameters", "generate_scenarios.jl"))

const SCEN = "no_RA_int_pr60_re100"


@testset "grids" begin

    @testset "sizes and bounds" begin
        g = build_grids(SCEN)
        s, b, Δb = g.state_grid
        κ, η = g.action_grid

        @test length(s)  == GRID_DEFAULTS.Ns
        @test length(b)  == GRID_DEFAULTS.Nb
        @test length(Δb) == GRID_DEFAULTS.NΔb
        @test length(κ)  == GRID_DEFAULTS.Nκ
        @test length(η)  == GRID_DEFAULTS.Nη

        # bounds hit exactly, matching the lo:step:hi idiom in src/run.jl
        @test (first(s),  last(s))  == (-4.0, 2.5)
        @test (first(b),  last(b))  == (-3.5, 0.5)
        @test (first(Δb), last(Δb)) == (-1.0, 1.25)
        @test (first(κ),  last(κ))  == (0.01, 0.99)
        @test (first(η),  last(η))  == (0.0, 0.3)

        @test g.δ == 0.95
        @test g.n_states == GRID_DEFAULTS.Ns * GRID_DEFAULTS.Nb * GRID_DEFAULTS.NΔb * 2 + 16
        @test g.n_actions == count(k + e <= 1 for k in κ, e in η)
        @test g.n_actions < length(κ) * length(η)   # the budget simplex bites
    end

    @testset "no insurance uses a degenerate premium grid" begin
        g = build_grids(SCEN; insurance = false)
        @test g.action_grid[2] == [0.0]
        @test g.sizes.Nη == 1                    # recorded honestly, not the ignored Nη
        @test g.bounds.ηmax == 0.0
        # state grid is untouched, so V has the same shape as the insurance solve
        @test g.state_grid == build_grids(SCEN).state_grid
    end

    @testset "sizes accept ARGS strings" begin
        a = build_grids(SCEN; Nκ = "10", Nη = "4", Ns = "6", Nb = "6", NΔb = "5", Nquad = "7")
        b = build_grids(SCEN; Nκ = 10, Nη = 4, Ns = 6, Nb = 6, NΔb = 5, Nquad = 7)
        @test a.state_grid == b.state_grid
        @test a.action_grid == b.action_grid
        @test a.sizes == b.sizes
    end

    @testset "grid sizes come from the config file" begin
        # the committed config is what build_grids falls back to
        @test read_grid_config() == GRID_DEFAULTS
        @test keys(GRID_DEFAULTS) == SIZE_KEYS
        @test all(v -> v isa Int && v >= 1, values(GRID_DEFAULTS))

        write_cfg(body) = (p = joinpath(mktempdir(), "config.toml"); write(p, body); p)
        full = """
               [sizes]
               "Nκ" = 10
               "Nη" = 4
               "Ns" = 6
               "Nb" = 6
               "NΔb" = 5
               Nquad = 7
               """

        @test read_grid_config(write_cfg(full)) ==
              (Nκ = 10, Nη = 4, Ns = 6, Nb = 6, NΔb = 5, Nquad = 7)
        # and those sizes are what the grid is actually built on
        @test build_grids(SCEN; read_grid_config(write_cfg(full))...).sizes ==
              (Nκ = 10, Nη = 4, Ns = 6, Nb = 6, NΔb = 5, Nquad = 7)

        # ASCII spellings, for a config written without quoting
        @test read_grid_config(write_cfg("""
              [sizes]
              Nkappa = 10
              Neta = 4
              Ns = 6
              Nb = 6
              NDeltab = 5
              Nquad = 7
              """)) == (Nκ = 10, Nη = 4, Ns = 6, Nb = 6, NΔb = 5, Nquad = 7)

        @test_throws ArgumentError read_grid_config(joinpath(mktempdir(), "nope.toml"))
        @test_throws ArgumentError read_grid_config(write_cfg("[other]\nx = 1\n"))
        # a partial config would silently mix configured and default sizes
        @test_throws ArgumentError read_grid_config(write_cfg("[sizes]\nNs = 6\n"))
        # a typo must not be ignored
        @test_throws ArgumentError read_grid_config(write_cfg(full * "Nsavings = 6\n"))
        # nor a size the grid constructors cannot use
        @test_throws ArgumentError read_grid_config(write_cfg(
            replace(full, "\"Nb\" = 6" => "\"Nb\" = 1")))
        @test_throws ArgumentError read_grid_config(write_cfg(
            replace(full, "\"Ns\" = 6" => "\"Ns\" = 6.5")))
        # an unparseable file is reported as such, with the quoting hint
        @test_throws ArgumentError read_grid_config(write_cfg("[sizes]\nNκ = 10\n"))
    end

    @testset "life history selection matches src/run.jl" begin
        @test life_history("no_RA_int_pr60_re100") === :intermediate
        @test life_history("no_RA_short_pr60_re100") === :short
        @test life_history("no_RA_long_pr60_re100")  === :long
        # every scenario in the sweep must resolve to :intermediate, otherwise
        # it would silently solve on the wrong b / Δb bounds
        for r in unique_scenarios()
            @test life_history(r.name) === :intermediate
        end
        @test state_bounds("x_short").bmax == 1.0
        @test state_bounds("x_long").bmin  == -4.5
    end

    @testset "grid too small is rejected" begin
        @test_throws ArgumentError build_grids(SCEN; Nb = 1)
        @test_throws ArgumentError gridrange(0.0, 1.0, 1)
    end

    @testset "metadata round trip" begin
        include(joinpath(@__DIR__, "..", "parameters", SCEN * ".jl"))   # defines p
        g = build_grids(SCEN)
        dir = mktempdir()
        path = joinpath(dir, "grid.toml")
        write_grid_meta(path, g, p, SCEN)
        meta = read_grid_meta(path)

        @test meta["scenario"] == SCEN
        @test meta["insurance"] == true
        @test meta["delta"] ≈ 0.95
        @test meta["life_history"] == "intermediate"
        @test meta["n_states"] == g.n_states

        # the 4x4 T survives the round trip, so analysis can confirm it is
        # reading the contract it thinks it is
        @test meta["parameters"]["T_size"] == [4, 4]
        @test reshape(meta["parameters"]["T"], 4, 4) ≈ p.T
        @test meta["parameters"]["p_low"]  ≈ 0.025
        @test meta["parameters"]["p_stay"] ≈ 0.5
        @test meta["parameters"]["cf"] ≈ p.cf
        @test meta["parameters"]["cv"] ≈ p.cv

        # rebuilding from the metadata reproduces the grid exactly
        g2 = grids_from_meta(meta)
        @test g2.state_grid  == g.state_grid
        @test g2.action_grid == g.action_grid
        @test g2.sizes == g.sizes
        @test assert_grid_match(meta, g) === nothing
    end

    @testset "assert_grid_match catches a mismatch" begin
        include(joinpath(@__DIR__, "..", "parameters", SCEN * ".jl"))
        g = build_grids(SCEN)
        meta = grid_meta(g, p, SCEN)
        @test assert_grid_match(meta, g) === nothing

        # a different state grid size must not be silently accepted
        @test_throws ErrorException assert_grid_match(meta, build_grids(SCEN; Nb = 20))
        # nor a different discount rate
        @test_throws ErrorException assert_grid_match(meta, build_grids(SCEN; δ = 0.9))
        # nor an insurance / no insurance mix up
        @test_throws ErrorException assert_grid_match(meta, build_grids(SCEN; insurance = false))
        # nor edited bound constants: simulate by tampering with the record
        tampered = deepcopy(meta)
        tampered["bounds"]["bmax"] = 0.75
        @test_throws ErrorException assert_grid_match(tampered, g)
    end

    @testset "grids drive a real DynamicProgram" begin
        include(joinpath(@__DIR__, "..", "parameters", SCEN * ".jl"))
        # small grid: this checks wiring, not the solution
        g = build_grids(SCEN; Nκ = 6, Nη = 4, Ns = 5, Nb = 6, NΔb = 4, Nquad = 5)
        X, Xmc = init_ranom_variables(p; Nquad = g.sizes.Nquad)
        prob = init(g.action_grid, g.state_grid, g.δ, p, X, index)

        @test prob isa ValueFunctionIterations.DynamicProgram
        # the value function carries the 4 joint mortality/index states
        @test size(prob.V.states, 1) == 5
        @test g.n_states == 5 * 6 * 4 * 2 + 16

        # F advances the state without error from an interior point
        node = X([0.0, 0.0, 0.0, 0.0, 1.0], p).nodes[:, 5]
        s1 = prob.F([0.0, 0.0, 0.0, 0.0, 1.0], [0.1, 0.1], node, p)
        @test length(s1) == 5
        @test all(isfinite, s1)
        @test s1[1] in (1.0, 2.0)          # bankruptcy indicator
        @test s1[5] in (1.0, 2.0, 3.0, 4.0) # joint mortality/index state
    end
end
