#####################################################################
#####################################################################
### Tests for the index insurance transition matrix.
###
### Run from the repository root:
###     julia --project=. index_model/test/test_transitions.jl
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

using Test, Random, LinearAlgebra, StatsBase

include(joinpath(@__DIR__, "..", "index_transitions.jl"))

# baseline mortality chain from parameters/universal_params.jl
const P_LOW  = 0.025
const P_STAY = 0.50

# the two arms of the sweep, sharing one precision ladder
const PR_LADDER = [1.00, 0.75, 0.60, 0.50]
const ARMS = vcat(
    [(arm = "iso_premium",   pr = pr, re = pr)  for pr in PR_LADDER],
    [(arm = "full_coverage", pr = pr, re = 1.0) for pr in PR_LADDER if pr < 1.0],
)


@testset "index transition matrix" begin

    @testset "validity over both arms" begin
        for a in ARMS
            T = joint_transition(P_LOW, a.pr, a.re, P_STAY)
            @test size(T) == (4, 4)
            @test all(0 .<= T .<= 1)
            # MarkovChain in src/random_variables.jl throws above 1e-6
            @test maximum(abs.(vec(sum(T, dims = 1)) .- 1)) < 1e-12
            # validate_joint_transition must accept every scenario we ship
            @test validate_joint_transition(T, P_LOW, P_STAY) === nothing
        end
    end

    @testset "mortality marginal is invariant to (pr, re)" begin
        T_M = [1-P_LOW 1-P_STAY; P_LOW P_STAY]
        for a in ARMS
            T = joint_transition(P_LOW, a.pr, a.re, P_STAY)
            marginal = vcat((T[1, :] .+ T[2, :])', (T[3, :] .+ T[4, :])')
            @test marginal ≈ T_M[:, [1, 1, 2, 2]] atol = 1e-12
        end
    end

    @testset "El is state independent and equals p_index" begin
        for a in ARMS
            T = joint_transition(P_LOW, a.pr, a.re, P_STAY)
            m = index_moments(a.pr, a.re, P_LOW, P_STAY)
            El = [sum(T[[2, 4], j]) for j in 1:4]
            @test all(El .≈ m.p_index)
            # strictly positive, so price_per_exposure never hits its nugget
            @test all(El .> 1e-12)
        end
    end

    @testset "precision and recall are recovered from M = 1" begin
        for a in ARMS
            d = index_diagnostics(joint_transition(P_LOW, a.pr, a.re, P_STAY))
            m = index_moments(a.pr, a.re, P_LOW, P_STAY)
            for j in 1:2                      # states with M = 1
                @test d.precision[j] ≈ a.pr atol = 1e-12
                @test d.recall[j]    ≈ a.re  atol = 1e-12
                @test d.b0[j]        ≈ m.b0  atol = 1e-12
                @test d.fpr[j]       ≈ m.fpr atol = 1e-12
            end
            for j in 3:4                      # M = 2, index uninformative by design
                @test d.precision[j] ≈ P_STAY      atol = 1e-12
                @test d.recall[j]    ≈ m.p_index   atol = 1e-12
            end
        end
    end

    @testset "arm properties" begin
        # arm A holds the premium at the onset rate and makes the two error
        # rates equal
        for pr in PR_LADDER
            m = index_moments(pr, pr, P_LOW, P_STAY)
            @test m.p_index ≈ P_LOW atol = 1e-12
            @test m.premium_ratio ≈ 1.0 atol = 1e-12
            @test m.fpr ≈ m.b0 atol = 1e-12
        end
        # arm B covers every event and prices at 1/pr times the perfect index
        for pr in PR_LADDER
            m = index_moments(pr, 1.0, P_LOW, P_STAY)
            @test m.b0 ≈ 0.0 atol = 1e-12
            @test m.premium_ratio ≈ 1 / pr atol = 1e-12
        end
    end

    @testset "corner cases" begin
        # perfect index: payout if and only if onset
        T = joint_transition(P_LOW, 1.0, 1.0, P_STAY)
        @test T[2, 1] ≈ 0.0 atol = 1e-12   # payout without event
        @test T[3, 1] ≈ 0.0 atol = 1e-12   # event without payout
        @test T[4, 1] ≈ P_LOW atol = 1e-12
        @test index_moments(1.0, 1.0, P_LOW, P_STAY).b0 ≈ 0.0 atol = 1e-12

        # uninformative index: the payout carries no information about M'
        p_index = P_LOW * (1 - P_LOW) / P_LOW
        b0 = P_LOW * (1 - 0.5) / (1 - 0.5)
        T_mx0 = [(1-b0) (1-P_STAY); b0 P_STAY]
        T_mx1 = [(1-P_LOW) (1-P_STAY); P_LOW P_STAY]
        @test T_mx0[:, 1] ≈ T_mx1[:, 1] atol = 1e-12
        Tu = joint_transition(P_LOW, P_LOW, 0.5, P_STAY)
        du = index_diagnostics(Tu)
        @test du.precision[1] ≈ P_LOW atol = 1e-12   # precision = base rate
        @test du.b0[1] ≈ P_LOW atol = 1e-12          # payout tells you nothing
    end

    @testset "notebook anchors" begin
        # test_intermidiate.ipynb cell 4
        T1 = joint_transition(0.05, 0.8, 0.8, 0.4)
        @test T1 ≈ [0.94 0.94 0.57 0.57;
                    0.01 0.01 0.03 0.03;
                    0.01 0.01 0.38 0.38;
                    0.04 0.04 0.02 0.02] atol = 1e-10

        # test_transition_matrices.ipynb cell 2
        T2 = joint_transition(0.05, 1.0, 0.5, 0.5)
        @test T2 ≈ [0.95  0.95  0.4875 0.4875;
                    0.0   0.0   0.0125 0.0125;
                    0.025 0.025 0.4875 0.4875;
                    0.025 0.025 0.0125 0.0125] atol = 1e-10
    end

    @testset "monte carlo recovery of pr and re" begin
        Random.seed!(20260804)
        pr, re = 0.60, 0.75
        T = joint_transition(P_LOW, pr, re, P_STAY)
        n = 2_000_000
        j = 1
        n_payout = 0; n_event = 0; n_both = 0; n_fromlow = 0
        for _ in 1:n
            from_low = j <= 2
            j = sample(1:4, Weights(T[:, j]))
            event = j >= 3; payout = iseven(j)
            if from_low                       # pr and re are conditional on M = 1
                n_fromlow += 1
                n_payout += payout; n_event += event; n_both += (event && payout)
            end
        end
        @test n_both / n_payout ≈ pr rtol = 0.05
        @test n_both / n_event  ≈ re rtol = 0.05
    end

    @testset "validation rejects bad designs" begin
        @test_throws ArgumentError validate_index_parameters(0.0, 0.5, P_LOW, P_STAY)
        @test_throws ArgumentError validate_index_parameters(1.5, 0.5, P_LOW, P_STAY)
        @test_throws ArgumentError validate_index_parameters(0.5, 0.0, P_LOW, P_STAY)
        @test_throws ArgumentError validate_index_parameters(0.5, 1.5, P_LOW, P_STAY)
        @test_throws ArgumentError validate_index_parameters(0.5, 0.5, P_LOW, 1.0)
        # p_index >= 1 is fatal, it is the denominator of b0
        @test_throws ArgumentError validate_index_parameters(0.01, 1.0, 0.5, P_STAY)
        # a non stochastic matrix must be caught
        @test_throws ArgumentError validate_joint_transition(fill(0.1, 4, 4), P_LOW, P_STAY)
    end

    @testset "set_index_T on a real parameter set" begin
        # the parameter files call mean_price, which needs baranov, so the
        # population model has to be loaded first.  src/run.jl relies on the
        # same ordering; the scenario files in index_model/parameters make it
        # explicit.
        include(joinpath(@__DIR__, "..", "population_model.jl"))
        include(joinpath(@__DIR__, "..", "..", "parameters", "no_RA_intermidiate.jl"))
        @test size(p.T) == (2, 2)
        @test p.T[2, 1] ≈ P_LOW
        @test p.T[2, 2] ≈ P_STAY

        q = set_index_T(p; pr = 0.60, re = 1.0)
        @test size(q.T) == (4, 4)
        @test maximum(abs.(vec(sum(q.T, dims = 1)) .- 1)) < 1e-12
        # p_low and p_stay are inherited from the base parameter file
        @test q.T ≈ joint_transition(P_LOW, 0.60, 1.0, P_STAY) atol = 1e-12
        # every other field survives the ComponentArray surgery untouched
        for k in propertynames(p)
            k == :T && continue
            @test getproperty(q, k) == getproperty(p, k)
        end
        # the original is not mutated
        @test size(p.T) == (2, 2)
    end
end
