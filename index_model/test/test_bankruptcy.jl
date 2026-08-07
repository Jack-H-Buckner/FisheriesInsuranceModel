#####################################################################
#####################################################################
### Tests for demand metric 3, the probability of bankruptcy.
###
### Solves one scenario at a tiny resolution and exercises the whole
### path: paired trajectories, the survival summary, the Poisson
### hazard, the bootstrap, and the CSV the driver writes.
###
### The properties checked here are the ones that would silently
### produce a wrong number rather than an error: reproducibility under
### threading, that the two arms really are paired, and that the
### censoring arithmetic in time_to_bankruptcy is right.
###
### Run from the repository root:
###     julia --project=. --threads=4 index_model/test/test_bankruptcy.jl
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

using Test

ENV["GKSwstype"] = "100"
using Plots, LaTeXStrings

include(joinpath(@__DIR__, "..", "analysis", "bankruptcy.jl"))

const SCEN = "no_RA_int_pr60_re100"
const TINY = (Nκ = 8, Nη = 6, Ns = 5, Nb = 5, NΔb = 4, Nquad = 5)

const WORKDIR = mktempdir()
let q = scenario_parameters(SCEN)
    solve_scenario(q, SCEN; insurance = true, plots = false,
                   outdir = joinpath(WORKDIR, "ins"), TINY...)
    solve_scenario(q, SCEN; insurance = false, plots = false,
                   outdir = joinpath(WORKDIR, "noins"), TINY...)
end

const RI  = rebuild_problem(SCEN; insurance = true,  dir = joinpath(WORKDIR, "ins"))
const RI0 = rebuild_problem(SCEN; insurance = false, dir = joinpath(WORKDIR, "noins"))

# Small enough to run in seconds, long enough that some trajectories fail.
const SMALL = (N = 6, T = 40, B = 200)


@testset "bankruptcy rate" begin

    @testset "time to bankruptcy censoring" begin
        # Ibt has T+1 entries; index 1 is the initial state, so a first hit at
        # index k means k-1 periods elapsed.  Off by one here would bias every
        # hazard rate in the table.
        @test time_to_bankruptcy([1.0, 1.0, 1.0], 2) == (exposure = 2.0, event = false)
        @test time_to_bankruptcy([1.0, 2.0, 2.0], 2) == (exposure = 1.0, event = true)
        @test time_to_bankruptcy([1.0, 1.0, 2.0], 2) == (exposure = 2.0, event = true)
        # a censored trajectory contributes the full window, not its length + 1
        @test time_to_bankruptcy(ones(11), 10).exposure == 10.0
    end

    @testset "starting state" begin
        Xmc = new_shocks(RI.p)
        s0 = stationary_start(RI.p, Xmc, 20; savings = 0.0)
        @test length(s0) == 5
        @test s0[1] == 1.0                     # solvent
        @test exp(s0[2]) ≈ -RI.p.s̄             # savings = 0, full borrowing capacity
        @test s0[5] in (1.0, 2.0)              # mortality state only

        @test exp(stationary_start(RI.p, Xmc, 5; savings = 1.0)[2]) ≈ 1.0 - RI.p.s̄
        # starting at or below the floor would start the firm bankrupt
        @test_throws ArgumentError stationary_start(RI.p, Xmc, 5; savings = RI.p.s̄)
    end

    @testset "paths are reproducible under threading" begin
        a = bankruptcy_paths(RI.prob, RI0.prob; N = SMALL.N, T = SMALL.T,
                             seed = 7, verbose = false)
        b = bankruptcy_paths(RI.prob, RI0.prob; N = SMALL.N, T = SMALL.T,
                             seed = 7, verbose = false)
        # This is what the per trajectory seeding buys, and what a shared
        # MCRandomVariable across threads would break.
        @test a.exposure == b.exposure
        @test a.event    == b.event
        @test a.exposure_no == b.exposure_no
        @test a.s0 == b.s0

        c = bankruptcy_paths(RI.prob, RI0.prob; N = SMALL.N, T = SMALL.T,
                             seed = 8, verbose = false)
        @test c.s0 != a.s0                       # a different seed is a different draw

        @test all(0 .<= a.exposure .<= SMALL.T)
        @test all(0 .<= a.exposure_no .<= SMALL.T)
        # exposure = T exactly when the trajectory survives the window
        @test a.exposure[.!a.event] == fill(float(SMALL.T), count(.!a.event))
        @test all(a.exposure[a.event] .< SMALL.T)
        @test all(a.s0[1, :] .== 1.0)
    end

    @testset "the two arms share a shock path when paired" begin
        # Mortality and payouts evolve independently of savings, so under common
        # random numbers the insured and uninsured runs see an identical
        # biological path and differ only through insurance.  A trajectory that
        # survives with insurance must then have survived at least as long
        # without... no: insurance costs premiums, so that is not guaranteed.
        # What IS guaranteed is that the same starting state is used, and that
        # dropping the pairing changes the uninsured arm but not the insured one.
        p = bankruptcy_paths(RI.prob, RI0.prob; N = SMALL.N, T = SMALL.T,
                             seed = 11, paired = true,  verbose = false)
        u = bankruptcy_paths(RI.prob, RI0.prob; N = SMALL.N, T = SMALL.T,
                             seed = 11, paired = false, verbose = false)
        @test p.s0 == u.s0                        # same starting states
        @test p.exposure == u.exposure            # insured arm runs first, unchanged
    end

    @testset "hazard rates and P10" begin
        paths = bankruptcy_paths(RI.prob, RI0.prob; N = SMALL.N, T = SMALL.T,
                                 seed = 3, verbose = false)
        res = hazard_rates(paths; B = SMALL.B)

        @test res.N == SMALL.N
        @test res.n_events == sum(paths.event)
        @test res.exposure == sum(paths.exposure)
        @test res.frac_bankrupt ≈ mean(paths.event)

        # Poisson MLE, and P10 as its constant hazard transform
        @test res.lambda ≈ sum(paths.event) / sum(paths.exposure)
        @test res.P10 ≈ 1 - exp(-10 * res.lambda)
        @test 0 <= res.P10 <= 1
        @test 0 <= res.P10_no_I <= 1
        # the bootstrap interval brackets the point estimate
        @test res.lambda_lo <= res.lambda <= res.lambda_hi
        @test res.P10_lo <= res.P10 <= res.P10_hi
        @test res.boot_dropped <= res.B

        # the reduction is the difference of the two reported probabilities and,
        # unlike the ratio, is always estimable
        @test res.P10_reduction ≈ res.P10_no_I - res.P10
        @test isfinite(res.P10_reduction_lo)
        @test res.P10_reduction_lo <= res.P10_reduction_hi

        # a different horizon rescales P10 but not the rate or its ratio
        r20 = hazard_rates(paths; B = 50, horizon = 20)
        @test r20.lambda ≈ res.lambda
        @test r20.P10 ≈ 1 - exp(-20 * res.lambda)
        @test isnan(res.lambda_ratio) ? isnan(r20.lambda_ratio) :
              r20.lambda_ratio ≈ res.lambda_ratio
    end

    @testset "degenerate arms are NA, not 0 or Inf" begin
        # No bankruptcies anywhere: the rate is 0, and the ratio is undefined
        # rather than silently 0/0.
        none = (event = falses(4), exposure = fill(50.0, 4),
                event_no = falses(4), exposure_no = fill(50.0, 4))
        r = hazard_rates(none; B = 20)
        @test r.lambda == 0.0
        @test r.P10 == 0.0
        @test isnan(r.lambda_ratio)
        @test r.boot_dropped == 20
        # no bankruptcies either way means no reduction, which IS a number
        @test r.P10_reduction == 0.0

        # Insurance eliminates bankruptcy but the uninsured arm fails: a real
        # ratio of 0, which must NOT be confused with the undefined case.
        one_sided = (event = falses(4), exposure = fill(50.0, 4),
                     event_no = BitVector([true, true, false, false]),
                     exposure_no = [10.0, 20.0, 50.0, 50.0])
        r2 = hazard_rates(one_sided; B = 20)
        @test r2.lambda_ratio == 0.0
        @test r2.P10_no_I ≈ 1 - exp(-10 * 2 / 130)
        @test r2.frac_bankrupt_no_I == 0.5
        # the ratio CI collapses here; the reduction is what carries the result
        @test isnan(r2.ratio_lo)
        @test r2.P10_reduction ≈ r2.P10_no_I
        @test isfinite(r2.P10_reduction_lo) && isfinite(r2.P10_reduction_hi)
    end

    @testset "mismatched solutions are refused" begin
        # two insurance solves, or the wrong base's no insurance solve
        @test_throws ArgumentError bankruptcy_paths(RI.prob, RI.prob;
                                                    N = 2, T = 5, verbose = false)
        @test_throws ArgumentError bankruptcy_paths(RI0.prob, RI0.prob;
                                                    N = 2, T = 5, verbose = false)
    end

    @testset "bankruptcy_rate assembles the reported table" begin
        res = bankruptcy_rate(RI.prob, RI0.prob; N = SMALL.N, T = SMALL.T,
                              B = SMALL.B, verbose = false)
        for k in (:P10, :P10_lo, :P10_hi, :P10_no_I, :lambda, :lambda_ratio,
                  :n_events, :exposure, :horizon, :T, :seed, :paired, :minutes)
            @test haskey(res, k)
        end
        @test res.T == SMALL.T
        @test res.paired
        @test res.minutes >= 0
        @test print_bankruptcy(res; io = devnull) === nothing
    end

    @testset "write_csv" begin
        rows = [(scenario = "a", arm = "iso_premium", P10 = 0.5, n = 1),
                (scenario = "b,c", arm = "full_coverage", P10 = NaN, n = 2)]
        path = joinpath(WORKDIR, "t.csv")
        write_csv(path, rows)
        lines = readlines(path)
        @test lines[1] == "scenario,arm,P10,n"
        @test lines[2] == "a,iso_premium,0.5,1"
        # a comma in a field is quoted, and NaN is written as NA so R and pandas
        # both read it as missing rather than as the string "NaN"
        @test lines[3] == "\"b,c\",full_coverage,NA,2"

        @test_throws ArgumentError write_csv(path, NamedTuple[])
        @test_throws ArgumentError write_csv(path, [(a = 1,), (b = 2,)])
    end
end


#####################################################################
### Merging the per scenario part tables the HPC jobs write
#####################################################################

# Standalone by design - no model, no VFI - so it is included on its own here
# rather than through the analysis stack above.
include(joinpath(@__DIR__, "..", "analysis", "merge_bankruptcy.jl"))

@testset "merge_bankruptcy" begin

    @testset "quoted fields survive the round trip" begin
        # the merge reads back what write_csv wrote, so the two have to agree on
        # quoting or a comma in a field would shift every later column
        @test split_csv_line("a,b,c") == ["a", "b", "c"]
        @test split_csv_line("a,\"b,c\",d") == ["a", "b,c", "d"]
        @test split_csv_line("\"say \"\"hi\"\"\",1") == ["say \"hi\"", "1"]
        @test split_csv_line("a,,b") == ["a", "", "b"]
        @test split_csv_line("a,b,") == ["a", "b", ""]
        @test_throws ArgumentError split_csv_line("a,\"unterminated")

        rows = [(scenario = "a,b", arm = "iso_premium", pr = 1.0)]
        f = joinpath(WORKDIR, "rt.csv")
        write_csv(f, rows)
        @test split_csv_line(readlines(f)[2]) == ["a,b", "iso_premium", "1.0"]
    end

    @testset "parts merge in sweep order" begin
        dir = joinpath(WORKDIR, "parts"); mkpath(dir)
        head = "arm,scenario,cost_tag,pr,P10"
        part(name, lines) = write(joinpath(dir, name), head * "\n" * join(lines, "\n") * "\n")

        # deliberately written out of order, and out of file name order
        part("c.csv", ["full_coverage,s50,fair,0.5,0.3"])
        part("a.csv", ["iso_premium,s60,fair,0.6,0.2"])
        part("b.csv", ["iso_premium,s100,fair,1.0,0.1",
                       "full_coverage,s100,fair,1.0,0.1"])
        # a job that died before writing anything, and one that wrote a header
        # only: reported and skipped, not fatal
        write(joinpath(dir, "d.csv"), "")
        write(joinpath(dir, "e.csv"), head * "\n")

        out = joinpath(WORKDIR, "merged.csv")
        @test main(["--parts=$dir", "--out=$out"]) == 0

        lines = readlines(out)
        @test lines[1] == head
        # arm A down the precision ladder, then arm B down the ladder
        @test [split_csv_line(l)[1:2] for l in lines[2:end]] ==
              [["iso_premium", "s100"], ["iso_premium", "s60"],
               ["full_coverage", "s100"], ["full_coverage", "s50"]]
    end

    @testset "parts from different code versions are refused" begin
        # silently merging mismatched columns would shift values between fields,
        # which is exactly the kind of error that survives into a figure
        dir = joinpath(WORKDIR, "parts_bad"); mkpath(dir)
        write(joinpath(dir, "a.csv"), "arm,scenario,pr\niso_premium,s,1.0\n")
        write(joinpath(dir, "b.csv"), "arm,scenario,pr,P10\niso_premium,s,1.0,0.1\n")
        @test_throws ErrorException main(["--parts=$dir",
                                          "--out=$(joinpath(WORKDIR, "bad.csv"))"])
    end

    @testset "nothing to merge is reported, not silently empty" begin
        @test main(["--parts=$(joinpath(WORKDIR, "does_not_exist"))"]) == 1
        empty_dir = mktempdir()
        @test main(["--parts=$empty_dir"]) == 1
        write(joinpath(empty_dir, "a.csv"), "")
        @test main(["--parts=$empty_dir"]) == 1
    end
end
