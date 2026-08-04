include("../src/model.jl"); include("../src/plotting.jl")
using Plots
using Random
using Statistics
using DataFrames, CSV


function calcualte_bankruptcy(pars)
    
    Nκ = "35"; Nη = "25"; Ns_grid = "24"; Nb_grid = "24"; NΔb_grid = "18"; Nquad = "50"
    include("../parameters/$pars.jl")

    ###################################
    ####  Define random variables #####
    ###################################
    X, Xmc = init_ranom_variables(p; Nquad = parse(Int,Nquad))

    # define state variables and check the grid 
    ######################################
    ####  Define state varible grids #####
    ######################################
    # check paramters to get min and max values
    bmax = 0.5; bmin = -3.5; Nb_grid = parse(Int,Nb_grid)
    Δbmax = 1.25; Δbmin = -1.0; NΔb_grid = parse(Int,NΔb_grid)
    if occursin("short", pars)
        bmax = 1.0; bmin = -3.0
        Δbmax = 2; Δbmin = -1.25
    elseif occursin("long", pars)
        bmax = 0.0; bmin = -4.5
        Δbmax = 0.75; Δbmin = -0.75
    end 
    # stock 

    states, m, f, h, y = simulate_trajectory(p,X,5000)
    Plots.scatter(states[1,:], states[2,:])
    b_grid = bmin:((bmax-bmin)/(Nb_grid-1)):bmax
    Δb_grid = Δbmin:((Δbmax-Δbmin)/(NΔb_grid-1)):Δbmax
    plt = Plots.scatter!(b_grid,[0.0])
    Plots.scatter!(zeros(length(Δb_grid)),Δb_grid)


    # savings 
    smax = 2.5; smin = -4.0; Ns_grid = parse(Int,Ns_grid)
    s_grid = smin:((smax-smin)/(Ns_grid-1)):smax
    plt = Plots.scatter(exp.(s_grid),s_grid, label = "Savings", xlabel = "Savings", ylabel = "log savings")

    state_grid = [s_grid,b_grid,Δb_grid]
    ###############################
    ####  Define action grid  #####
    ###############################
    # no inusruacne
    κmin = 0.01; κmax = 0.99; Nκ = parse(Int,Nκ)
    κ_grid = κmin:((κmax-κmin)/(Nκ-1)):κmax
    ηmin = 0.0; ηmax=1e-10
    η_grid = ηmin:(ηmax-ηmin):ηmax
    no_insurance_action_grid = [κ_grid,η_grid]


    # no inusruacne
    κmin = 0.01; κmax = 0.99
    κ_grid = κmin:((κmax-κmin)/(Nκ-1)):κmax
    ηmin = 0.0; ηmax=0.2; Nη = parse(Int,Nη)
    η_grid = ηmin:((ηmax-ηmin)/(Nη-1)):ηmax
    action_grid = [κ_grid,η_grid]

    #######################################
    ####  Strike and index functions  #####
    #######################################
    function strike(s,p)
        return 0.5*(p.b_target +p.b_limit) 
    end

    function index(s,X)
        return exp(s[1])+X[1]
    end 

    ##################################
    ####  Set up dynamic program  ####
    ##################################
    δ = 0.95
    prob, El = init(action_grid,state_grid,δ,p,X,index,strike)
    prob_no_insurance, El = init(no_insurance_action_grid,state_grid,δ,p,X,index,strike)

    ##################################
    ####  load in solutions       ####
    ##################################
    load_solution(prob,"results/$pars/sol.jld2")
    load_solution(prob_no_insurance,"results/no_insurance_$pars/sol.jld2")


    N = 60
    T = 250
    times = collect(1:(T+1))
    # Per-trajectory exposure time (periods observed until bankruptcy or
    # censoring at T) and event indicator (1 = bankrupt within window).
    exposure              = zeros(N)
    exposure_no_insurance = zeros(N)
    event              = falses(N)
    event_no_insurance = falses(N)
    # Per-iteration RNG seeding. Julia's default RNG is thread-local, so a
    # single Random.seed! before a Threads.@threads loop only affects one
    # thread — results would depend on JULIA_NUM_THREADS and scheduling.
    # Seeding inside the loop makes each i reproducible regardless of
    # threading. Hashing (123, i) avoids correlations between adjacent seeds.
    Threads.@threads for i in 1:N
        Random.seed!(hash((123, i)) % UInt32)
        print(i, " ")
        sim, m, f, h, y = simulate_trajectory(prob.p,Xmc,T)
        x = simulation(prob,El,Xmc,T;with_policy = false, s0 = vcat([1.0,log(-1*prob.p.s̄)],sim[:,end]))
        x_no_insurnace = simulation(prob_no_insurance,El,Xmc,T;with_policy = false, s0 = vcat([1.0,log(-1*prob.p.s̄)],sim[:,end]))

        # times[k] = k indexes the state vector x.Ibt of length T+1, where
        # index 1 is the initial state (t=0, never bankrupt). An event first
        # detected at index k means k-1 periods of exposure elapsed before
        # bankruptcy. With no event, the trajectory is censored at T.
        filtered_times = times[x.Ibt .> 1.0]
        if isempty(filtered_times)
            exposure[i] = T
            event[i]    = false
        else
            exposure[i] = filtered_times[1] - 1
            event[i]    = true
        end

        filtered_times_no_insurnace = times[x_no_insurnace.Ibt .> 1.0]
        if isempty(filtered_times_no_insurnace)
            exposure_no_insurance[i] = T
            event_no_insurance[i]    = false
        else
            exposure_no_insurance[i] = filtered_times_no_insurnace[1] - 1
            event_no_insurance[i]    = true
        end
    end

    # Poisson MLE hazard rate λ̂ = (# events) / (total person-time). Under
    # a constant-hazard / exponential time-to-bankruptcy assumption this is
    # both the MLE and an unbiased rate estimator.
    n_events      = sum(event)
    n_events_no_I = sum(event_no_insurance)
    λ̂      = n_events      / sum(exposure)
    λ̂_no_I = n_events_no_I / sum(exposure_no_insurance)
    ratio  = λ̂ / λ̂_no_I

    # Paired nonparametric bootstrap over the N trajectories. Pairing
    # preserves the within-iteration correlation between matched
    # with/without-insurance runs (both use the same Xmc per i), so the
    # ratio CI does not over-state the uncertainty.
    B = 2000
    Random.seed!(456)
    λ_boot         = zeros(B)
    λ_no_boot      = zeros(B)
    logratio_boot  = zeros(B)
    for b in 1:B
        idx = rand(1:N, N)
        λb    = sum(event[idx])              / sum(exposure[idx])
        λb_no = sum(event_no_insurance[idx]) / sum(exposure_no_insurance[idx])
        λ_boot[b]        = λb
        λ_no_boot[b]     = λb_no
        logratio_boot[b] = log(λb / λb_no)
    end

    # 95% percentile CIs. With few events per case, bootstrap resamples can
    # contain zero events (λ = 0) or zero no-insurance events, producing
    # 0 / Inf / NaN log-ratios — drop those before quantiling. Log-scale CI
    # for the ratio is back-transformed via exp.
    λ̂_lower,      λ̂_upper      = quantile(λ_boot,    [0.025, 0.975])
    λ̂_lower_no_I, λ̂_upper_no_I = quantile(λ_no_boot, [0.025, 0.975])
    logratio_clean = filter(isfinite, logratio_boot)
    ratio_lower, ratio_upper =
        isempty(logratio_clean) ? (NaN, NaN) :
        exp.(quantile(logratio_clean, [0.025, 0.975]))

    # 10-year bankruptcy probability under the constant-hazard model:
    # P10 = 1 - exp(-10 λ). Monotone in λ, so the bootstrap CI endpoints on
    # λ transform directly to a valid CI on P10. P10 is reported as the
    # interpretable marginal probability; ratios remain on the instantaneous
    # λ scale (rate ratios are scale-invariant in horizon, P10 ratios are not).
    horizon = 10
    P10            = 1 - exp(-horizon * λ̂)
    P10_lower      = 1 - exp(-horizon * λ̂_lower)
    P10_upper      = 1 - exp(-horizon * λ̂_upper)
    P10_no_I       = 1 - exp(-horizon * λ̂_no_I)
    P10_lower_no_I = 1 - exp(-horizon * λ̂_lower_no_I)
    P10_upper_no_I = 1 - exp(-horizon * λ̂_upper_no_I)

    println("events (with / no insurance): ", n_events, " / ", n_events_no_I, " out of ", N)
    println("P10      = ", round(P10,     digits=4), " (", round(P10_lower,     digits=4), ", ", round(P10_upper,     digits=4), ")")
    println("P10_no_I = ", round(P10_no_I,digits=4), " (", round(P10_lower_no_I,digits=4), ", ", round(P10_upper_no_I,digits=4), ")")
    println("λ ratio  = ", round(ratio, digits=4),
            " (", round(ratio_lower, digits=4), ", ", round(ratio_upper, digits=4), ")")

    return P10, P10_lower, P10_upper,
           P10_no_I, P10_lower_no_I, P10_upper_no_I,
           ratio, ratio_lower, ratio_upper,
           n_events, n_events_no_I
end


cases = ["no_RA_short", "no_RA_intermidiate", "no_RA_long"]

results = DataFrame(
    case                   = String[],
    n_events_with          = Int[],
    n_events_no_insurance  = Int[],
    P10_with_insurance     = Float64[],
    P10_with_lower         = Float64[],
    P10_with_upper         = Float64[],
    P10_no_insurance       = Float64[],
    P10_no_lower           = Float64[],
    P10_no_upper           = Float64[],
    lambda_ratio_with_over_no = Float64[],
    lambda_ratio_lower     = Float64[],
    lambda_ratio_upper     = Float64[],
)

for case in cases
    P10, P10_lower, P10_upper,
    P10_no_I, P10_lower_no_I, P10_upper_no_I,
    ratio, r_lower, r_upper,
    n_events, n_events_no_I = calcualte_bankruptcy(case)

    push!(results, (case,
        n_events, n_events_no_I,
        P10, P10_lower, P10_upper,
        P10_no_I, P10_lower_no_I, P10_upper_no_I,
        ratio, r_lower, r_upper))
end

CSV.write("manuscript/table_4.csv", results)
results