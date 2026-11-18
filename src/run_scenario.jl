#####################################################################
### Solve one scenario of the 3 life histories x {insurance, no insurance}
### x {disaster relief, no relief} comparison.
###
### Usage:
###   julia src/run_scenario.jl <params> <insurance 0|1> <relief 0|1> \
###         Nκ Nη Ns_grid Nb_grid NΔb_grid Nquad maxiter
###
###   insurance = 0 -> premium spending is fixed at zero (η grid = [0.0])
###   relief    = 1 -> model_with_relief.jl (8 discrete states, p.Relief paid
###                    after a declared disaster); 0 -> model.jl (4 states)
###
### Everything else (state grids, κ grid, η upper bound, discount rate,
### insurer parameters, strike and index) is identical across scenarios so
### that value functions and policies are directly comparable.
### Output: results/V2/scenarios/<params>_<ins|noins>_<rel|norel>/
#####################################################################
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
using Plots
params, insurance, relief_flag, Nκ, Nη, Ns_grid, Nb_grid, NΔb_grid, Nquad, maxiter = ARGS
insurance = parse(Int, insurance) == 1
with_relief = parse(Int, relief_flag) == 1

include(joinpath(@__DIR__, with_relief ? "model_with_relief.jl" : "model.jl"))
include(joinpath(@__DIR__, "plotting.jl"))
include(joinpath(@__DIR__, "../parameters/V2/$params.jl"))
include(joinpath(@__DIR__, "../parameters/V2/insurer_capital_model.jl"))
if with_relief
    include(joinpath(@__DIR__, "../parameters/V2/disaster_relief.jl")) # pR
    p = ComponentArray(merge(NamedTuple(p), (pR = pR,)))
    # plotting.jl calls the two argument encode_discrete(Ib, IM); with relief the
    # plots show the slice with no disaster declared (ID = 0).
    encode_discrete(Ib,IM) = encode_discrete(Ib,0,IM)
end

tag = "$(params)_$(insurance ? "ins" : "noins")_$(with_relief ? "rel" : "norel")"
outdir = joinpath(@__DIR__, "..", "results", "V2", "scenarios", tag)
mkpath(outdir)

function main(params, Nκ, Nη, Ns_grid, Nb_grid, NΔb_grid, Nquad, maxiter)
    println("scenario $tag: insurance = $insurance, relief = $with_relief")
    ###################################
    ####  Define random variables #####
    ###################################
    nq = parse(Int, Nquad)
    X, Xmc = with_relief ? init_ranom_variables_relief(p; Nquad = nq) :
                           init_ranom_variables(p; Nquad = nq)

    ######################################
    ####  Define state varible grids #####
    ######################################
    bmax = 0.5; bmin = -3.5; Nb_grid = parse(Int, Nb_grid)
    Δbmax = 1.25; Δbmin = -1.0; NΔb_grid = parse(Int, NΔb_grid)
    if occursin("short", params)
        bmax = 1.0; bmin = -3.0
        Δbmax = 2; Δbmin = -1.25
    elseif occursin("long", params)
        bmax = 0.0; bmin = -4.5
        Δbmax = 0.75; Δbmin = -0.75
    end
    b_grid = bmin:((bmax-bmin)/(Nb_grid-1)):bmax
    Δb_grid = Δbmin:((Δbmax-Δbmin)/(NΔb_grid-1)):Δbmax
    smax = 2.5; smin = -4.0; Ns_grid = parse(Int, Ns_grid)
    s_grid = smin:((smax-smin)/(Ns_grid-1)):smax
    state_grid = [s_grid, b_grid, Δb_grid]

    ###############################
    ####  Define action grid  #####
    ###############################
    κmin = 0.01; κmax = 0.99; Nκ = parse(Int, Nκ)
    κ_grid = κmin:((κmax-κmin)/(Nκ-1)):κmax
    ηmin = 0.0; ηmax = η_max; Nη = parse(Int, Nη)
    η_grid = insurance ? collect(ηmin:((ηmax-ηmin)/(Nη-1)):ηmax) : [0.0]
    action_grid = [κ_grid, η_grid]

    #######################################
    ####  Strike and index functions  #####
    #######################################
    strike(s, p) = 0.5 * (p.b_target + p.b_limit)
    index(s, X) = exp(s[1]) + X[1]

    ##################################
    ####  Set up dynamic program  ####
    ##################################
    δ = p.δ
    prob, El = init(action_grid, state_grid, δ, p, X, index, strike)
    println("solving (maxiter = $maxiter)")
    flush(stdout)
    t0 = time()
    ValueFunctionIterations.solve!(prob, maxiter = parse(Int, maxiter))
    println("solved in $(round(time() - t0, digits = 1)) s")
    save_solution(prob, joinpath(outdir, "sol.jld2"))

    ##################################
    ####  Plot solutions.         ####
    ##################################
    # start solvent, typical mortality (and no disaster declared with relief)
    s0 = with_relief ? [encode_discrete(0,0,0), 0.0, log(0.3), 0.0] :
                       [encode_discrete(0,0),   0.0, log(0.3), 0.0]
    x = simulation(prob, X, Xmc, 250; with_policy = false, s0 = s0)
    savefig(plot_simulation(x, max_USD = 6.0, max_Bt = 1.5, T = 200), joinpath(outdir, "simulation.png"))
    savefig(plot_value_function(prob), joinpath(outdir, "value_function.png"))
    if insurance
        savefig(plot_premiums(prob), joinpath(outdir, "premiums.png"))
        savefig(plot_loading(El, b_grid, Δb_grid, X, prob.p, index, strike), joinpath(outdir, "loading.png"))
    end
    println("wrote scenario outputs to $outdir")
end
main(params, Nκ, Nη, Ns_grid, Nb_grid, NΔb_grid, Nquad, maxiter)
