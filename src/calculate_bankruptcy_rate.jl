using DataFrames
using CSV
using Plots
using Random 

include("../src/model.jl")
include("../src/plotting.jl")

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

    states, m, f, h, y = simulate_trajectory(p,X,5)
    # Plots.scatter(states[1,:], states[2,:])
    b_grid = bmin:((bmax-bmin)/(Nb_grid-1)):bmax
    Δb_grid = Δbmin:((Δbmax-Δbmin)/(NΔb_grid-1)):Δbmax
    # plt = Plots.scatter!(b_grid,[0.0])
    # Plots.scatter!(zeros(length(Δb_grid)),Δb_grid)


    # savings 
    smax = 2.5; smin = -4.0; Ns_grid = parse(Int,Ns_grid)
    s_grid = smin:((smax-smin)/(Ns_grid-1)):smax
    # plt = Plots.scatter(exp.(s_grid),s_grid, label = "Savings", xlabel = "Savings", ylabel = "log savings")

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


    Random.seed!(123)
    N = 2 #100
    T = 2 #100
    times = collect(1:(T+1))
    T_bank = zeros(N)
    T_bank_no_insurance = zeros(N)
    Threads.@threads for i in 1:N
        print(i, " ")
        sim, m, f, h, y = simulate_trajectory(prob.p,Xmc,100)
        x = simulation(prob,El,Xmc,T;with_policy = false, s0 = vcat([1.0,log(-1*prob.p.s̄)],sim[:,end]))
        x_no_insurnace = simulation(prob_no_insurance,El,Xmc,T;with_policy = false, s0 = vcat([1.0,log(-1*prob.p.s̄)],sim[:,end]))
        
        filtered_times = times[x.Ibt.>1.0]
        if length(filtered_times) == 0
            T_bank[i] = T
        else      
            T_bank[i] = filtered_times[1]
        end

        filtered_times_no_insurnace  = times[x_no_insurnace.Ibt.>1.0]
        if length(filtered_times_no_insurnace) == 0
            T_bank_no_insurance[i] = T
        else      
            T_bank_no_insurance[i] = filtered_times_no_insurnace[1]
        end
    end 

    p̂ = sum(T_bank.<T)/(sum(T_bank))
    σ̂p = sqrt(p̂*(1-p̂)/(sum(T_bank)))
    p̂_lower = p̂ - 2*σ̂p
    p̂_upper = p̂ + 2*σ̂p

    p̂_no_I = sum(T_bank_no_insurance.<T)/(sum(T_bank_no_insurance))
    σ̂p_no_I = sqrt(p̂*(1-p̂)/(sum(T_bank_no_insurance)))
    p̂_lower_no_I = p̂ - 2*σ̂p
    p̂_upper_no_I = p̂ + 2*σ̂p

    println(round(p̂,digits=4)," (", round(p̂_lower,digits = 4), ", ", round(p̂_upper,digits = 4), ")")
    println(round(p̂_no_I,digits=4)," (", round(p̂_lower_no_I,digits = 4), ", ", round(p̂_upper_no_I,digits = 4), ")")

    return p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I
end

df =  DataFrame(params=String[],p=Float64[], p_sigma=Float64[], p_lower=Float64[], p_upper=Float64[],
            p_no_I=Float64[], p_no_I_sigma=Float64[], p_no_I_lower=Float64[],  p_no_I_upper= Float64[])

### Cases in table 4 ###
results = calcualte_bankruptcy("no_RA_long")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,["no_RA_long",p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I])

print("ran 1!")
results = calcualte_bankruptcy("no_RA_intermidiate")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,(params="no_RA_intermidiate",p=p̂, p_sigma=σ̂p, p_lower =p̂_lower, p_upper =p̂_upper,
            p_no_I=p̂_no_I, p_no_I_sigma= σ̂p_no_I, p_no_I_lower=p̂_lower_no_I,  p_no_I_upper= p̂_upper_no_I))
print(DataFrame(df))
throw()
results = calcualte_bankruptcy("no_RA_long")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,(params="no_RA_long",p=p̂, p_sigma=σ̂p, p_lower =p̂_lower, p_upper =p̂_upper,
            p_no_I=p̂_no_I, p_no_I_sigma= σ̂p_no_I, p_no_I_lower=p̂_lower_no_I,  p_no_I_upper= p̂_upper_no_I))

results = calcualte_bankruptcy("no_RA_intermidiate_lower_blimit")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,(params="no_RA_intermidiate_lower_blimit",p=p̂, p_sigma=σ̂p, p_lower =p̂_lower, p_upper =p̂_upper,
            p_no_I=p̂_no_I, p_no_I_sigma= σ̂p_no_I, p_no_I_lower=p̂_lower_no_I,  p_no_I_upper= p̂_upper_no_I))
            
### Probabilitieis  table S5 ###
# Short life history
results = calcualte_bankruptcy("no_RA_short_higher_frequency")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,(params="no_RA_short_higher_frequency",p=p̂, p_sigma=σ̂p, p_lower =p̂_lower, p_upper =p̂_upper,
            p_no_I=p̂_no_I, p_no_I_sigma= σ̂p_no_I, p_no_I_lower=p̂_lower_no_I,  p_no_I_upper= p̂_upper_no_I))
            
results = calcualte_bankruptcy("no_RA_short_higher_severity")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,(params="no_RA_short_higher_severity",p=p̂, p_sigma=σ̂p, p_lower =p̂_lower, p_upper =p̂_upper,
            p_no_I=p̂_no_I, p_no_I_sigma= σ̂p_no_I, p_no_I_lower=p̂_lower_no_I,  p_no_I_upper= p̂_upper_no_I))
            
results = calcualte_bankruptcy("no_RA_short_longer_duration")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,(params="no_RA_short_longer_duration",p=p̂, p_sigma=σ̂p, p_lower =p̂_lower, p_upper =p̂_upper,
            p_no_I=p̂_no_I, p_no_I_sigma= σ̂p_no_I, p_no_I_lower=p̂_lower_no_I,  p_no_I_upper= p̂_upper_no_I))
            
# Intermediate life history 
results = calcualte_bankruptcy("no_RA_intermediate_lower_frequency")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,(params="no_RA_intermediate_lower_frequency",p=p̂, p_sigma=σ̂p, p_lower =p̂_lower, p_upper =p̂_upper,
            p_no_I=p̂_no_I, p_no_I_sigma= σ̂p_no_I, p_no_I_lower=p̂_lower_no_I,  p_no_I_upper= p̂_upper_no_I))
            
results = calcualte_bankruptcy("no_RA_intermediate_higher_frequency")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,["no_RA_intermediate_higher_frequency",p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I])

results = calcualte_bankruptcy("no_RA_intermediate_higher_severity")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,["no_RA_intermediate_higher_severity",p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I])

results = calcualte_bankruptcy("no_RA_intermediate_longer_duration")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,["no_RA_intermediate_longer_duration",p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I])

# Long life history 
results = calcualte_bankruptcy("no_RA_long")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,["no_RA_long",p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I])

results = calcualte_bankruptcy("no_RA_long_higher_frequency")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,["no_RA_long_higher_frequency",p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I])

results = calcualte_bankruptcy("no_RA_long_higher_severity")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,["no_RA_long_higher_severity",p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I])

results = calcualte_bankruptcy("no_RA_long_longer_duration")
p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I = results
append!(df,["no_RA_long_longer_duration",p̂, σ̂p, p̂_lower, p̂_upper,p̂_no_I, σ̂p_no_I, p̂_lower_no_I,  p̂_upper_no_I])


df = data.frame(df)
CSV.write(df,"results/bankrptcy_probabilties.csv")