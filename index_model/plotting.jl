#####################################################################
#####################################################################
### This file creates plotting functions for the fisheries insurnace 
### demand model.
###
### The full model is defined in index_model/model.jl. I have chosen 
### to split the code in this way to keep individual files from being 
### too long.  
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################



function plot_value_function(prob)
    s_grid,b_grid,Δb_grid,m_grid = prob.V.states
    I_grid = unique(prob.V.states[1,:])
    s_grid = sort(unique(prob.V.states[2,:]))
    b_grid = sort(unique(prob.V.states[3,:]))
    Δb_grid = unique(prob.V.states[4,:])
    m_grid = unique(prob.V.states[5,:])

    Mt = 1.0; Ib = 1.0
    V = broadcast(b -> prob.V([Ib,s_grid[1],b,0.0,Mt]), b_grid)
    colors = cgrad(:viridis, length(s_grid))
    plt = Plots.plot(exp.(b_grid),V, label = string("s = ",s_grid[1]), color = colors[1], linewidth = 2)
    i = 0
    for s in s_grid[2:end]
        i +=1
        V = broadcast(b -> prob.V([Ib,s,b,0.0,Mt]), b_grid)
        plt = Plots.plot!(exp.(b_grid),V, label = string("s = ",s), color = colors[round(Int,256*i/length(s_grid))], linewidth = 2)
    end 
    Plots.plot!(plt,xlabel = "Stock (biomass)", ylabel = "Premiums", size = (425,275),legend = :none)
    return plt
end


function simulation(prob,El,Xmc,T;with_policy = false, s0 = nothing)
    p = prob.p
    # Run simulation 
    states,actions,rewards,vals,randos=simulate(prob,Xmc,T,with_policy=with_policy, s0 = s0)

    # Get variables and derived quantities
    Ibt = states[1,:]
    st = exp.(states[2,1:(end-1)]).+p.s̄
    bt = exp.(states[3,1:(end-1)]) 
    bt1 = exp.(states[3,2:end]) # leading Stock
    κt = actions[1,:] # consumption
    ηt = actions[2,:] # premium spending
    Elt = broadcast(t->El(states[:,t],prob.p),1:(size(states)[2]-1)) # Expected losses 
    η̄t = broadcast(t->price_per_exposure(Elt[t],p.cf,p.cv), 1:(size(states)[2]-1)) # insurance premiums
    Zϵt = broadcast(t -> randos[1,t], 1:(size(states)[2]-1))
    Zt= broadcast(t -> index(states[:,t-1],randos[:,t-1],p), 2:size(states)[2]) # index 
    #Zct = broadcast(t -> strike(states[3:5,t],p),1:(size(states)[2]-1)) # claims strike
    wt = broadcast(t -> index_insurance(Zt[t],ηt[t],η̄t[t]), 1:(size(states)[2]-2)) # claims
    # wt = Zt .< Zct
    # wt3 = ηt .> 0
    ft = broadcast(t -> control_rule(bt[t],p.b_target,p.f_target,p.b_limit),1:(size(states)[2]-1))
    mt = broadcast(t -> environment_mortaltiy_index(states[5,t],p.m̄,p.m̲),2:size(states)[2])
    ht = broadcast(t -> baranov(bt[t],mt[t],ft[t]), 1:(size(states)[2]-1))
    pt = broadcast(t -> randos[2,t], 1:(size(states)[2]-1))
    c_et = broadcast(t -> randos[3,t], 1:(size(states)[2]-1))
    yt = broadcast(t -> income(ht[t], ft[t], pt[t], c_et[t]),1:(size(states)[2]-1))

    wt2 = st[2:end] .-update_savings.(st[1:(end-1)],κt[1:(end-1)], ηt[1:(end-1)], yt[1:(end-1)], 0.0, prob.p.y0, prob.p.k)

    return (Ibt = Ibt, st=st, bt=bt, bt1=bt1, κt=κt, ηt=ηt, Elt=Elt, η̄t=η̄t, Zt=Zt,  wt=wt, ft=ft, mt=mt, ht=ht, pt=pt, c_et=c_et, yt=yt)
end


function mask_to_vspans(x::AbstractVector, mask::AbstractVector{<:Integer})
    @assert length(x) == length(mask) "x and mask must be the same length"
    spans = []
    n = length(mask)
    i = 1
    while i <= n
        if mask[i] == 1
            j = i
            while j <= n && mask[j] == 1
                j += 1
            end
            # run of ones is i:(j-1)
            push!(spans, [x[i], x[j-1]+1.0])
            i = j
        else
            i += 1
        end
    end
    return spans
end

function plot_simulation(x;T=length(x.wt),t0=1,max_USD = 4.0, max_Bt = 1.5)
    xt = 1:length(x.mt)
    spans1 = mask_to_vspans(xt[t0:T], x.mt[t0:T] .> 0.5)
    spans2 = mask_to_vspans(xt[t0:T], x.bt[t0:T] .< 0.3)

    p1 = Plots.scatter(xt[t0:T],x.bt[t0:T], color = :black, label = "Stock       ", ylabel = "Biomass")
    Plots.scatter!(xt[t0:T],max_Bt *(x.mt[t0:T].>0.5), color = :red, label = string("High ",L"M_t"), alpha = 0.75)
    Plots.plot!(xt[t0:T],x.bt[t0:T], label = "", color = :black, ylims = (-0.05,1.05*max_Bt))
    Plots.hline!([0.3], label = "Threshold", color = :grey, yticks = [0.0,0.2,0.4,0.6,0.8,1.0])

    p2 = Plots.scatter(xt[t0:T],x.wt[t0:T], label = "Claims", ylims = (-0.1,max_USD), color = :black, ylabel = "USD")
    Plots.plot!(xt[t0:T],x.wt[t0:T], label = "", color = :black)
    Plots.scatter!(xt[t0:T],x.yt[t0:T],  color = :grey, label = "Income   ", ylabel = "USD", alpha = 0.5)
    Plots.plot!(xt[t0:T],x.yt[t0:T], label = "",  color = :grey)

    p3 = Plots.scatter(xt[t0:T],x.st[t0:T].-p.s̄.-x.ηt[t0:T].-x.κt[t0:T], label = string("Savings"), color = :black, ylabel = "USD")
    Plots.plot!(xt[t0:T],x.st[t0:T].-p.s̄.-x.ηt[t0:T].-x.κt[t0:T], label = "", color = :black)
    Plots.scatter!(xt[t0:T],x.ηt[t0:T], label = "Premiums", ylims = (-0.1,max_USD), color = :grey, alpha = 0.5)
    Plots.plot!(xt[t0:T],x.ηt[t0:T], label = "",  color = :grey, xlabel = "Time")


    Plots.plot(p1,p2,p3, layout = (3,1), markersize = 3.0, linewidth = 0.75, size = (600,400), 
        legend=:outerright,legend_foreground_color=:transparent,legendfontsize = 10,
        tickfontsize = 10, guidefontsize = 14)

end
