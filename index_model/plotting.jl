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
    # The axes are recovered from the solved value function itself, so this
    # cannot disagree with the grid the problem was solved on.
    s_grid = sort(unique(prob.V.states[2,:]))
    b_grid = sort(unique(prob.V.states[3,:]))

    # Sliced at Δb = 0 and the typical mortality / no payout state.
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
    Plots.plot!(plt,xlabel = "Stock (biomass)", ylabel = "Value", size = (425,275),legend = :none)
    return plt
end


function simulation(prob,El,Xmc,T;with_policy = false, s0 = nothing)
    p = prob.p
    # Run simulation
    states,actions,rewards,vals,randos=simulate(prob,Xmc,T,with_policy=with_policy, s0 = s0)

    # simulate returns states of size (5, T+1) but actions and randos of size
    # (., T), so every per period series below has length nT = T.  Ibt is the
    # one exception: it keeps all T+1 states, because the bankruptcy indicator
    # of the final state is what the bankruptcy rate calculation reads.
    nT = size(states,2) - 1

    # Get variables and derived quantities
    Ibt = states[1,:]
    st = exp.(states[2,1:nT]).+p.s̄
    bt = exp.(states[3,1:nT])
    bt1 = exp.(states[3,2:end]) # leading Stock
    κt = actions[1,:] # consumption
    ηt = actions[2,:] # premium spending
    Elt = broadcast(t->El(states[:,t],prob.p),1:nT) # Expected losses
    η̄t = broadcast(t->price_per_exposure(Elt[t],p.cf,p.cv), 1:nT) # insurance premiums
    # The index is drawn from the state and shock at t, and is realised entering
    # t+1; Zt[t] pairs with the premium ηt[t] bought at t.
    Zt = broadcast(t -> index(states[:,t],randos[:,t],p), 1:nT) # index
    wt = broadcast(t -> index_insurance(Zt[t],ηt[t],η̄t[t]), 1:nT) # claims
    ft = broadcast(t -> control_rule(bt[t],p.b_target,p.f_target,p.b_limit),1:nT)
    # Mortality applied over period t is carried in state t+1 (update_stock puts
    # the drawn Mt into s_next), so mt[t] pairs with bt[t].
    mt = broadcast(t -> environment_mortaltiy_index(states[5,t+1],p.m̄,p.m̲),1:nT)
    ht = broadcast(t -> baranov(bt[t],mt[t],ft[t]), 1:nT)
    pt = broadcast(t -> randos[2,t], 1:nT)
    c_et = broadcast(t -> randos[3,t], 1:nT)
    yt = broadcast(t -> income(ht[t], ft[t], pt[t], c_et[t]),1:nT)

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

# p is passed in explicitly: it used to be read as a global, which silently
# picked up whatever parameter set happened to be in scope.
function plot_simulation(x,p;T=length(x.wt),t0=1,max_USD = 4.0, max_Bt = 1.5)
    xt = 1:length(x.mt)

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
