#####################################################################
#####################################################################
### This file creates plotting functions for the fisheries insurnace 
### demand model.
###
### The full model is defined in model.jl. I have chosen to split the 
### code in this way to keep individual files from being too long.  
### Jack H. Buckner, Oregon State University, 2025
#####################################################################
#####################################################################

using LaTeXStrings, CairoMakie

function plot_value_function(prob)
    s_grid,b_grid,Δb_grid,m_grid = prob.V[1].grid
    Mt = 1.0; Ib = 1.0
    V = broadcast(b -> prob.V([Ib,s_grid[1],b,0.0,Mt]), b_grid)
    colors = cgrad(:viridis, length(s_grid))
    plt = Plots.plot(b_grid,V, label = string("s = ",s_grid[1]), color = colors[1], linewidth = 2)
    i = 0
    for s in s_grid[2:end]
        i +=1
        V = broadcast(b -> prob.V([Ib,s,b,0.0,Mt]), b_grid)
        plt = Plots.plot!(b_grid,V, label = string("s = ",s), color = colors[round(Int,265*i/length(s_grid))], linewidth = 2)
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
    Elt = broadcast(t->El(states[3,t],states[4,t],states[5,t]),1:(size(states)[2]-1)) # Expected losses 
    η̄t = broadcast(t->price_per_exposure(Elt[t],p.cf,p.cv), 1:(size(states)[2]-1)) # insurance premiums
    Zϵt = broadcast(t -> randos[1,t], 1:(size(states)[2]-1))
    Zt= broadcast(t -> index(states[3:5,t],0.0), 2:size(states)[2]) # index 
    Zct = broadcast(t -> strike(states[3:5,t],p),1:(size(states)[2]-1)) # claims strike
    wt = broadcast(t -> index_insurance(Zt[t],Zct[t],ηt[t],η̄t[t]), 1:(size(states)[2]-2)) # claims
    # wt = Zt .< Zct
    # wt3 = ηt .> 0
    ft = broadcast(t -> control_rule(bt[t],p.b_target,p.f_target,p.b_limit),1:(size(states)[2]-1))
    mt = broadcast(t -> environment_mortaltiy(states[5,t],p.m̄,p.m̲),2:size(states)[2])
    ht = broadcast(t -> baranov(bt[t],mt[t],ft[t]), 1:(size(states)[2]-1))
    pt = broadcast(t -> randos[2,t], 1:(size(states)[2]-1))
    c_et = broadcast(t -> randos[3,t], 1:(size(states)[2]-1))
    yt = broadcast(t -> income(ht[t], ft[t], pt[t], c_et[t]),1:(size(states)[2]-1))

    wt2 = st[2:end] .-update_savings.(st[1:(end-1)],κt[1:(end-1)], ηt[1:(end-1)], yt[1:(end-1)], 0.0, prob.p.y0, prob.p.k)

    return (Ibt = Ibt, st=st, bt=bt, bt1=bt1, κt=κt, ηt=ηt, Elt=Elt, η̄t=η̄t, Zt=Zt, Zct=Zct, wt=wt, ft=ft, mt=mt, ht=ht, pt=pt, c_et=c_et, yt=yt)
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

function plot_premiums(prob;Δblow = -0.5, Δbhigh = 0.5)
    s_grid = prob.V.Bslines[1].grid[1]
    b_grid = prob.V.Bslines[1].grid[2]
    m_Δb_grid = [Δblow  Δblow Δbhigh Δbhigh; 1 2 1 2]
    P =zeros(length(s_grid),length(b_grid),size(m_Δb_grid)[2])
    i = 0
    for s in s_grid 
        i +=1; j = 0
        for b in b_grid
            j += 1
            for k in 1:size(m_Δb_grid)[2]
                P[i,j,k] = prob.P([1,s,b,m_Δb_grid[1,k],m_Δb_grid[2,k]])[2]
            end
        end
    end

    b_grid = exp.(b_grid); s_grid = exp.(s_grid).+p.s̄
    p1 = Plots.heatmap(b_grid,s_grid,P[:,:,1],  xlabel = "Stock", ylabel = "Savings", color = :viridis, title = "low ΔB; typical M")
    p2 = Plots.heatmap(b_grid,s_grid,P[:,:,2],  xlabel = "Stock", ylabel = "Savings", color = :viridis, title = "low ΔB; high M")
    p3 = Plots.heatmap(b_grid,s_grid,P[:,:,3],  xlabel = "Stock", ylabel = "Savings", color = :viridis, title = "high ΔB; typical M")
    p4 = Plots.heatmap(b_grid,s_grid,P[:,:,4],  xlabel = "Stock", ylabel = "Savings", color = :viridis, title = "high ΔB; high M")
    
    return Plots.plot(p1,p2,p3,p4)
end 

function plot_value_function(prob;Δblow = -0.5, Δbhigh = 0.5)
    s_grid = prob.V.Bslines[1].grid[1]
    b_grid = prob.V.Bslines[1].grid[2]
    m_Δb_grid = [Δblow  Δblow Δbhigh Δbhigh; 1 2 1 2]
    titles = ["low ΔB; typical M", "low ΔB; high M", "high ΔB; typical M", "high ΔB; high M"]
    colors = cgrad(:viridis, length(s_grid))
    i = 0;plts = []
    for i in 1:size(m_Δb_grid)[2]
        j = 0; mt = m_Δb_grid[1,i]; Δbt = m_Δb_grid[2,i] 
        V = broadcast(b -> prob.V([1.0,s_grid[1],b,Δbt,mt]), b_grid)
        plt = Plots.plot(exp.(b_grid),V,color = colors[1], linewidth = 2, label = "", 
                        xlabel = "Stock", ylabel = "Utility present value", title = titles[i])
        for s in s_grid
            j +=1
            V = broadcast(b -> prob.V([1.0,s,b,Δbt,mt]), b_grid)
            Plots.plot!(plt,exp.(b_grid),V, color = colors[round(Int,256*j/length(s_grid))],linewidth = 2,label = "")
        end 
        push!(plts,plt)
    end
    Plots.plot(plts...,layout = (2,2))
end 


function plot_policy(coverage_1,coverage_2,coverage_3,title;
        titles = ("Fast life history", "Intermediate life history", "Slow life history"))

    η_1, η̄_1, claims_1, state_1 = coverage_1
    η_2, η̄_2, claims_2, state_2 = coverage_2
    η_3, η̄_3, claims_3, state_3 = coverage_3

    coverage_1 = η_1./ η̄_1
    coverage_2 = η_2./ η̄_2
    coverage_3 = η_3./ η̄_3

    f = Figure(size = (800,450))
    ax1 = CairoMakie.Axis(f[1, 1], title = titles[1], ylabel = "Savings",
                xticklabelsize = 0, xticksize = 0, titlesize = 18,ylabelsize = 18,yticklabelsize = 16)
    ax2 = CairoMakie.Axis(f[1, 2], title = titles[2],
                xticklabelsize = 0, xticksize = 0,yticklabelsize=0,yticksize=0, titlesize = 18)
    ax3 = CairoMakie.Axis(f[1, 3], title = titles[3],
                xticklabelsize = 0, xticksize = 0,yticklabelsize=0,yticksize=0, titlesize = 18)

    centers_x = collect(exp.(state_1[1]))
    centers_y = collect(exp.(state_1[2]))
    ind = 2

    data1 = η_1[:,:,2,1].+0.01
    data2 = η_2[:,:,2,1].+0.01
    data3 = η_3[:,:,2,1].+0.01

    scale = ReversibleScale(x -> x, x -> x)
    hm = CairoMakie.heatmap!(ax1, centers_y,  centers_x, data1';colorscale = scale,colorrange = (0.0,0.3))
    hm = CairoMakie.heatmap!(ax2, centers_y,  centers_x, data2';colorscale = scale,colorrange = (0.0,0.3))
    hm = CairoMakie.heatmap!(ax3, centers_y,  centers_x, data3';colorscale = scale,colorrange = (0.0,0.3))
    CairoMakie.Colorbar(f[1, 4], hm,  label = "Premiums", ticks = [0.0, 0.1, 0.2, 0.3], labelsize = 18, ticklabelsize = 16)


    ax1 = CairoMakie.Axis(f[2, 1], xlabel = "Stock", ylabel = "Savings",
                ylabelsize = 18,xlabelsize = 18,xticklabelsize = 16,yticklabelsize = 16)
    ax2 = CairoMakie.Axis(f[2, 2],  xlabel = "Stock",
                yticklabelsize=0,yticksize=0,ylabelsize = 18,xlabelsize = 18,xticklabelsize = 16)
    ax3 = CairoMakie.Axis(f[2, 3],  xlabel = "Stock",
                yticklabelsize=0,yticksize=0,ylabelsize = 18,xlabelsize = 18,xticklabelsize = 16)

    ind = 2

    data1 = coverage_1[:,:,ind ,1].*claims_1[:,:,ind ,1].+0.1
    data2 = coverage_2[:,:,ind ,1].*claims_2[:,:,ind ,1].+0.1
    data3 = coverage_3[:,:,ind ,1].*claims_3[:,:,ind ,1].+0.1

    hm = CairoMakie.heatmap!(ax1, centers_y,  centers_x, data1';colorscale = scale,colorrange = (0.0,10.0))
    hm = CairoMakie.heatmap!(ax2, centers_y,  centers_x, data2';colorscale = scale,colorrange = (0.0,10.0))
    hm = CairoMakie.heatmap!(ax3, centers_y,  centers_x, data3';colorscale = scale,colorrange = (0.0,10.0))
    CairoMakie.Colorbar(f[2, 4], hm,  label = "Claims", labelsize = 18, ticklabelsize = 16)
    #CairoMakie.save("../../../documents/policy_figure_$title.png",f)
    return f
end 