
using Pkg, Random, Plots, LaTeXStrings
using Plots.Measures: mm

include("../src/model.jl")
include("../src/plotting.jl")

# Define state grid to set up problem
bmax = 0.25; bmin = -3.5; Nb_grid = 12
b_grid = bmin:((bmax-bmin)/(Nb_grid-1)):bmax
Δbmax = 1.25; Δbmin = -1.0; NΔb_grid = 7
Δb_grid = Δbmin:((Δbmax-Δbmin)/(NΔb_grid-1)):Δbmax
smax = 2.5; smin = -4.0; Ns_grid = 10
s_grid = smin:((smax-smin)/(Ns_grid-1)):smax
state_grid = [s_grid,b_grid,Δb_grid]

# Define action grid
κmin = 0.01; κmax = 0.99; Nκ = 30
κ_grid = κmin:((κmax-κmin)/(Nκ-1)):κmax
ηmin = 0.0; ηmax=0.2; Nη = 30
η_grid = ηmin:((ηmax-ηmin)/(Nη-1)):ηmax
action_grid = [κ_grid,η_grid]

# Define strike and index functions
function strike(s,p)
    return 0.5*(p.b_target +p.b_limit)
end

function index(s,X)
    return exp(s[1])+X[1]
end

δ = 0.95
seed = 232
T_sim = 101
s0 = [1.0, 0.0, log(0.3), 0.0, 1.0]

# Bankruptcy-utility ("no_RA") intermediate life history under three
# fixed-cost regimes. y0 is the exogenous income offset, so lower y0 = higher
# effective fixed costs and vice versa.
pars = "no_RA_intermidiate_lower_y0"
include("../parameters/$pars.jl")
X, Xmc = init_ranom_variables(p; Nquad = 15)
prob, El = init(action_grid, state_grid, δ, p, X, index, strike)
load_solution(prob, "results/$pars/sol.jld2")
Random.seed!(seed)
x_A = simulation(prob, El, Xmc, T_sim; with_policy = false, s0 = s0)
s̄_A = p.s̄

pars = "no_RA_intermidiate"
include("../parameters/$pars.jl")
X, Xmc = init_ranom_variables(p; Nquad = 15)
prob, El = init(action_grid, state_grid, δ, p, X, index, strike)
load_solution(prob, "results/$pars/sol.jld2")
Random.seed!(seed)
x_B = simulation(prob, El, Xmc, T_sim; with_policy = false, s0 = s0)
s̄_B = p.s̄

pars = "no_RA_intermidiate_higher_y0"
include("../parameters/$pars.jl")
X, Xmc = init_ranom_variables(p; Nquad = 15)
prob, El = init(action_grid, state_grid, δ, p, X, index, strike)
load_solution(prob, "results/$pars/sol.jld2")
Random.seed!(seed)
x_C = simulation(prob, El, Xmc, T_sim; with_policy = false, s0 = s0)
s̄_C = p.s̄

# Build a 3×3 grid: columns = fixed-cost regime (A, B, C); rows = biomass, USD, USD.
T = 100; t0 = 1

panels = [
    (label = "A", title = "Higher fixed costs", x = x_A, s̄ = s̄_A),
    (label = "B", title = "Base fixed costs",   x = x_B, s̄ = s̄_B),
    (label = "C", title = "Lower fixed costs",  x = x_C, s̄ = s̄_C),
]

# Compute global y-axis limits per row from the data so no points are clipped,
# while keeping axes comparable across columns.
pad(lo, hi) = (lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))

bt_lo = minimum(minimum(p.x.bt[t0:T]) for p in panels)
bt_hi = maximum(maximum(p.x.bt[t0:T]) for p in panels)
max_Bt = max(1.0, bt_hi)                       # keep ≥ 1 so 0.3 threshold sits inside axis
ylim_row1 = (min(-0.05, bt_lo - 0.05), 1.05 * max_Bt)

row2_lo = minimum(min(minimum(p.x.wt[t0:T]), minimum(p.x.yt[t0:T])) for p in panels)
row2_hi = maximum(max(maximum(p.x.wt[t0:T]), maximum(p.x.yt[t0:T])) for p in panels)
ylim_row2 = pad(min(-0.1, row2_lo), row2_hi)

savings_all = [p.x.st[t0:T] .- p.s̄ .- p.x.ηt[t0:T] .- p.x.κt[t0:T] for p in panels]
row3_lo = min(minimum(minimum(s) for s in savings_all),
              minimum(minimum(p.x.ηt[t0:T]) for p in panels))
row3_hi = max(maximum(maximum(s) for s in savings_all),
              maximum(maximum(p.x.ηt[t0:T]) for p in panels))
ylim_row3 = pad(min(-0.1, row3_lo), row3_hi)

row1 = []; row2 = []; row3 = []

for (i, panel) in enumerate(panels)
    x  = panel.x
    s̄ = panel.s̄
    xt = 1:length(x.mt)

    is_left  = i == 1
    is_right = i == 3
    legend_pos = is_right ? :outerright : :none
    col_title = string(panel.label, ") ", panel.title)

    # Drop data and shade x-ranges where Ibt > 1 (firm is bankrupt).
    bankrupt_mask = x.Ibt[t0:T] .> 1.0
    bankrupt_spans = mask_to_vspans(xt[t0:T], bankrupt_mask)
    bankrupt_xs = Float64[]
    for s in bankrupt_spans
        push!(bankrupt_xs, s[1], s[2])
    end
    mask_nan(v) = [bankrupt_mask[k] ? NaN : float(v[k]) for k in eachindex(v)]

    bt_p      = mask_nan(x.bt[t0:T])
    mt_flag_p = mask_nan(max_Bt .* (x.mt[t0:T] .> 0.5))
    wt_p      = mask_nan(x.wt[t0:T])
    yt_p      = mask_nan(x.yt[t0:T])
    savings_p = mask_nan(x.st[t0:T] .- s̄ .- x.ηt[t0:T] .- x.κt[t0:T])
    ηt_p      = mask_nan(x.ηt[t0:T])

    # Initialize each plot empty so the bankruptcy vspan sits behind the data.
    xlim = (xt[t0], xt[T] + 1)

    # Row 1: stock biomass + high mortality flag
    p1 = Plots.plot(; title = col_title, legend = legend_pos,
        ylabel = is_left ? "Biomass" : "",
        xlims = xlim, ylims = ylim_row1)
    if !isempty(bankrupt_xs)
        Plots.vspan!(p1, bankrupt_xs, color = :grey, alpha = 0.25,
            linewidth = 0, label = is_right ? "Bankrupt" : "")
    end
    Plots.scatter!(p1, xt[t0:T], bt_p, color = :black,
        label = is_right ? "Stock       " : "")
    Plots.scatter!(p1, xt[t0:T], mt_flag_p, color = :red,
        label = is_right ? string("High ", L"M_t") : "", alpha = 0.75)
    Plots.plot!(p1, xt[t0:T], bt_p, label = "", color = :black)
    Plots.hline!(p1, [0.3], label = is_right ? "Threshold" : "",
        color = :grey)

    # Row 2: claims and income (USD)
    p2 = Plots.plot(; legend = legend_pos,
        ylabel = is_left ? "USD" : "",
        xlims = xlim, ylims = ylim_row2)
    if !isempty(bankrupt_xs)
        Plots.vspan!(p2, bankrupt_xs, color = :grey, alpha = 0.25,
            linewidth = 0, label = is_right ? "Bankrupt" : "")
    end
    Plots.scatter!(p2, xt[t0:T], wt_p, color = :black,
        label = is_right ? "Claims" : "")
    Plots.plot!(p2, xt[t0:T], wt_p, label = "", color = :black)
    Plots.scatter!(p2, xt[t0:T], yt_p, color = :grey,
        label = is_right ? "Income   " : "", alpha = 0.5)
    Plots.plot!(p2, xt[t0:T], yt_p, label = "", color = :grey)

    # Row 3: savings and premiums (USD)
    p3 = Plots.plot(; legend = legend_pos,
        ylabel = is_left ? "USD" : "", xlabel = "Time",
        xlims = xlim, ylims = ylim_row3)
    if !isempty(bankrupt_xs)
        Plots.vspan!(p3, bankrupt_xs, color = :grey, alpha = 0.25,
            linewidth = 0, label = is_right ? "Bankrupt" : "")
    end
    Plots.scatter!(p3, xt[t0:T], savings_p, color = :black,
        label = is_right ? "Savings" : "")
    Plots.plot!(p3, xt[t0:T], savings_p, label = "", color = :black)
    Plots.scatter!(p3, xt[t0:T], ηt_p, color = :grey,
        label = is_right ? "Premiums" : "", alpha = 0.5)
    Plots.plot!(p3, xt[t0:T], ηt_p, label = "", color = :grey)

    push!(row1, p1); push!(row2, p2); push!(row3, p3)
end

fig = Plots.plot(row1..., row2..., row3...;
    layout = (3, 3),
    markersize = 3.0, linewidth = 0.75,
    size = (1400, 750),
    legend_foreground_color = :transparent,
    legendfontsize = 9, tickfontsize = 9,
    guidefontsize = 12, titlefontsize = 13,
    left_margin = 8mm, bottom_margin = 6mm,
    top_margin = 3mm, right_margin = 3mm)

Plots.savefig(fig, "manuscript/figure_S10.png")
fig
