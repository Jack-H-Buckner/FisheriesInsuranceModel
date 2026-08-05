#####################################################################
#####################################################################
### Policy summaries read off a solved index insurance problem.
###
### The full model is defined in index_model/model.jl.
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################


"""
    coverage_levels(prob)

Optimal premium spending and the actuarially fair price per unit exposure over
the non bankrupt state grid of a solved problem.

Returns a NamedTuple with

  η        premium spending, size (Ns, Nb, NΔb, 4)
  η̄        price per unit exposure, same size
  coverage η ./ η̄, the units of exposure bought
  grid     the (s, b, Δb, M) axes, taken from the value function itself

The axes come from `prob.V`, not from a grid passed in by the caller, so they
cannot disagree with the grid the problem was solved on.  For this model `El`
is state independent and equals the payout rate p_index, so η̄ is flat and all
of the variation in `coverage` is demand, not price.
"""
function coverage_levels(prob)
    # non bankrupt states only; the bankrupt branch carries a degenerate grid
    states = prob.V.states[:, prob.V.states[1, :] .== 1]

    # RegularGridBspline stores the axis ranges in `grid` and their lengths in
    # `dims` (see ValueFunctions.jl); `dims` is what reshape needs.
    bspline = prob.V.Bslines[1]
    dims = bspline.dims
    grid = bspline.grid

    u_opt = ValueFunctionIterations.policy(states, prob.u, prob.X, prob.R,
                                           prob.F, prob.p, prob.δ, prob.V)
    η_s = reshape(u_opt[2, :], dims...)
    η̄_s = mapslices(s -> price_per_exposure(El(s, prob.p), prob.p.cf, prob.p.cv),
                    states, dims = 1)
    η̄_s = reshape(η̄_s, dims...)

    return (η = η_s, η̄ = η̄_s, coverage = η_s ./ η̄_s, grid = grid)
end
