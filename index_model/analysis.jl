function coverage_levels(prob)
    # states and dimensions 
    states = prob.V.states[:,prob.V.states[1,:].==1]
    dims = prob.V.Bslines[1].dims

    # policy 
    u_opt = ValueFunctionIterations.policy(states,prob.u,prob.X,prob.R,prob.F,prob.p,prob.δ,prob.V)
    η_s = reshape(u_opt[2,:],dims...)
    η̄_s =  mapslices(s -> price_per_exposure(El(s,prob.p),prob.p.cf,prob.p.cv),states,dims = 1)
    η̄_s = reshape(η̄_s,dims...)

    return η_s, η̄_s, state_grid
end 