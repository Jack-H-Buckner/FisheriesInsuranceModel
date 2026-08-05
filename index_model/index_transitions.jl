#####################################################################
#####################################################################
### Builds the 4 state joint transition matrix for the index
### insurance model from interpretable contract quality parameters
### and swaps it into a base parameter set.
###
### The transition algebra itself lives in derived_parameters.jl;
### this file only handles validation, ComponentArray surgery and
### diagnostics.
### Jack H. Buckner, Oregon State University, 2025
### generated with Claude Code
#####################################################################
#####################################################################

using ComponentArrays, LinearAlgebra

include("derived_parameters.jl")

# The joint state is indexed j = 2*(M-1) + x with
#
# state | 1 | 2 | 3 | 4 |
# ------|---|---|---|---|
# M     | 1 | 1 | 2 | 2 |
# x     | 1 | 2 | 1 | 2 |
#
# M = 1 typical mortality, M = 2 high mortality
# x = 1 no payout,         x = 2 payout
#
# The index informs transitions *into* the high mortality state from the
# low mortality state.  Conditional on already being in the high mortality
# state it carries no information, so the precision and recall below are
# defined conditional on M = 1.


# Derived quantities implied by a contract design.  Shared by the
# validation, the diagnostics and the scenario generator so the algebra
# is written down exactly once.
#
# p_index  P(payout)          = the actuarially fair premium rate
# b0       P(event | no payout) = false negative leakage
# fpr      P(payout | no event) = false positive rate
function index_moments(pr, re, p_low, p_stay)
    p_index = re * p_low / pr
    b0      = p_low * (1 - re) / (1 - p_index)
    fpr     = (p_index - re * p_low) / (1 - p_low)
    return (
        p_index       = p_index,
        b0            = b0,
        fpr           = fpr,
        return_period = 1 / p_index,
        premium_ratio = p_index / p_low, # premium relative to a perfect index
        p_low         = p_low,
        p_stay        = p_stay,
        pr            = pr,
        re            = re,
    )
end


# Check a contract design before it is used to build a transition matrix.
# p_index >= 1 is fatal: it is the denominator of b0.  p_index above
# max_p_index is legal but is not really insurance, so warn.
function validate_index_parameters(pr, re, p_low, p_stay; max_p_index = 0.25)
    (0 < pr <= 1)     || throw(ArgumentError("precision pr must be in (0,1], got $pr"))
    (0 < re <= 1)     || throw(ArgumentError("recall re must be in (0,1], got $re"))
    (0 < p_low < 1)   || throw(ArgumentError("p_low must be in (0,1), got $p_low"))
    (0 <= p_stay < 1) || throw(ArgumentError("p_stay must be in [0,1), got $p_stay"))

    p_index = re * p_low / pr
    p_index < 1 || throw(ArgumentError(
        "p_index = re*p_low/pr = $p_index must be < 1; " *
        "raise precision or lower recall"))
    if p_index > max_p_index
        @warn "payout rate p_index = $(round(p_index, digits=4)) " *
              "(1 in $(round(1/p_index, digits=1)) years) exceeds max_p_index = $max_p_index" pr re p_low
    end
    return nothing
end


# Check the assembled matrix.  Column stochasticity is the important one:
# MarkovChain in src/random_variables.jl throws above a 1e-6 error, and a
# matrix that is subtly off produces silently wrong expectations.
function validate_joint_transition(T, p_low, p_stay; tol = 1e-12)
    size(T) == (4, 4) || throw(ArgumentError("T must be 4x4, got $(size(T))"))
    all(0 .<= T .<= 1) || throw(ArgumentError("T has entries outside [0,1]"))

    colsums = vec(sum(T, dims = 1))
    err = maximum(abs.(colsums .- 1))
    err <= tol || throw(ArgumentError(
        "columns of T must sum to 1, max error $err > $tol; colsums = $colsums"))

    # The mortality marginal must be free of the contract design so that
    # event frequency and duration are held fixed across the sweep.
    T_M      = [1-p_low 1-p_stay; p_low p_stay]
    marginal = vcat((T[1, :] .+ T[2, :])', (T[3, :] .+ T[4, :])')
    expected = T_M[:, [1, 1, 2, 2]]
    merr = maximum(abs.(marginal .- expected))
    merr <= tol || throw(ArgumentError(
        "mortality marginal of T does not match [1-p_low 1-p_stay; p_low p_stay], " *
        "max error $merr > $tol"))

    return nothing
end


# Convert a ComponentArray back to a NamedTuple so that a field can be
# replaced with an array of a different size.
recursive_nt(x) = x
recursive_nt(x::ComponentArray) =
    NamedTuple{propertynames(x)}(map(k -> copy(recursive_nt(getproperty(x, k))),
                                     ComponentArrays.valkeys(x)))


"""
    set_index_T(p; pr, re, p_low = p.T[2,1], p_stay = p.T[2,2])

Return a copy of the parameter set `p` with the inherited 2x2 mortality
chain replaced by the 4x4 joint chain over the index insurance state
`2*(M-1) + x`.

  pr      precision, P(high mortality next | payout), given M = 1
  re      recall,    P(payout | high mortality next), given M = 1
  p_low   onset rate,   P(M'=2 | M=1)
  p_stay  persistence,  P(M'=2 | M=2)

`p_low` and `p_stay` default to the onset and persistence rates already in
`p.T`, so the base parameter file remains the single source of truth for
event frequency and duration.
"""
function set_index_T(p; pr, re, p_low = p.T[2, 1], p_stay = p.T[2, 2], max_p_index = 0.25)
    validate_index_parameters(pr, re, p_low, p_stay; max_p_index = max_p_index)
    T = joint_transition(p_low, pr, re, p_stay)
    validate_joint_transition(T, p_low, p_stay)
    return ComponentArray(merge(recursive_nt(p), (T = T,)))
end


# Stationary distribution of a column stochastic transition matrix.
function stationary_distribution(T)
    ns = nullspace(T - I)
    size(ns, 2) == 1 || throw(ArgumentError(
        "T has a $(size(ns,2))-dimensional stationary space; chain is reducible"))
    π = vec(real.(ns[:, 1]))
    return π ./ sum(π)
end


safediv(a, b) = b == 0 ? NaN : a / b

# Realised contract properties read back out of an assembled T, reported per
# current joint state.  Precision and recall match (pr, re) only in states 1
# and 2 (M = 1); in states 3 and 4 the index is uninformative by design, so
# precision collapses to the base rate p_stay.
function index_diagnostics(T)
    j = 1:4
    payout    = [T[2, i] + T[4, i] for i in j] # P(x'=2) = El
    event     = [T[3, i] + T[4, i] for i in j] # P(M'=2)
    precision = [safediv(T[4, i], T[2, i] + T[4, i]) for i in j]
    recall    = [safediv(T[4, i], T[3, i] + T[4, i]) for i in j]
    b0        = [safediv(T[3, i], T[1, i] + T[3, i]) for i in j]
    fpr       = [safediv(T[2, i], T[1, i] + T[2, i]) for i in j]
    return (
        state      = collect(j),
        M          = [1, 1, 2, 2],
        x          = [1, 2, 1, 2],
        payout     = payout,
        event      = event,
        precision  = precision,
        recall     = recall,
        b0         = b0,
        fpr        = fpr,
        stationary = stationary_distribution(T),
    )
end


# Human readable dump of index_diagnostics, printed by the run scripts so
# every solve is self documenting in its job log.
function print_index_diagnostics(T; io = stdout)
    d = index_diagnostics(T)
    println(io, "index diagnostics (columns = current joint state)")
    println(io, rpad("  state", 14), join(lpad.(string.(d.state), 10)))
    println(io, rpad("  M", 14),     join(lpad.(string.(d.M), 10)))
    println(io, rpad("  x", 14),     join(lpad.(string.(d.x), 10)))
    for (name, v) in (("El = P(payout)", d.payout), ("P(event)", d.event),
                      ("precision", d.precision), ("recall", d.recall),
                      ("b0", d.b0), ("fpr", d.fpr), ("stationary", d.stationary))
        println(io, rpad("  " * name, 14), join(lpad.(string.(round.(v, digits = 5)), 10)))
    end
    return nothing
end
