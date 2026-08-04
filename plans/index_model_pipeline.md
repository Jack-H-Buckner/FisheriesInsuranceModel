# Index-insurance model: HPC solve pipeline + demand analysis

## Context

`index_model/` is a variant of the base model in `src/`. Instead of an index built from
observed biomass plus basis-risk noise, the contract pays a fixed \$1 per unit exposure
whenever a 4-state Markov chain lands in a payout state. Basis risk arises because the
payout state is only imperfectly correlated with the high-natural-mortality state.

Right now the index model exists only as a notebook (`index_model/test_intermidiate.ipynb`):
`T` is hand-typed in a cell, the solve runs interactively (~24 min at coarse resolution),
and nothing is saved. There is no way to sweep contract designs, no HPC entry point, and
no welfare / coverage / bankruptcy analysis for this model.

Goal: a reproducible pipeline that (a) defines index-contract scenarios as files, (b) solves
them on the CQLS HPC via the existing `Makefile`/`hqsub` pattern, (c) syncs solutions locally,
and (d) computes the three demand metrics — average coverage, change in welfare, probability
of bankruptcy — for the `no_RA_*` parameter sets, all averaged over the stationary
distribution of biological states at a fixed wealth level, exactly as the base model does.

---

## Decisions already made

| Question | Choice |
|---|---|
| Scenario `T` parameterization | `joint_transition(p_low, pr, re, p_stay)` — precision / recall |
| Scenario list | Two arms over `pr ∈ {1.00, 0.75, 0.60, 0.50}`: `re = pr` and `re = 1` |
| Scheduler | CQLS `hqsub` — extend the existing `Makefile` |
| Layout | Self-contained `index_model/` tree; solutions in `results_index/` |
| Solve coverage | `{RA, no_RA} × {insurance, no-insurance}` per scenario; metrics on `no_RA` |
| Grid resolution | Match base model: `Nκ=35 Nη=25 Ns=24 Nb=24 NΔb=18 Nquad=50` |

---

## The transition matrix: `joint_transition(p_low, pr, re, p_stay)`

`index_model/derived_parameters.jl` now builds `T` from four interpretable scalars. The
index informs only the *transition into* the high-mortality state; conditional on already
being in the high-mortality state it carries no information.

| symbol | meaning |
|---|---|
| `p_low` | `P(M'=2 \| M=1)` — onset rate of a high-mortality event |
| `p_stay` | `P(M'=2 \| M=2)` — persistence |
| `re` | **recall** = `P(payout \| onset)` |
| `pr` | **precision** = `P(onset \| payout)` |

with the two derived quantities

```
p_index = re·p_low / pr          # unconditional payout rate  = the fair premium rate
b0      = p_low(1-re)/(1-p_index) # false-negative leakage: P(onset | no payout)
```

Four properties of this parameterization that drive the plan:

1. **The marginal mortality chain is exactly `T_M = [1-p_low  1-p_stay; p_low  p_stay]`,
   independent of `(pr, re)`.** Verified: `(1-p_index)·b0 + p_index·pr = p_low(1-re) + re·p_low
   = p_low`. So event frequency and duration are held fixed while basis risk varies — the
   sweep is clean.
2. **`El` is state-independent and equals `p_index`.** `El(s,p) = sum(p.T[[2,4], M])` collapses
   to `T_x[2,x]·colsum(T_mx1)[M] = p_index` for all four joint states. The fair premium
   `η̄ = (p_index + cf)/(1-cv)` is therefore flat, so coverage heatmaps over `(s, b)` show
   pure demand variation, not price variation. (Contrast the base model, where `El` varies
   with `b`, `Δb`, `M`.)
3. **The perfect index is reachable.** `pr = re = 1` ⇒ `p_index = p_low`, `b0 = 0` — payout
   if and only if onset. This was *not* reachable under the previous `(tpr, fpr)` form.
4. **The uninformative benchmark is `pr = p_low`.** Then `p_index = re`, `b0 = p_low`, and
   `T_mx0`'s first column equals `T_mx1`'s — the payout carries zero information about `M'`.
   Demand should go to zero here at `cv = 0`. This is the key falsification test.

**Anchor check (already verified analytically):** `joint_transition(0.05, 0.8, 0.8, 0.4)`
reproduces the hand-built blocks in `test_intermidiate.ipynb` cell 4 exactly —
`T_mx0 = [0.94/0.95 0.6; 1-0.94/0.95 0.4]`, `T_mx1 = [0.2 0.6; 0.8 0.4]`,
`T_x = [0.95 0.95; 0.05 0.05]`. And `joint_transition(0.05, 1.0, 0.5, 0.5)` reproduces the
matrix printed in `test_transition_matrices.ipynb` cell 2. Both become regression tests.

**Feasibility.** `set_index_T` must assert `0 < pr ≤ 1`, `0 < re ≤ 1`, and
`p_index = re·p_low/pr < 1` (strict — it is the denominator of `b0`). In practice constrain
`p_index ≤ 0.25`: a contract paying out more than one year in four is not insurance. This
binds only in the low-precision / high-recall corner.

**Note:** the 3-argument `joint_transition(T_mx0, T_mx1, T_x)` was renamed
`joint_transition_from_matrix`, so `test_intermidiate.ipynb` cell 4 no longer runs as
written. Not worth fixing — the notebook is superseded by the scenario files below.

---

## Proposed precision × recall grid

`p_low` and `p_stay` are **read from the base parameter file's inherited 2×2 `T`**
(`p.T[2,1]`, `p.T[2,2]`), so an index scenario layered on `no_RA_intermidiate.jl` gets
`(0.025, 0.50)` and one layered on `..._higher_frequency.jl` automatically gets `(0.05, 0.50)`.
No duplication of frequency/duration constants.

Rather than a factorial, the sweep is **two arms sharing one precision ladder**

```
pr ∈ {1.00, 0.75, 0.60, 0.50}
```

**Arm A — iso-premium diagonal (`re = pr`).** Payout frequency is pinned to event frequency,
so `p_index = p_low = 0.025` (40-year return period) in every case. Cost is held fixed and
only accuracy moves. This isolates the pure basis-risk effect.

**Arm B — full coverage (`re = 1.00`).** Every event is covered; declining precision buys
that coverage with false positives, so the premium rises. This isolates the pure price effect
of basis risk.

The two arms meet at `(pr, re) = (1.00, 1.00)` — the perfect index — so there are **7 unique
scenarios**, not 8.

| arm | `pr` | `re` | `p_index` | return period | premium vs. perfect | `fpr` | `b0` |
|---|---|---|---|---|---|---|---|
| A & B | 1.00 | 1.00 | 0.0250 | 40 yr | 1.00× | 0 | 0 |
| A | 0.75 | 0.75 | 0.0250 | 40 yr | 1.00× | 0.0064 | 0.0064 |
| A | 0.60 | 0.60 | 0.0250 | 40 yr | 1.00× | 0.0103 | 0.0103 |
| A | 0.50 | 0.50 | 0.0250 | 40 yr | 1.00× | 0.0128 | 0.0128 |
| B | 0.75 | 1.00 | 0.0333 | 30 yr | 1.33× | 0.0085 | 0 |
| B | 0.60 | 1.00 | 0.0417 | 24 yr | 1.67× | 0.0171 | 0 |
| B | 0.50 | 1.00 | 0.0500 | 20 yr | 2.00× | 0.0256 | 0 |

(`fpr` = `P(payout | no event)`, `b0` = `P(event | no payout)`.)

Why this precision ladder rather than even spacing on `[0.5, 1.0]`: on arm B the premium
multiplier is `1/pr`, so `{1.00, 0.75, 0.60, 0.50}` puts the *premium* on an evenly spaced
`{1.00, 1.33, 1.67, 2.00}×` ladder. Since arm B is defined by the premium response, spacing
the economically active variable evenly is the better design. Using the same `pr` values on
arm A then lets you read the two arms against each other at matched precision — "hold cost
fixed and lose recall" vs. "hold recall and pay more" — which is the comparison that
actually answers whether basis risk bites through hedge quality or through price.

Two properties worth knowing when reading results:

- **On arm A, `fpr == b0` exactly.** With `p_index = p_low`, both reduce to
  `p_low(1-pr)/(1-p_low)`. The diagonal is a symmetric-error ladder: each step adds equal
  false-positive and false-negative rates. Clean to interpret, and a good numerical check.
- **On arm B, `b0 = 0` exactly** for every case — there are no uncovered events at all, so
  any decline in demand along arm B is attributable to price alone.

One scenario outside both arms anchors the bottom of the quality range:

- **`pr = p_low = 0.025, re = 0.50`** ⇒ `T_mx0[:,1] == T_mx1[:,1]`, an uninformative index.
  Demand must be ≈0 at `cv = 0`. Falsification test, not a design point.

Your prior that `pr < 0.5` gives zero optimal coverage is plausible but untested — arm B at
`pr = 0.50` already doubles the premium for the same protection. If coverage is still
substantial there, extending the ladder downward costs 2 jobs per added point, so leave that
decision until the first results land.

**Job count.** 7 scenarios × 2 base parameter sets = 14 insurance solves, plus the
uninformative pair = 16, plus 2 no-insurance solves = **18 jobs**. Submit in one go.

**Compute saving — solve the no-insurance problem once, not once per cell.** With
`η_grid = [0.0]`, `index_insurance(·, 0, η̄) = 0`, so the index never enters `F`; and because
`T_x` has identical columns, the `x`-dimension is non-persistent and `M' | M` is the marginal
`T_M` regardless of `(pr, re)`. The no-insurance value function therefore depends only on
`(p_low, p_stay)` and is **identical across both arms** — 2 solves, not 16. Step 6 asserts
this empirically rather than assuming it (solve two scenarios' no-insurance problems, check
`V` matches to `1e-8`).

---

## New / changed files

```
index_model/
  index_transitions.jl        NEW  set_index_T wrapper + validation + index_diagnostics
  grids.jl                    NEW  build_grids(params) — SINGLE source of truth
  parameters/
    index_universal_params.jl NEW  contract cost structure (cf, cv)
    generate_scenarios.jl     NEW  writes the precision x recall grid
    <scenario>.jl             GEN  one file per (base, pr, re) cell
  run.jl                      NEW  solve with insurance      -> results_index/<s>/
  run_no_insurance.jl         NEW  solve with η_grid = [0.0] -> results_index/no_insurance_<s>/
  analysis/
    setup.jl                  NEW  rebuild_problem(scenario) + load_solution
    coverage_levels.jl        NEW  metric 1
    value_diffs.jl            NEW  metric 2
    bankruptcy.jl             NEW  metric 3
    run_analysis.jl           NEW  driver -> manuscript_index/*.csv
  test/test_transitions.jl    NEW  T validity + tpr/fpr recovery + notebook anchor
  model.jl                    FIX  environment_mortaltiy -> _index (2 sites); error text
  plotting.jl                 FIX  environment_mortaltiy -> _index; bare global p.s̄; wt length
  analysis.jl                 FIX  coverage_levels returns undefined global state_grid

Makefile                      EDIT add setup_index / run_index / run_index_params / sync_index
.gitignore                    EDIT ignore results_index/ (results/ is already 309 MB untracked)
```

Nothing under `src/`, `parameters/`, `results/`, `analysis/`, or `manuscript/` is touched.

---

## Step 1 — `index_model/index_transitions.jl`

Thin wrapper over the existing `joint_transition` — it does the `ComponentArray` surgery and
the validation, nothing more. The transition algebra stays in `derived_parameters.jl`.

```julia
include("derived_parameters.jl")   # joint_transition, joint_transition_from_matrix

"""
    set_index_T(p; pr, re, p_low = p.T[2,1], p_stay = p.T[2,2])

Replace the inherited 2x2 mortality chain in `p` with the 4x4 joint chain over
state 2*(M-1)+x, M in {1,2} mortality (1 typical, 2 high), x in {1,2} payout.

`p_low` and `p_stay` default to the onset and persistence rates already in `p.T`,
so the base parameter file remains the single source for event frequency and duration.
"""
function set_index_T(p; pr, re, p_low = p.T[2,1], p_stay = p.T[2,2]) ... end
```

Validation, all cheap and all worth having (a bad `T` fails deep inside `MarkovChain` with
an opaque message otherwise):

- `0 < pr ≤ 1`, `0 < re ≤ 1`, `0 < p_low < 1`, `0 ≤ p_stay < 1`
- `p_index = re*p_low/pr < 1` strictly (it is `b0`'s denominator), and warn if `> 0.25`
- `all(sum(T, dims=1) .≈ 1)` to `1e-12` — `MarkovChain` (`src/random_variables.jl:19`)
  throws above `1e-6`
- `0 .≤ T .≤ 1`
- the M-marginal `[T[1,:]+T[2,:]; T[3,:]+T[4,:]]` collapses to `[1-p_low 1-p_stay; p_low p_stay]`

Rebuilding the `ComponentArray` with a 4×4 `T` in place of the inherited 2×2 needs the
`recursive_nt` trick from `test_intermidiate.ipynb` cell 3 — lift it here so scenario files
stay one-liners:

```julia
recursive_nt(x) = x
recursive_nt(x::ComponentArray) =
    NamedTuple{propertynames(x)}(map(k -> copy(recursive_nt(getproperty(x, k))), valkeys(x)))
```

Also `index_diagnostics(T)` → `p_index`, `b0`, realised precision/recall per current state,
`El` per state, implied `fpr = (p_index - re·p_low)/(1-p_low)`, and the stationary
distribution of the 4-state chain. `run.jl` prints it into the job log so every solve is
self-documenting.

## Step 2 — scenario files

`index_model/parameters/index_universal_params.jl` — contract cost structure only
(`p_low`/`p_stay` now come from the base parameter file, so nothing else belongs here):

```julia
index_cf = 1e-6   # base-model default; the notebook used cf = 0.0
index_cv = 0.0    # base-model default; the notebook used cv = 0.1
```

Each scenario file, e.g. `index_model/parameters/no_RA_int_pr50_re75.jl`:

```julia
# parameters/*.jl call mean_price -> baranov, so the population model must be
# loaded first. src/run.jl relies on this ordering implicitly; make it explicit
# here so a scenario file can be included standalone (e.g. by a test).
include(joinpath(@__DIR__, "..", "population_model.jl"))
include(joinpath(@__DIR__, "..", "..", "parameters", "no_RA_intermidiate.jl"))  # defines p
include(joinpath(@__DIR__, "index_universal_params.jl"))
include(joinpath(@__DIR__, "..", "index_transitions.jl"))

p = set_index_T(p; pr = 0.50, re = 0.75)   # p_low, p_stay inherited from p.T
p.cf = index_cf; p.cv = index_cv
```

These are mechanical, so generate them rather than hand-writing 16 files —
`index_model/parameters/generate_scenarios.jl` holds the precision ladder, the two arm
definitions and the base list in one place and writes one file per `(base, pr, re)` scenario,
deduplicating the `(1.00, 1.00)` cell the arms share. Rerunning it after you edit the ladder
is the supported way to change the sweep.

Naming: `{RA|no_RA}_int_pr<100·pr>_re<100·re>.jl`, e.g. `no_RA_int_pr50_re100.jl`; the
uninformative case is `no_RA_int_uninformative.jl`, and the two regression anchors are
`no_RA_int_anchor_notebook.jl` (`p_low=0.05, pr=0.8, re=0.8, p_stay=0.4`) and
`no_RA_int_anchor_transitions.jl` (`p_low=0.05, pr=1.0, re=0.5, p_stay=0.5`). The
`no_insurance_` prefix is added by the runner, never by a filename, so `ls parameters/`
maps 1:1 to scenarios.

**`occursin` hazard:** `src/run.jl:26-31` picks grid bounds by substring-matching `"short"` /
`"long"` in the scenario name. Scenario tags must avoid those substrings; `grids.jl` will
`@assert` this rather than silently using the wrong bounds.

## Step 3 — `index_model/grids.jl`

The base model's single biggest reproducibility hazard: only `V` and `P` are serialized
(`save_solution` in `ValueFunctionIterations/src/DynamicPrograms.jl:277`), so every consumer
must rebuild the problem with byte-identical grids — and today `src/run.jl` (20/20/15),
`src/analysis.jl` (20/2/24/24/18/50) and `analysis/table_4.jl` (35/25/24/24/18/50) all
disagree, silently.

Fix for the index model: one function used by the runners *and* every analysis script.

```julia
function build_grids(params::String; Nκ, Nη, Ns, Nb, NΔb, insurance::Bool)
    # bmax/bmin = 0.5/-3.5, Δbmax/Δbmin = 1.25/-1.0, smax/smin = 2.5/-4.0
    # ηmax = insurance ? 0.2 : 0.0
    return (action_grid = [κ_grid, η_grid], state_grid = [s_grid, b_grid, Δb_grid], δ = 0.95)
end
```

`run.jl` additionally writes `results_index/<scenario>/grid.toml` recording all six sizes,
the bounds, `δ`, `Nquad`, and the resolved 4×4 `T`. Analysis scripts read it back and
`@assert` the rebuilt grid matches. This is what makes the pipeline reproducible.

## Step 4 — runners

`index_model/run.jl`, modelled on `src/run.jl` with four deliberate changes:

```julia
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))      # absolute, not Pkg.activate(".")
const ROOT = normpath(joinpath(@__DIR__, ".."))
params, Nκ, Nη, Ns, Nb, NΔb, Nquad = ARGS
include(joinpath(@__DIR__, "model.jl"))                 # @__DIR__, not "../src/..."
include(joinpath(@__DIR__, "plotting.jl"))
include(joinpath(@__DIR__, "grids.jl"))
include(joinpath(@__DIR__, "parameters", "$params.jl"))

outdir = joinpath(ROOT, "results_index", params)
mkpath(outdir)                                          # mkpath, not mkdir (race-safe)

g = build_grids(params; ..., insurance = true)
X, Xmc = init_ranom_variables(p; Nquad = parse(Int, Nquad))
prob = init(g.action_grid, g.state_grid, g.δ, p, X, index)   # index model: no strike arg
solve!(prob, maxiter = round(Int, 5/(1-g.δ)))
save_solution(prob, joinpath(outdir, "sol.jld2"))
write_grid_meta(outdir, g, p, Nquad)
```

1. No `Pkg.instantiate()` in the worker — 40+ concurrent jobs racing on the depot while
   fetching the unregistered `ValueFunctionIterations.jl` git dep is a real failure mode.
   Moved to a one-shot `make setup_index`.
2. All paths derived from `@__DIR__`, so cwd no longer matters.
3. `init(...)` takes 6 args (no `strike`) and returns `prob` alone, not `(prob, El)` —
   `index_model/model.jl:70` vs `src/model.jl:69`.
4. Writes `grid.toml` and prints `index_diagnostics(p.T)` to stdout.

`index_model/run_no_insurance.jl` is identical with `insurance = false`. Because the
no-insurance solution is invariant to `(pr, re)` (see the grid section), it is keyed on the
**base** parameter set and the mortality chain only, not the scenario:
`outdir = results_index/no_insurance_<base>_plow<...>_pstay<...>/`. `run_no_insurance.jl`
takes the base name (`no_RA_intermidiate`) rather than a scenario name and is submitted
once per base, not once per cell. `setup.jl` resolves a scenario to its matching
no-insurance directory by reading `p_low`/`p_stay` back out of the scenario's `grid.toml`.

Both runners take the same 7 positional args (the base `Makefile`'s `run:` target passes
only 1 to `run_no_insurance.jl`, which errors — the new targets pass all 7).

## Step 5 — bug fixes in existing `index_model/` files

These are real defects, not cosmetics — the dynamics in `update_stock` are correct, but
every diagnostic path uses the wrong mortality mapping:

- `index_model/model.jl:35` and `:52` — `simulate_trajectory` calls
  `environment_mortaltiy(states[3,t+1], p.m̄, p.m̲)`, the **2-state** version from
  `src/population_model.jl`. Under the 4-state labelling it returns `m̲` (high mortality)
  for `M=2`, which is a *typical*-mortality-with-payout state. Must be
  `environment_mortaltiy_index` (`index_model/population_model.jl:26`). This corrupts the
  `m`/`h`/`y` series used to choose the `b`/`Δb` grid.
- `index_model/plotting.jl:59` — same call inside `simulation`, corrupting the reported
  `mt` series that `bankruptcy.jl` and the simulation figures consume.
- `index_model/analysis.jl:12` — `coverage_levels(prob)` returns `state_grid`, which is
  never defined in scope; it only worked because the notebook had a matching global.
  Change to `coverage_levels(prob, state_grid)`.
- `index_model/model.jl:117` — error text still says "Mt_1 must be one or two" though the
  guard now admits 1–4.
- `index_model/plotting.jl:107-108` — `plot_simulation` reads a bare global `p.s̄`.
- `index_model/plotting.jl` — `wt` has length `n-2` while `ηt`/`Zt` have `n-1`, so the
  default `T = length(x.wt)` silently truncates plots by one period.

## Step 6 — analysis (`no_RA_*` scenarios)

All three metrics share one stationary-state sample, matching `src/analysis.jl:288`:

```julia
Random.seed!(123)
sim_states, m, f, h, y = simulate_trajectory(p, Xmc, 2000)   # 3 x 2001 = (log b, Δlog b, M∈1:4)
```

**`analysis/setup.jl`** — `rebuild_problem(scenario; insurance)` reads `grid.toml`, rebuilds
via `build_grids`, asserts a match, calls `init`, and `load_solution`s the right `.jld2`.
Every metric script goes through it, so a grid mismatch is impossible rather than silent.

**Metric 1 — average coverage levels** (`coverage_levels.jl`). Reuse
`index_model/analysis.jl :: coverage_levels` (fixed per Step 5) for the full-grid heatmaps,
and add the stationary-weighted scalar the demand table needs:

```julia
u_opt = ValueFunctionIterations.policy(S, prob.u, prob.X, prob.R, prob.F, prob.p, prob.δ, prob.V)
η̄    = mapslices(s -> price_per_exposure(El(s, prob.p), prob.p.cf, prob.p.cv), S, dims = 1)
mean_coverage = mean(u_opt[2,:] ./ vec(η̄))
```
where `S = [vcat([1.0, log(s0 - p.s̄)], sim_states[:,t]) for t in 1:2000]` at each fixed
wealth level. Coverage is `η / η̄` throughout (premium spend ÷ price per unit exposure).

**Metric 2 — change in welfare** (`value_diffs.jl`). Port `src/analysis.jl:205`,
substituting `rebuild_problem` for its hardcoded grids:

```julia
certianty_equivelent(prob, state) = inverse_utility(prob.V(state)*(1-prob.δ), prob.p.γ)

function average_diff(s, states)
    st   = log(s - p.s̄)
    V    = mean(mapslices(x -> certianty_equivelent(prob_ins,    vcat([1.0,st], x)), states, dims=1))
    V_no = mean(mapslices(x -> certianty_equivelent(prob_no_ins, vcat([1.0,st], x)), states, dims=1))
    return (V - V_no) / V_no
end
st_vals = [0.25, 0.5, 1.0, 2.0, 3.5] .+ p.s̄
```

Percent increase in certainty-equivalent consumption, averaged over the stationary
biological distribution, at five fixed wealth levels (`s = 1.0 + s̄` is the headline row).
`certianty_equivelent` touches only `prob.V`, so the action grids differ correctly between
the two problems without affecting the welfare number.

**Metric 3 — probability of bankruptcy** (`bankruptcy.jl`). Port `analysis/table_4.jl:8`
(the production version — `src/calculate_bankruptcy_rate.jl` is a broken stub with
`N=2, T=2`, an unconditional `throw()`, and an R-ism `data.frame(df)`). Unchanged method:
`N=60` paired trajectories × `T=250` periods, per-`i` `Random.seed!(hash((123,i)) % UInt32)`
inside `Threads.@threads`, Poisson MLE hazard `λ̂ = events / Σ exposure`, `B=2000` paired
nonparametric bootstrap, reported as `P10 = 1 - exp(-10λ̂)` with a percentile CI. Only the
`init`/`load_solution` calls change (no `strike`, `results_index/` paths).

**`analysis/run_analysis.jl`** loops the `no_RA_*` scenarios and writes one tidy CSV per
metric plus a combined demand table to `manuscript_index/`:

```
manuscript_index/index_coverage.csv    scenario, pr, re, p_index, savings, mean_coverage
manuscript_index/index_value_diff.csv  scenario, pr, re, p_index, savings, pct_increase_CE
manuscript_index/index_bankruptcy.csv  scenario, pr, re, p_index, P10, P10_lo, P10_hi, P10_no_I, ..., lambda_ratio
manuscript_index/index_demand.csv      one row per scenario at s = 1.0 + s̄, all three metrics
```

Every table carries `arm` (`"iso_premium"` / `"full_coverage"`), `pr`, `re` and the derived
`p_index`, so both arms plot against the shared precision ladder straight from the CSV. The
`(1.00, 1.00)` scenario is solved once but emitted under both arm labels so each curve has
its perfect-index endpoint.

## Step 7 — Makefile + sync

Add alongside the existing targets (leaving them untouched):

```make
CLUSTER_HOST=bucknejo@<cqls-host>
INDEX_PARAMS=$(DIRPATH)FisheriesInsuranceModel/index_model/parameters
IDX_NKAPPA=35
IDX_NETA=25
IDX_NS=24
IDX_NB=24
IDX_NDELTAB=18
IDX_NQUAD=50

setup_index:
	julia --project=$(DIRPATH)FisheriesInsuranceModel -e 'using Pkg; Pkg.instantiate()'

IDX_BASES=no_RA_intermidiate RA_intermidiate

# one no-insurance solve per BASE parameter set (invariant to pr, re)
run_index_no_insurance:
	@for base in $(IDX_BASES); do \
	  hqsub -P $(PROCS) "julia --project=$(DIRPATH)FisheriesInsuranceModel --threads=$(PROCS) $(DIRPATH)FisheriesInsuranceModel/index_model/run_no_insurance.jl $${base} $(IDX_NKAPPA) $(IDX_NETA) $(IDX_NS) $(IDX_NB) $(IDX_NDELTAB) $(IDX_NQUAD)" -r solve-index-no-insurance-$${base}-StdOut -q ceoas@$(NODE); \
	done;

# one insurance solve per scenario cell
run_index:
	@for params in `ls -1 $(INDEX_PARAMS) | grep '\.jl$$' | grep -v -e universal -e generate | sed 's/\.jl//g'`; do \
	  hqsub -P $(PROCS) "julia --project=$(DIRPATH)FisheriesInsuranceModel --threads=$(PROCS) $(DIRPATH)FisheriesInsuranceModel/index_model/run.jl $${params} $(IDX_NKAPPA) $(IDX_NETA) $(IDX_NS) $(IDX_NB) $(IDX_NDELTAB) $(IDX_NQUAD)" -r solve-index-$${params}-StdOut -q ceoas@$(NODE); \
	done;

run_index_params:
	hqsub ... index_model/run.jl ${PARAMS} ...
	hqsub ... index_model/run_no_insurance.jl ${PARAMS} ...

sync_index:
	rsync -avz --progress $(CLUSTER_HOST):$(DIRPATH)FisheriesInsuranceModel/results_index/ ./results_index/
```

Deliberately fixed relative to the existing `run:` target: the missing `/` in the params
path, `--project=` so the env resolves off-cwd, `grep -v universal` so the shared params
file is not submitted as a scenario (the base model has a stray empty `results/universal_params/`
from exactly this), all 7 args passed to the no-insurance runner, and **tabs not spaces**
in `run_index_params` (the existing `run_params` target is space-indented and won't parse
under GNU make).

`sync_index` is one-directional cluster → local. Add `results_index/` to `.gitignore` —
`results/` is already 309 MB and untracked.

**Compute estimate.** 24×24×18×4 = 41,472 continuous states, ~2× the base model. At the
notebook's 5,392 states / 24 min on 8 threads, one solve is roughly 3–4 h on 40 threads.
The whole sweep is 18 jobs (see the scenario section). Run
`make run_index_params PARAMS=no_RA_int_anchor_notebook` alone first and check the wall
clock before submitting the rest.

---

## Verification

1. **Transition matrix** — `julia --project=. index_model/test/test_transitions.jl`, over the
   whole `(pr, re)` grid:
   - every column of `T` sums to 1 to within `1e-12` (`MarkovChain` in
     `src/random_variables.jl:19` throws above `1e-6`), and `0 ≤ T ≤ 1`;
   - the M-marginal recovers `[1-p_low 1-p_stay; p_low p_stay]` exactly, i.e. frequency and
     duration really are invariant to `(pr, re)`;
   - simulating 10^6 steps recovers the target `pr` and `re` to Monte-Carlo error;
   - **anchor 1:** `joint_transition(0.05, 0.8, 0.8, 0.4)` reproduces
     `[0.94 0.94 0.57 0.57; 0.01 0.01 0.03 0.03; 0.01 0.01 0.38 0.38; 0.04 0.04 0.02 0.02]`
     (`test_intermidiate.ipynb` cell 4);
   - **anchor 2:** `joint_transition(0.05, 1.0, 0.5, 0.5)` reproduces
     `[0.95 0.95 0.4875 0.4875; 0 0 0.0125 0.0125; 0.025 0.025 0.4875 0.4875; 0.025 0.025 0.0125 0.0125]`
     (`test_transition_matrices.ipynb` cell 2);
   - **corner cases:** `pr = re = 1` gives `b0 = 0` and payout iff onset; `pr = p_low` gives
     `T_mx0[:,1] == T_mx1[:,1]` (uninformative);
   - `El(s,p) = sum(p.T[[2,4], M]) == p_index` for all four states and is strictly positive,
     so `price_per_exposure` never falls back to its `1e-12` nugget.
2. **Local smoke solve** — run `index_model/run.jl` at notebook resolution
   (`15 10 12 14 8 10`) for one scenario. Confirm it writes `sol.jld2` + `grid.toml`,
   and reproduce the notebook's cell-20 check:
   `prob.F([0,0,0,0,1.0], [0.1,0.1], X_node_5, p) ≈ [1.0, 1.197, 0.2085, 0.2085, 1.0]`.
3. **Grid-contract check** — deliberately edit `grid.toml`, confirm `rebuild_problem`
   raises instead of returning wrong numbers.
4. **No-insurance invariance** — the load-bearing assumption behind the compute saving.
   Solve the no-insurance problem at two different `(pr, re)` cells with the same
   `(p_low, p_stay)` and assert the value functions match to `1e-8`. If this fails, the
   reuse is wrong and the plan reverts to one no-insurance solve per cell.
5. **Economic sanity** — the `pr = p_low` scenario (uninformative index) should give ≈0 mean
   coverage and ≈0 welfare gain at `cv = 0`. `no_RA` is risk-neutral (`γ_no_RA = 0.0`), so
   demand there comes from the borrowing constraint and bankruptcy risk, not risk aversion.
   Both arms should decline monotonically in falling `pr`, and both must agree exactly at
   `(1.00, 1.00)` — that shared endpoint is a free consistency check across the two arms.
   Whether arm A or arm B falls faster is the open empirical question the sweep exists to
   answer, so do not treat either ordering as a bug.
6. **End-to-end** — `make setup_index && make run_index_no_insurance && make run_index` on
   the cluster (the no-insurance jobs are independent, so order does not matter), `make sync_index`
   locally, then `julia --project=. index_model/analysis/run_analysis.jl`, and confirm all
   four CSVs land in `manuscript_index/` with one row per `no_RA_*` scenario.

## Open items

- Both arms hold `(p_low, p_stay) = (0.025, 0.50)` from `universal_params.jl`. To also vary
  event frequency or duration, layer the same two arms on
  `no_RA_intermidiate_higher_frequency.jl` / `..._longer_duration.jl` — `set_index_T` picks
  the new `p_low`/`p_stay` up automatically, at 2–3× the job count.
- If arm B still shows substantial coverage at `pr = 0.50`, extend the ladder to 0.40 / 0.30
  before concluding where demand dies. 2 jobs per added point.
- `cv = 0.0` (actuarially fair) is the default. The notebook used `cv = 0.1`. A loaded
  premium would sharpen the precision/recall trade-off considerably; worth one extra column
  of scenarios once the fair-price surface is understood.
