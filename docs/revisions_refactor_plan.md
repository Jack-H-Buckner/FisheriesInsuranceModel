# Revisions refactor plan

## Context

Two model changes are planned. The user will hand-write much of the code, so this document is a map of where changes land, not an implementation.

1. **Fisher profit model.** Today effort `f` is a *deterministic* harvest control rule (HCR) of biomass. Prices and costs are constants, so profit does not drive behaviour. The new process in `notebooks/fishermen_profit_max.ipynb` has the fisher choose `f` to satisfy the profit FOC (`fisher_response`). The HCR then acts as a cap.
2. **Insurer capital costs.** Today the premium is a constant loading, `(EL+cf)/(1-cv)`, with `cv=0.33` for the "50% loading" scenarios. The new process in `notebooks/insurer_capital_model.ipynb` makes the loading an additive, state-dependent marginal capital cost, `price = EL + R'(Q*)`.

Line numbers refer to the repo as of 2026-11-16 (uncommitted changes included) and will drift.

## Current structure (what matters for the refactor)

- `index_model/` `include`s `src/population_model.jl`, `src/income_model.jl`, `src/derived_parameters.jl` and `src/random_variables.jl`. An edit to those `src/` files changes both the base model and the index model. Decide up front whether the base `src/` model gets the revision as well (see Open questions).
- Fisher profit and effort:
  - `src/population_model.jl:13-19` `control_rule(bt,b_target,f_target,b_limit)`.
  - `:23-25` `baranov`.
  - `:56-97` `update_stock`, which returns `(next_state, ht, ft)`.
  - `index_model/population_model.jl:57-70` overrides `update_stock`.
  - `src/income_model.jl:12` `income(ht,ft,pt,c_et) = pt*ht - c_et*ft`.
- Price and cost are constants:
  - `μ_P` comes from `mean_price`, `μ_Ce` from `cost_By0`, and `f_target` from `f_b_target` (`src/derived_parameters.jl:13-24`). These are set in `parameters/RA_intermidiate.jl:17-34`.
  - `src/random_variables.jl:16-39` fixes the price and effort-cost nodes at those means.
- Dynamic program: `src/model.jl:69-161` and `index_model/model.jl:74-188`.
  - The transition `F` calls `update_stock`, then `income`, then the claim, then `update_savings`.
  - The premium is computed at `src/model.jl:132` and `index_model/model.jl:153`.
- Premium: `price_per_exposure(EL,cf,cv)` is defined in `index_model/insurance_contracts.jl:36-40` and, as a legacy copy, in `src/insurnace_contracts.jl:29-33`. `El(s,p)` at `index_model/insurance_contracts.jl:47-50` depends only on the mortality state.

## Part A: Fisher profit-max (effort chosen by FOC, HCR as cap)

### A1. Core functions to add
- Port `harvest`, `dharvest_df`, `profit_foc`, `harvest_constaint` and `fisher_response` from the notebook into `src/population_model.jl`, with `using Roots`.
  - `baranov` already equals the notebook's `harvest`. Reuse it and drop the duplicate.
  - Fix the `harvest_constaint` typo while porting.
- Add `Roots` to `Project.toml`. It is not currently a dependency.
- `control_rule` stays. It becomes the cap `f_allowed(b)`, and `fisher_response` takes `min` over harvest. Resolve the cap-on-harvest versus cap-on-effort question in the notebook first (see Open questions).

### A2. Performance: do not root-find inside the Bellman loop
- `F` is evaluated for every state, action and quadrature node. A root find per call is too slow.
- Price and cost are constants. The optimal effort depends only on `b` (and on the mortality assumption used for the decision).
- **Recommended:** precompute `f*(b)` on a fine biomass grid at `init` or parameter-build time. Interpolate with `Interpolations` (already a dependency). The interpolant replaces `control_rule` inside `update_stock`.
- This keeps the DP structure unchanged. `update_stock` still returns `(next_state, ht, ft)`.

### A3. Places that must change
- `src/population_model.jl:56-97` and `index_model/population_model.jl:57-70`: swap `control_rule` for the interpolated `f*(b)`. Check the realized-mortality `m` against the decision `m`: the fisher chooses `f` before the mortality draw, but `ht` uses the draw.
- `src/income_model.jl:12` `income`: no change needed unless price or cost become stochastic. The Part A model needs no change here.
- Calibration, `src/derived_parameters.jl:13-24`, `parameters/universal_params.jl` and the per-scenario `parameters/*.jl`:
  - `f_b_target`, `cost_By0` and `mean_price` currently pin the parameters so that income at the target equals `y_target`.
  - With an FOC-driven `f`, `c` and `p` determine `f*(b)`. Decide which of `f_target`, `c` and `p` stays a primitive, and recalibrate (see Open questions).
  - About 30 files in `parameters/*.jl`, plus the generated `index_model/parameters/*.jl` from `index_model/parameters/generate_scenarios.jl`, inherit this.
- `src/random_variables.jl:16-39`: no change unless price or `Ce` become random.
- `src/model.jl` and `index_model/model.jl`: the `u`, `F` and `R` functions should not need structural changes. Retest the budget-constraint logic.
- `src/model_order_2.jl` and `src/population_model_order_2.jl` are variant copies. Decide whether they are still used and update them or retire them.

### A4. Downstream consumers that re-derive `ft`, `ht` or `yt`
- `src/plotting.jl:45-75` `simulation` rebuilds `ft`, `ht` and `yt` itself. `index_model/plotting.jl:37-71` has the same logic. Update both to use the same `f*(b)` call so they cannot drift from the model.
- Analysis code that calls `simulation`: `index_model/analysis/bankruptcy.jl` (lines 207 and 213), `src/calculate_bankruptcy_rate.jl`, `analysis/table_4.jl`, `analysis/table_S5.jl`, `analysis/figure_3.jl`, `figure_S2.jl`, `figure_S7.jl`, `figure_S9.jl`, `figure_S10.jl` and `figure_S11.jl`.
- `src/natural_capital.jl` has its own `period_payoff` (line 52), a payout-free DP that calls `update_stock`. It picks up the new effort rule automatically but needs rechecking.
- Scripts that hard-code biomass grids by scenario name (`short` or `long`) may need revisiting if `f*(b)` changes the stock dynamics range: `src/run*.jl` and `index_model/grids.jl` and `config.toml`.

### A5. New figure
- Promote the notebook plot (`h`, `f` and profit against `b`, one line per `c`) to a script or a supplementary figure.

## Part B: Insurer capital costs

### B1. What the capital model provides
`notebooks/insurer_capital_model.ipynb` defines these pieces.
- `lognormal_mean_cv`, `sample_L0` and `tvar`.
- `S` (SST capital requirement, `TVaR_α(L) - E[L]`), `R` (`τ·(S(Q)-S(0))`) and the derivative `dR`.
- `equilibrium_residual` and `equilibrium_Q(K,R,dR)`, which solve for the market equilibrium exposure `Q*` with `Roots.find_zero`.
- `capital_costs(s0, strike, pars, K, n, Q_base, CV, α, τ; Nmc, seed)`, which returns `R'(Q*)`, the marginal capital cost per unit of exposure.

The new premium is `price = EL + R'(Q*)`. It is additive and depends on state. The existing `(EL+cf)/(1-cv)` is multiplicative and constant.

### B2. Port into the repo
- Move the capital functions out of the notebook into a source file, for example `src/insurer_capital.jl` and `index_model/insurer_capital.jl`. Keep the notebook for exploration.
- Known issue to fix while porting: `capital_costs` overwrites its `α` and `τ` arguments with 0.99 and 0.06.

### B3. The index model does not match the notebook's payout
- The notebook builds the payout from the continuous stock, using `sample_stock` and `parametric_payout(s, strike) = max(strike - s, 0)`.
- The index model's payout is binary and drawn jointly with the mortality state (`index_model/insurance_contracts.jl:22-31`, `index_model/index_transitions.jl`). `El(s,p)` is the probability of a payout state, and it depends on `M` only.
- The capital cost therefore needs a new payout sampler matching the index contract: a Bernoulli payout with probability `El(M)`, correlated with the stock if that matters.
- Decide the state dependence. If the premium depends only on `M`, there are two values of `R'(Q*)` to precompute. If it depends on biomass too, tabulate on the state grid.
- Tabulate once per scenario at setup or `init` time. Do not call Monte Carlo inside `F`.

### B4. Places that must change
- Replace or wrap `price_per_exposure`:
  - `index_model/insurance_contracts.jl:36-40`
  - `src/insurnace_contracts.jl:29-33`
- Call sites:
  - `index_model/model.jl:153`
  - `src/model.jl:132`
  - `src/model_order_2.jl:132`
  - `index_model/analysis.jl:43`
  - `index_model/plotting.jl:56`
  - `src/plotting.jl:58`
  - `src/analysis.jl:117` and `:218`
  - `index_model/analysis/coverage.jl:134`
  - `index_model/analysis/run_bankruptcy.jl:181`
  - `index_model/index_transitions.jl:39,51` (`p_index`, `premium_ratio`)
- Parameters and scenarios:
  - `cf` and `cv` live in `parameters/universal_params.jl:16-17` and `parameters/RA_intermidiate.jl:38-39`.
  - `index_model/parameters/generate_scenarios.jl:52-56` defines `COST_VARIANTS` (`""` and `cv33`). Replace or extend with a capital-cost variant carrying `K`, `τ`, `α` and the `n`, `Q_base`, `CV` for the base portfolio.
  - Regenerate the `index_model/parameters/*_cv33.jl` files. The legacy `*_cv50.jl` files can be deleted if unused.
  - `index_model/grids.jl:295-296` writes `cf` and `cv` into `grid.toml`.
  - `index_model/analysis/setup.jl:48-49` applies the cf and cv overrides, and `:75-78` checks them (`assert_parameters_match`). Both need the new fields.
- Downstream code keyed on `cv`, `cv33` or the "0% versus 50% loading" label:
  - `index_model/analysis/run_composite_figure.jl:134-219`
  - `index_model/analysis/run_utility_figures.jl:37` and `:216-224`
  - `index_model/analysis/run_bankruptcy_figures.jl:236-296`
  - `index_model/analysis/figure_support.jl:61-67`
  - `index_model/run_no_insurance.jl:17` and `Makefile:138` (comments saying the no-insurance solution is invariant to `cf` and `cv`; still true, but reword)
  - `Makefile`: cheap and full run targets that name `cv33` scenarios
- `notebooks/insurer_model.ipynb` is a separate uncertainty-loading model (`λ` as loading for capital costs). It does not need to change, but its `λ` could be replaced by `R'(Q*)` if the two models are meant to agree.

## Tests

All tests are in `index_model/test/`. `test_coverage.jl` and `test_bankruptcy.jl` are not included in `make test_index`.

Add tests for the following:
- `fisher_response`:
  - FOC is about 0 at an interior solution.
  - The cap binds when `f*` exceeds the HCR.
  - Effort is 0 when the FOC is non-negative at zero effort.
  - Effort is 0 below `b_limit`.
- `dharvest_df` against a finite difference of `harvest`. The original had a wrong formula and an undefined `M`.
- The interpolated `f*(b)` against direct root finding.
- Capital cost:
  - `equilibrium_Q` residual is about 0.
  - `R'(Q*)` is non-negative and monotone in `K` and `τ`.
  - Premium equals `EL` plus the capital cost.
- Update `test_grids.jl` (metadata at lines 196-197) and `test_analysis.jl` for the new parameter fields.
- Consider adding `test_coverage.jl` and `test_bankruptcy.jl` to `make test_index`.

## Suggested order of work
1. Settle the open questions below.
2. Part A: port the fisher functions, tests, then the `f*(b)` interpolant, then swap into `update_stock`.
3. Recalibrate the parameters. Rerun a single small scenario end to end.
4. Part B: port the capital functions, tests, then the state-dependent premium table, then swap into `price_per_exposure`.
5. Update `generate_scenarios.jl`, regenerate the scenario files, and update analysis and figure code.
6. Rerun the full pipeline (`make` targets), regenerate figures and compare with the previous results.

## Verification
- `make test_index` (plus the two tests it omits) and the new unit tests.
- Run one small scenario end to end with both the new and old model, for example `no_RA_int_pr60_re100` or the cheap targets. Confirm the value function, simulation and bankruptcy table are sensible.
- Set `c` so that `f*` equals the old HCR effort. The new model should then reproduce the old results, which is a useful regression check.
- Check that the new capital premium with `K` very small approaches the fair price.

## Open questions (decide before coding)
1. **Scope.** Does the revision apply to the base model in `src/` and `analysis/` as well as `index_model/`, or only `index_model/`? The shared `src/` includes mean any `src/` edit changes both unless forked.
2. **Cap semantics.** The notebook caps *harvest* and then back-solves effort. Should the DP cap *effort* directly (`min(f*, f_allowed)`), which is simpler and monotone?
3. **Mortality in the decision.** Does the fisher choose `f` using expected mortality (`m̄`) or the realized mortality draw? This changes whether `f*` depends on `M`.
4. **Calibration.** Which of `f_target`, `c` and `p` stays a calibrated primitive? Today `μ_P` and `μ_Ce` are calibrated from `y_target` and `f_target`.
5. **Capital premium state dependence.** Does the capital cost depend on mortality state only or on biomass as well? The binary index payout differs from the notebook's continuous payout.
6. **Loading.** Does the capital cost replace the `cv33` multiplicative loading entirely, or do `cf` and `cv` stay for sensitivity analysis? Note the repo uses `cv = 0.33` for the "50% loading". The "0.3333" appears only in a comment.
7. **Legacy files.** Is the `*_order_2` code path still in use?
