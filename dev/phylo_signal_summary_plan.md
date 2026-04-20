# Plan: `phylo_signal_summary()` for BACE

## Context

BACE imputes missing values in phylogenetic comparative datasets using per-variable phylogenetic mixed models (MCMCglmm). Imputation quality depends on how much phylogenetic information each variable carries — yet BACE currently gives users no way to check this before committing to a full multi-hour run. Users who invoke BACE on a dataset with weak phylogenetic signal and weak covariate coupling may be surprised by low imputation accuracy without understanding why.

The goal is a diagnostic tool (`phylo_signal_summary()`) that fits fast univariate phylogenetic models and reports per-variable phylogenetic heritability H² (plus classical λ / K / D metrics where applicable). A `phylo_signal = TRUE` argument to `bace()` runs the diagnostic and returns it without executing the full imputation loop, so users can inspect the signal table and decide whether to proceed, adjust priors, drop variables, or abandon the run. Critically, H² formulas must track BACE's single-RE vs dual-RE configuration (`species = FALSE` vs `species = TRUE`), and the formulas must be validated by simulation — not merely asserted correct.

**Scope note (reviewer caveat to include in docs):** phylogenetic signal is *necessary but not sufficient* for good imputation. A high-H² trait can still impute poorly under MNAR missingness; a low-H² trait can still impute well if tightly coupled to an observed predictor. The table sets expectations; it does not predict performance.

## Design

### New function

```r
phylo_signal_summary(
  data,
  tree,
  species_col     = "Species",
  variables       = NULL,         # default: all columns except species_col
  species         = FALSE,        # mirrors bace(): FALSE = single RE, TRUE = dual RE
  methods         = "auto",       # "auto" picks metrics by type; or c("H2","lambda","K","D")
  ovr_categorical = FALSE,        # per-level OVR H² for multinomial variables
  nitt            = NULL,         # NULL = use type-specific defaults (see below)
  burnin          = NULL,
  thin            = NULL,
  quick           = FALSE,        # halve nitt and burnin for rough first-pass exploration
  prior           = NULL,         # NULL -> .make_prior() defaults
  n_sim           = 999,          # for K p-value
  min_ess         = 1000,         # flag fits with effective size < this on any variance component
  keep_models     = FALSE,
  verbose         = TRUE
)
```

Returns an S3 `phylo_signal` object with:
- `$table` — data.frame, one row per variable, columns: `variable`, `type`, `n_obs`, `H2_mean`, `H2_lo`, `H2_hi`, `R_species_mean`, `R_species_lo`, `R_species_hi` (populated only when `species = TRUE`), `lambda`, `lambda_p`, `K`, `K_p`, `D`, `D_random_p`, `D_BM_p`, `interpretation`, `flag` (e.g. `"low_n"` when n_obs < 20).
- `$models` — named list of fitted MCMCglmm objects (only if `keep_models = TRUE`).
- `$call`, `$species` (single vs dual RE), `$n_sim`.

### Base models fit per variable

For each variable `v` in `variables`, fit one univariate MCMCglmm with `v` as response and an intercept only as fixed effect. The random-effect structure, family, prior, and MCMC defaults vary by variable type and by `species`. All models use `ginverse = list(Species = A)` built from the tree via `MCMCglmm::inverseA(tree)`. Dual-RE fits require a `Species2` identity-ginverse copy column (created internally, not persisted).

**Random-effect structure:**

| `species` | MCMCglmm `random` | G-structure blocks |
|---|---|---|
| `FALSE` | `~ Species` | one (phylo) |
| `TRUE` | `~ Species + Species2` | two (phylo + non-phylo species) |

`Species2` is a string-identical copy column — MCMCglmm needs distinct term names to fit two species-level REs; identity covariance comes from no ginverse entry for `Species2`.

**Family, residual prior, and G prior by type** (reusing `.make_prior()` / `.list_of_G()` from [R/model_functions.R:299](R/model_functions.R#L299)):

| Variable type | `family` | `R` prior (residual) | `G` prior (each block) |
|---|---|---|---|
| gaussian | `"gaussian"` | `list(V = 1, nu = 2)` | `list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25²)` (par_expand = TRUE) |
| poisson | `"poisson"` | `list(V = 1, nu = 2)` | same as gaussian |
| binary | `"threshold"` | `list(V = 1, fix = 1)` | same as gaussian |
| ordered (≥3 levels) | `"threshold"` | `list(V = 1, fix = 1)` | same as gaussian |
| multinomial (K levels) | `"categorical"` | `list(V = (I_{K-1} + J_{K-1})/K, fix = 1)` (IJ fix) | `list(V = diag(K-1), nu = K-1, alpha.mu = rep(0, K-1), alpha.V = diag(K-1) * 25²)` (par_expand) |
| OVR multinomial (per level) | `"threshold"` (J fits per variable) | `list(V = 1, fix = 1)` | same as gaussian |

These match what BACE itself uses per variable — the diagnostic is literally "what H² does BACE's own per-variable model see?" — which is why the answer is directly relevant. Also `slice = TRUE` is set automatically for threshold/ordered/categorical inside `.model_fit()`.

**MCMC defaults — type-specific (the reviewer concern):**

Hadfield's MCMCglmm course notes (§3 "Categorical Traits", §4 "Running a Chain") flag that variance components for threshold and multinomial models mix far slower than gaussian, and that H² is disproportionately sensitive to variance-component autocorrelation. If `nitt`/`burnin`/`thin` are `NULL`, use type-specific defaults balanced to deliver **1 500 post-burn draws** for every type, with thinning scaled to the mixing speed so effective sample size on each variance component targets **ESS > 1 000**:

| Type | nitt | burnin | thin | post-burn draws | burnin fraction |
|---|---|---|---|---|---|
| gaussian | 20 000 | 5 000 | 10 | 1 500 | 25 % |
| poisson | 30 000 | 7 500 | 15 | 1 500 | 25 % |
| binary / ordered | 50 000 | 12 500 | 25 | 1 500 | 25 % |
| multinomial (K level) | 80 000 | 20 000 | 40 | 1 500 | 25 % |

Higher thin for slower-mixing models reduces autocorrelation between saved draws, so ESS approaches 1 500 when variance components mix adequately post-burn. Gaussian chains typically achieve ESS near 1 500 (near-independent draws); threshold / multinomial may sit at ESS 1 000–1 400 on a typical dataset — still above the flag threshold. The 25 % burnin fraction matches Hadfield's rule of thumb.

A single-species gaussian fit (~60 species) completes in ~30 s on a modern laptop; multinomial ~2–5 min. For a 10-variable dataset with mixed types this is ~10–20 min of diagnostic compute — still one to two orders of magnitude faster than a full BACE run, and firmly in "cheap relative to what it prevents".

**Quick mode (`quick = TRUE`):** for a rough first-pass exploration, halve `nitt` and `burnin` (thin unchanged). Post-burn draws drop to ~750 and ESS typically sits around 500 on slow-mixing components — so the `low_ess` flag will trigger more often, which is the intended signal to the user that quick mode is just a preview. Emit a single warning at function start naming this.

| Type | quick `nitt` | quick `burnin` | `thin` | post-burn draws |
|---|---|---|---|---|
| gaussian | 10 000 | 2 500 | 10 | 750 |
| poisson | 15 000 | 3 750 | 15 | 750 |
| binary / ordered | 25 000 | 6 250 | 25 | 750 |
| multinomial | 40 000 | 10 000 | 40 | 750 |

**Tree-size warnings (adaptive, at function start):**

The defaults above are calibrated for the typical comparative-dataset regime (50–500 tips). Both extremes need user awareness:

- **`n_species < 30` → "small_tree" warning.** Phylogenetic variance σ²_phylo is near-unidentifiable when few tips are available; the posterior H² collapses toward the prior mean regardless of true signal. A tiny-tree fit with `sim_bace(phylo_signal = 0.2, n_species = 15, ...)` will return a posterior mean of H² that drifts back toward ~0.5 with a 95 % HPD spanning roughly (0.05, 0.93). The ESS flag will NOT catch this — it's a prior-sensitivity problem, not a mixing problem. Warning text: "Tree has <30 species. Phylogenetic variance is weakly identified; H² posterior is prior-sensitive. Interpret with caution, and consider the HPD width more than the posterior mean."

- **`n_species > 1000` → "large_tree" warning.** Wall clock scales roughly linearly with `n_species` (phylogenetic kernel evaluation cost). A multinomial fit on 3 000 tips at default `nitt = 80 000` takes 30–60 min per variable. Five multinomial variables becomes a 3–5 h diagnostic. Warning text: "Tree has >1000 species. Default MCMC may underrun; consider increasing nitt by √(n_species/500) ≈ {computed factor} or using quick=TRUE for a rough first pass."

Both warnings are informational, not blocking — the function still runs. They appear in the printed output header and are also attached to `$warnings` on the returned object so users can programmatically inspect them.

**Convergence diagnostic per fit (mandatory, not optional):**

After each fit compute `coda::effectiveSize()` on `V_A`, `V_S` (if dual), `V_R`, and each diagonal entry for multinomial. Store per-fit `min_ess` in the output table. If `min_ess < min_ess` threshold (default 1 000), set the row's `flag` column to `"low_ess"` and emit a warning naming the variable. Users should rerun that variable with higher nitt. Do NOT silently return biased H² estimates from unconverged fits — this is the failure mode the user is worried about and the table must surface it.

Also compute `coda::geweke.diag()` on each variance component and flag `"geweke_fail"` if |z| > 2 on any component (Hadfield course notes §4). A fit that fails either diagnostic gets `flag = "unreliable"` and its H² cell gets an asterisk in `print()`.

### H² formulas (all latent-scale to match `sim_bace()` ground truth)

Let `V_A = σ²_phylo` (VCV column named after phylo ginverse key), `V_S = σ²_species` (VCV column for the Species2 identity-ginverse RE, dual-RE only), `V_R = σ²_residual` (VCV `units`).

| Type | Single RE (`species = FALSE`) | Dual RE (`species = TRUE`) |
|---|---|---|
| gaussian | `V_A / (V_A + V_R)` | `V_A / (V_A + V_S + V_R)` |
| poisson (log link, latent scale) | `V_A / (V_A + V_R)` | `V_A / (V_A + V_S + V_R)` |
| binary / ordered threshold (probit, R fixed = 1) | `V_A / (V_A + 1)` | `V_A / (V_A + V_S + 1)` |
| multinomial (K levels, R = IJ fix, par_expand = TRUE) | `tr(G_phylo) / (tr(G_phylo) + tr(R))` | `tr(G_phylo) / (tr(G_phylo) + tr(G_species) + tr(R))` |

Per-draw computation over posterior samples of VCV → posterior mean + 95% HPD. For `species = TRUE`, also report `R_species = V_S / total`.

**Why trace-based for multinomial:** under MCMCglmm's `R = IJ, fix = 1` parameterization for multinomial probit (Hadfield course-notes §5), the K-1 latent liabilities are exchangeable under the same phylo process. Trace H² is the direct generalization of the scalar case and aligns with sim_bace's generating model (K-1 independent liabilities sharing the same σ²_phylo). Per-dim averaging inflates uncertainty; determinant-based is near-singular for K=3. References: de Villemereuil et al. 2016 *Genetics* 204:1281; Nakagawa & Schielzeth 2013 *Methods Ecol Evol* 4:133; Hadfield & Nakagawa 2010 *J Evol Biol* 23:494.

**par_expand invariance:** under parameter expansion, the rescaling multiplier cancels between numerator and denominator in trace-H², so the formula is invariant (Gelman 2006).

**OVR path (`ovr_categorical = TRUE`):** fit J binary threshold models per multinomial variable (mirroring BACE's OVR path), report a per-level H² column and a mean H² across levels.

### Classical metrics (auxiliary columns)

- λ (Pagel 1999 *Nature* 401:877) — continuous only, via `phytools::phylosig(method = "lambda", test = TRUE)`
- K (Blomberg, Garland & Ives 2003 *Evolution* 57:717) — continuous only, via `phytools::phylosig(method = "K", test = TRUE, nsim = n_sim)`
- D (Fritz & Purvis 2010 *Conservation Biology* 24:1042) — binary only, via `caper::phylo.d()`

Populated as `NA` for types where the metric is inapplicable. Skip with a message if the package is missing.

### `bace()` integration

New argument `phylo_signal = FALSE` in `bace()` [R/bace.R:94](R/bace.R#L94), slotted after `species`. When `TRUE`:

1. Call `phylo_signal_summary(data, tree, species_col, variables = all RHS + LHS, species = species, ...)`.
2. Print the table.
3. Return the `phylo_signal` object with a top-level class attribute `c("phylo_signal", "bace_preview")`.
4. Emit a message: `"phylo_signal = TRUE: returning signal-only preview. Rerun with phylo_signal = FALSE to proceed with imputation."`
5. Do **not** invoke `bace_imp()` / `bace_final_imp()` / `pool_posteriors()`.

### `print.phylo_signal` method

Match existing `print.bace_complete` style (`cat` + `sprintf`, hierarchical `===`/`---` headers). Include legend explaining metrics, interpretation bands (H² < 0.2 low, 0.2–0.5 moderate, > 0.5 high), a caveat footer stating signal is necessary-but-not-sufficient, and a `low_n` flag for variables with n_obs < 20.

## Files

### Create

- [R/phylo_signal_summary.R](R/phylo_signal_summary.R) — main function + internal helpers `.compute_H2_gaussian`, `.compute_H2_threshold`, `.compute_H2_multinomial`, `.compute_H2_poisson`, `.compute_lambda_K`, `.compute_D`.
- [R/print.phylo_signal.R](R/print.phylo_signal.R) — S3 print method (or co-located in the main file).
- [tests/testthat/test-phylo_signal_summary.R](tests/testthat/test-phylo_signal_summary.R) — fast smoke test (see validation plan below).
- [dev/validate_phylo_signal_summary.R](dev/validate_phylo_signal_summary.R) — full Monte Carlo validation grid, manually run.
- [man/phylo_signal_summary.Rd](man/phylo_signal_summary.Rd) — roxygen-generated.

### Modify

- [R/bace.R:94](R/bace.R#L94) — add `phylo_signal = FALSE` argument; short-circuit branch at function top after input validation.
- [DESCRIPTION](DESCRIPTION) — add `phytools`, `caper` to `Suggests`. Soft-fail with a message if missing when classical metrics are requested.
- [NAMESPACE](NAMESPACE) — regenerated via `devtools::document()` to export `phylo_signal_summary` + `print.phylo_signal`.
- [CLAUDE.md](CLAUDE.md) — add a short architecture entry under "Key internal functions" pointing at the new file.

### Reuse (do not recreate)

- `.make_prior()` / `.list_of_G()` from [R/model_functions.R:315](R/model_functions.R#L315) — use with `n_rand = 1` or `2` depending on `species`.
- `.detect_type()` from [R/prep_functions.R](R/prep_functions.R) — per-variable type detection.
- `.data_prep()` conventions for the Species2 copy column — see [R/bace_imp.R:298](R/bace_imp.R#L298).
- Ginverse construction mirroring [R/model_functions.R:73](R/model_functions.R#L73).
- `sim_bace()` from [R/simulate_simBACE.R:1](R/simulate_simBACE.R#L1) — generator with known `phylo_signal` ground truth for validation.

## Simulation validation plan

Ground truth in `sim_bace()` (confirmed at [R/simulate_simBACE.R:291](R/simulate_simBACE.R#L291)): `phylo_signal = p` sets `σ²_phylo = p · (σ²_species + σ²_residual) / (1 − p)`, giving latent-scale H² = p exactly. Our posterior H² formula is on the same latent scale by construction, so recovery is apples-to-apples.

### Two-tier strategy

**Fast smoke (`tests/testthat/test-phylo_signal_summary.R`, `skip_on_cran`, target <6 min):** 4 cells, 1 replicate each, `n_species = 50`, `n_cases = 150`. MCMC settings use the **function's type-specific defaults** (not reduced) so the smoke test also confirms that defaults actually deliver adequate ESS — catches the failure mode the user flagged (low-nitt → biased H²):

| # | type | RE | phylo_signal | σ²_species | nitt/burnin/thin |
|---|---|---|---|---|---|
| 1 | gaussian | single | 0.6 | 0 | 20 000 / 5 000 / 10 |
| 2 | binary | single | 0.6 | 0 | 50 000 / 12 500 / 25 |
| 3 | categorical-3 | single | 0.6 | 0 | 80 000 / 20 000 / 40 |
| 4 | gaussian | dual | 0.6 | 1 (3 reps/sp) | 20 000 / 5 000 / 10 |

Per-cell assertions:
- Posterior mean of H² within ±0.20 of truth.
- 95% HPD contains truth.
- `min_ess ≥ 1 000` on every variance component (matches function default).
- Table shape, column types, and `print()` run without error.

Test 3 (categorical) is the single expensive cell (~3–4 min). Gate it behind `skip_if(Sys.getenv("BACE_SKIP_SLOW") == "true")` so developers can skip locally but CI still runs it. A separate ultra-fast structural test (no MCMC recovery check, `nitt = 2000`, just asserts class/columns/print) runs in ~10 s on every `devtools::test()`.

**Full validation (`dev/validate_phylo_signal_summary.R`, manual, ~4–6 h):** 4 types × 2 RE × 2 signals × B = 20 reps = 320 fits at user's proposed settings (`n_species = 80`, `n_cases = 240`, `replicates_per_species = 3` for dual, `nitt = 13000, burnin = 3000, thin = 10`).

| types | RE | signals |
|---|---|---|
| gaussian, binary, categorical-3, poisson | {single, dual} | {0.2, 0.7} |

Pass criteria:
- |mean(posterior mean − truth)| ≤ 0.05 per cell
- 95% HPD coverage ≥ 0.90 per cell
- RMSE across replicates ≤ 0.15

Prune `threshold-3` (redundant with binary under scalar-latent treatment). Keep categorical-3 (only trace-aggregation test). Dual-RE cells REQUIRE ≥ 3 replicates per species — with 1 replicate, σ²_species is confounded with σ²_residual and the Species2 posterior collapses toward the prior.

Output: `dev/benchmark_results/phylo_signal_validation/run_<tag>/` with a `summary.csv` (bias, coverage, RMSE per cell) and a `validation_report.md`. Gate: if any cell fails, investigate and fix the formula before merging.

## Dependencies

- **New `Suggests`:** `phytools` (λ, K), `caper` (D). Soft-fail with informative message if missing when a classical metric is requested. H² column (the headline) depends only on already-present `MCMCglmm` + `ape`.

## Verification

1. `devtools::document()` — regenerate NAMESPACE + man pages.
2. `devtools::test(filter = "phylo_signal")` — runs fast smoke in <6 min. All cells pass.
3. `devtools::check()` — 0 errors, 0 warnings, existing notes only.
4. `Rscript dev/validate_phylo_signal_summary.R` — full Monte Carlo validation. Inspect `summary.csv`; all cells meet pass criteria.
5. Manual end-to-end: on AVONET-style test data, run `bace(..., phylo_signal = TRUE)` and confirm the table prints, the return object has class `phylo_signal`, and no full BACE run executes.
6. Regression: rerun `dev/00_benchmark_AVONET.R` and `dev/01_benchmark_simulated.R`; confirm prior benchmark metrics are unchanged.

## Open questions (not yet resolved)

1. **ESS threshold for `low_ess` flag** — 200 (permissive), 500 (strict), or two-tier yellow/red (500/200)?
2. **Document tree-size limitations** — defaults assume 50–500 species; wider trees mix worse, smaller trees have near-unidentified variance components. Worth a dedicated caveats section in roxygen?
3. **Opt-in fast mode** — 5 multinomial variables × 80 k iterations ≈ 15 min wall clock. Expose a `quick = TRUE` that halves nitt with a loud warning that flags will be triggered more often?
