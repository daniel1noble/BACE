# Plan: BACE posterior predictive calibration

**Author**: Daniel Noble
**Date**: 2026-04-19
**Status**: For collaborator review before any code change

## TL;DR

While building a simulation-based imputation benchmark for BACE (`dev/simulation_imputation_quality.R`), we found that BACE's **95% posterior predictive intervals under-cover by 40–50 percentage points** on gaussian and poisson traits. The cause is not MCMC budget or number of imputations — it is a structural choice in `.predict_bace()` where `sample = TRUE` returns a draw of the **posterior mean** `E[y | X, u]` rather than a draw of the **posterior predictive** `y | X, u`. Missing the observation-level variance shrinks every predictive interval below its nominal width.

The scope is broader than "just gaussian and poisson". The same posterior-mean collapse happens (via a different mechanism) in **every non-gaussian, non-poisson branch** — binary, ordered, threshold, multinomial via OVR, and multinomial via probit — because the prediction helpers (`.pred_threshold`, `.pred_threshold_forward`, `.pred_cat`, `.pred_cat_forward`) all aggregate MCMC draws with `colMeans()` before any sampling happens. Since `ovr_categorical = TRUE` is the default, multinomial traits flow through `.pred_threshold_forward` and inherit this issue.

This matters because multiple imputation (Rubin 1987) requires imputations to be draws from the posterior predictive — not point estimates — for downstream variance estimates to be valid. The fix for the gaussian and poisson branches is small (~10 lines each). The fix for the categorical family is more invasive (return type changes in four helpers) but shares a single conceptual change: defer the averaging until after per-iteration sampling.

---

## 1. Problem statement

`bace()` returns `imputed_datasets`, a list of `n_final` imputed data frames. Each is drawn with `sample = TRUE` inside `.predict_bace()`. The standard interpretation — and the one assumed by any downstream analysis using Rubin's rules or credible-interval summaries — is that these are independent draws from the posterior predictive distribution of the missing data given the observed data.

If that interpretation holds, then for a continuous trait the sample interval

```r
[quantile(imputed_cell_i, 0.025),  quantile(imputed_cell_i, 0.975)]
```

should cover the true (hidden) value of cell `i` in approximately 95% of cells, up to finite-`n_final` noise. We set up a controlled simulation (`dev/simulation_imputation_quality.R`) to test this directly, because empirical datasets don't give us ground truth.

**Observed coverage on the demo replicate is 35–46%.** The size and direction of the miss are literature-consistent with the `.predict_bace` design choice below.

---

## 2. Evidence

### 2.1 Coverage does not improve with bigger MCMC

I re-ran the demo at two MCMC budgets — one minimal, one production-grade — holding everything else fixed:

| Setting                    | Small budget | Production budget |
|----------------------------|--------------|-------------------|
| `nitt`                     | 3000         | 50000             |
| `thin`                     | 5            | 25                |
| `burnin`                   | 1000         | 10000             |
| retained samples/formula   | ~400         | ~1600             |
| `runs` (convergence phase) | 4            | 10                |
| `n_final`                  | 6            | 20                |
| **y coverage95**           | 0.38         | **0.46**          |
| **x3 coverage95** (poisson)| 0.354        | **0.354**         |

With `n_final = 20` i.i.d. draws from the true posterior predictive, the sample-quantile coverage should hit roughly 90–95% (minor deflation below 95% from finite `n_final`, per standard order-statistics results — see van Buuren 2018 §2.5). Observed coverage at 35–46% is **40–50 percentage points below that**. Increasing MCMC budget 16.7× and imputations 3.3× moved y coverage by 8 pp and x3 coverage by **zero**. Monte Carlo noise alone cannot produce that pattern — the issue is structural.

### 2.2 Inspection of `.predict_bace()` explains the gap

`R/model_functions.R:521-547` (gaussian branch, `sample = TRUE`):

```r
eta     <- Sol[, common, drop = FALSE] %*% t(W[, common, drop = FALSE])
K       <- min(3L, nrow(eta))
i_samps <- sample.int(nrow(eta), K)
pred_values <- apply(eta[i_samps, , drop = FALSE], 2L, stats::median) * sd_val + mean_val
```

What this is: the median of 3 random posterior draws of `X·β + Z·u`. What it is **not**: a draw of `X·β + Z·u + ε` where `ε ~ N(0, σ²_units)`. The residual term is never added.

Consequences:
- The interval across `n_final` imputations reflects **parameter uncertainty only**. For models with adequate data, parameter uncertainty is far smaller than observation noise, so intervals are systematically too narrow.
- Median-of-3 further shrinks within-draw variance. This was a deliberate robustness measure against chain excursions, but it compounds the coverage problem.

`R/model_functions.R:550-561` (poisson branch, `sample = TRUE`):

```r
liab <- as.matrix(model$Liab)
i_samp <- sample.int(nrow(liab), 1L)
pred_values <- round(exp(as.numeric(liab[i_samp, ])), digits = 0)
```

This captures parameter uncertainty via the MCMC draw of `Liab`, but never samples the Poisson observation `Y ~ Poisson(exp(Liab))`. The predictive interval misses the entire Poisson sampling variance. That x3 coverage is *identical* at both MCMC budgets is the predicted signature of this issue: increasing MCMC cannot add variance the code never generated.

`R/model_functions.R:629-688` (`.pred_cat`, multinomial probit branch):

```r
prop_results[[i]] <- colMeans(exp_liab_list[[i]] / exp_sum)
```

The probability matrix is the **posterior mean** of class probabilities across MCMC iterations. `.impute_levels(..., sample = TRUE)` then draws a class from those mean probabilities. This is one averaging step too many: posterior uncertainty in the probabilities is smoothed out before sampling, so the `n_final` imputations cluster more tightly than the true posterior predictive.

`R/model_functions.R:842-904` (`.pred_threshold_forward`, used for binary/ordered/threshold AND for the J binary sub-models inside OVR):

```r
beta_mean <- colMeans(as.matrix(model$Sol))    # posterior MEAN of fixed effects
cp_all    <- c(0, colMeans(as.matrix(model$CP))) # posterior MEAN of cut-points
# ... species BLUPs also taken from beta_mean ...
eta       <- X_all %*% beta_k + species_BLUPs_mean
prob      <- pnorm(t_high - eta) - pnorm(t_low - eta)
```

This is the most severe collapse of the four helpers. Every probability is computed from posterior-mean coefficients, posterior-mean cut-points, AND posterior-mean species BLUPs. There is no MCMC iteration indexing at all. The returned probabilities are deterministic given the posterior mean of the model parameters. When `.impute_levels(sample = TRUE)` draws from these probabilities, every one of the `n_final` imputations samples from the *same* probability row per cell — the only variability across imputations is the multinomial draw itself, not any parameter uncertainty.

`R/model_functions.R:803-828` (`.pred_threshold`, the "non-forward" helper used when `formula + data_full` are not supplied):

```r
all_probs[[j]] <- colMeans(pnorm(t_high, mean = liab, sd = 1) -
                           pnorm(t_low,  mean = liab, sd = 1))
```

Slightly less collapsed than the forward version — it averages *after* the nonlinear `pnorm` — but still returns a deterministic probability row per cell after `colMeans`. Same downstream problem.

### 2.3 Why OVR inherits all of this

`ovr_categorical = TRUE` is the default in `bace()`. Under OVR, a k-class multinomial becomes k binary threshold models, each fit with MCMCglmm and each predicted via `.pred_threshold_forward`. The OVR aggregator at `R/model_functions.R:480-486` then row-normalises the k "yes" probabilities to form a multinomial probability row. Because every input to this aggregation is already a posterior-mean probability (from `.pred_threshold_forward`), the aggregated multinomial row inherits the collapse directly. So **multinomial-via-OVR is not insulated from the issue**; it compounds it by stacking k collapsed binary predictions.

Brier scores computed against these probability rows reflect the *point* probability estimate only, not the posterior uncertainty around it — so "good Brier" under the current code does not imply "well-calibrated class probabilities".

---

## 3. Why this conflicts with the multiple imputation literature

Multiple imputation (Rubin 1987) is built on the rule that each of the `m` imputations must be a **draw from the posterior predictive distribution** of the missing data given the observed data. Rubin's variance-combining formulas (the within-between variance decomposition) assume this directly; if the imputations are from the posterior of the mean rather than the posterior predictive, the between-imputation variance is too small and the pooled variance estimates understate uncertainty.

Direct statements to this effect:
- Rubin (1987, *Multiple Imputation for Nonresponse in Surveys*, §3) — the original formalism.
- van Buuren (2018, *Flexible Imputation of Data* 2e, §3.2, §5.1.5) — explicit: "Each imputation is a draw from the posterior predictive distribution of the missing data." §5.1.5 calls mean-based imputations "improper" and shows they produce invalid variance estimates.
- Little & Rubin (2019, *Statistical Analysis with Missing Data* 3e, §5.4) — proper MI coverage analysis requires posterior predictive draws.
- Gelman & Hill (2007, *Data Analysis Using Regression and Multilevel/Hierarchical Models*, §25.5) — in the context of `mi`, emphasises that imputations must include random error.
- Reference implementation `mice` (van Buuren & Groothuis-Oudshoorn 2011, *JSS* 45(3)): each elementary imputation method (`norm`, `logreg`, `polyreg`, `pmm`) explicitly draws parameters AND observation noise.

So the literature agreement is strong. The existing `.predict_bace` sampling path is what van Buuren calls "improper imputation".

---

## 4. Proposed fixes

### 4.1 Gaussian branch — REQUIRED

Change `.predict_bace()` so `sample = TRUE` returns proper posterior predictive draws:

```r
# Draw one MCMC iteration
i_samp <- sample.int(nrow(model$Sol), 1L)

# Linear predictor for this iteration
Sol_i <- model$Sol[i_samp, common, drop = FALSE]
eta_i <- as.numeric(Sol_i %*% t(W[, common, drop = FALSE]))

# Residual variance at the SAME iteration (critical for coherence)
# MCMCglmm stores residual variance in the 'units' VCV column.
sigma2_i <- model$VCV[i_samp, "units"]

# Posterior PREDICTIVE draw = linear predictor + residual noise
eps <- rnorm(length(eta_i), 0, sqrt(sigma2_i))
pred_values <- (eta_i + eps) * sd_val + mean_val
```

Keep the ±5σ safety clip (applied **after** the residual is added, so it only catches pathological chain excursions, not normal observation noise).

Drop the median-of-3 step. With proper per-iteration sampling, the occasional extreme draw is handled by the clip; the median-of-3 was a workaround for a problem (sampling `η` without `ε`) that goes away once residuals are included. If we're uncomfortable dropping the shrinkage entirely, a compromise is to median across K draws of `η_i + ε_i` — but that still shrinks residual variance and doesn't solve the coverage problem, so I'd argue against.

### 4.2 Poisson branch — REQUIRED

```r
i_samp <- sample.int(nrow(model$Liab), 1L)
rate   <- exp(as.numeric(model$Liab[i_samp, ]))
pred_values <- rpois(length(rate), lambda = rate)
```

The extra `rpois()` adds the Poisson sampling variance that the current `round(exp())` drops. No clip needed because `rpois` can't produce non-integer or negative values by construction.

### 4.3 Categorical-family branches — RECOMMENDED, larger surface, wide reach

**Affects:** binary (x1 in the benchmark), ordered / threshold (x4), multinomial via OVR (x2, default behaviour), multinomial via probit (non-default path). All four flow through helpers that collapse to posterior-mean probabilities before `.impute_levels(sample = TRUE)`.

Cleanest fix — one conceptual change, four helpers:

- `.pred_threshold_forward(model, ...)` currently returns an `n × k` posterior-mean probability matrix. Change it to return an `n × k × m` array where `m` is the number of retained MCMC iterations, built from per-iteration `Sol[i, ]`, `CP[i, ]`, and species BLUPs at iteration `i`. No cross-iteration averaging.
- `.pred_threshold(model, ...)` same: return per-iteration probabilities built from `Liab[i, ]` directly, without the outer `colMeans`.
- `.pred_cat(model, ...)` same: return per-iteration class probabilities from `Liab[i, ]` / `exp_sum[i]`, not `colMeans(exp_liab / exp_sum)`.
- `.pred_cat_forward(model, ...)` same.

`.impute_levels(pred_prob_array, levels_var, sample = TRUE)` then, per cell per requested draw:
1. Picks a random MCMC iteration `i`.
2. Reads the per-iteration probability row `pred_prob_array[cell, , i]`.
3. Samples a class with `sample(levels_var, 1, prob = row_i)`.

For OVR (`.fit_predict_ovr` at `R/model_functions.R:448-490`): each binary sub-model now returns a per-iteration probability array; the row-normalisation step aggregates them **at the same iteration** — so the J binary predictions are coupled by iteration, not by posterior mean. That preserves the joint MCMC uncertainty across the J models rather than collapsing each independently.

Why this is more invasive than (4.1) / (4.2):
- Return type changes in four helpers.
- `.impute_levels()` signature changes from matrix to 3-array input. It currently takes a 2D matrix; we'd either add a new arm or migrate wholesale.
- Downstream `prob_preds` object stored in `bace_final` and consumed by `pool_cat()` in the AVONET benchmark currently expects a cell × class matrix. Decision point: either `prob_preds` stays as a posterior-mean summary (computed ex post from the new 3-array), preserving downstream compatibility; or it migrates to the 3-array too.
- Memory: a per-iteration probability array is `n_cells × k × m`. At `n_cells = 500`, `k = 3`, `m = 1600`, that's 2.4M doubles per formula per run per imputation (~20 MB). With `n_final = 20` imputations it's manageable. For AVONET-scale (~11k rows), it's ~440 MB per formula per run — tractable but not negligible. If that's a concern, we can store only the per-cell posterior-mean row for `prob_preds` while retaining the per-iteration array for `sample = TRUE` draws internally.

Alternatives we are explicitly rejecting (for the record so collaborators can push back):

- **Jittering posterior-mean probabilities on the logit scale** to mimic iteration-to-iteration variability. Requires estimating a per-cell variance of `logit(prob)` from the MCMC output, then sampling `logit(prob_mean) + rnorm(0, sd)`. Workable but ad hoc; the full fix is not much harder and is what the MI literature actually recommends.
- **Doing nothing for categorical** on the grounds that point-estimate accuracy is fine. This would ship BACE with *known* miscalibrated class probabilities — which we document in `?bace_final_imp` as being posterior predictive. I'd argue against.

---

## 5. Backwards-compatibility impact

This change modifies the behaviour of `bace_final_imp()` (the function that calls `.predict_bace(sample = TRUE)`) and therefore `bace()`. Specifically:

- Imputed point values will have higher variance across the `n_final` datasets — this is the whole point.
- Point-estimate metrics (accuracy, NRMSE on the posterior mean, correlation) should change only marginally; they depend on the central tendency, not the spread.
- Any user analysis that treats individual `imputed_datasets[[i]]` as a deterministic best guess will see a shift (now genuinely stochastic). For categorical traits under OVR the shift will be largest, because currently `n_final` imputations are essentially re-draws from the same posterior-mean probability row; after the fix they are genuinely different posterior samples.
- `prob_preds` (stored on `bace_final` / `bace_complete`) currently contains posterior-mean class probabilities. Decision needed: either keep `prob_preds` as a posterior-mean summary (recomputed from the new per-iteration arrays, preserving the downstream signature used by `pool_cat()` in the AVONET benchmark), or migrate it to a per-iteration 3-array. My default is the former, with a new `prob_preds_array` slot added for users who want the full posterior.
- Pooled posterior models (`pool_posteriors`) are unaffected — they already stack full MCMC chains.

I think the change is a **correctness fix that supersedes the current design**, rather than a breaking change. But it's worth asking whether Szymek and Shinichi agree with that framing before we commit.

---

## 6. Test plan

Before and after each fix:

1. **Continuous coverage (gaussian, poisson)**: on a single simulated dataset with known truth, compute 95% PI coverage at (n_final = 20, production MCMC). Expected: 0.85–0.95 per variable after fix; ~0.35–0.45 before fix (current behaviour).
2. **Categorical probability calibration**: for binary / ordered / multinomial, compute per-class Brier and a reliability-diagram check (observed frequency vs predicted probability in probability bins; see Gneiting & Raftery 2007 for proper-scoring-rule treatment, or Dawid 1982 *JASA* 77:605 for the original "well-calibrated Bayesian" formalism). Expected after fix: Brier slightly higher (honest uncertainty), reliability diagram closer to the diagonal. Before fix: over-confident probabilities because posterior-mean rows have too little variance.
3. **Between-imputation variance is non-trivial for categorical traits**: for x2 (multinomial via OVR), measure `sd(imputed_class_i)` across `n_final` imputations for each hidden cell. Currently most cells will have sd = 0 (always the same imputed class because the probability row is deterministic). After fix, sd should be > 0 for cells whose posterior class probabilities are not concentrated on one level.
4. **AVONET benchmark regression**: re-run `avonet_mixed_benchmark.R`. Expected: point-estimate metrics within noise of current values; posterior pooling summaries unchanged.
5. **Point metrics on the simulation benchmark**: accuracy, NRMSE, correlation should not move more than ~2–3 pp. If they move more, something unintended happened.
6. **Convergence diagnostics** (`assess_convergence`): should be unaffected — the sampling change is only at the prediction-extraction step, not inside `.model_fit`.
7. **CRAN check**: `devtools::check()` passes with warnings = errors.

---

## 7. Decisions we need from collaborators

1. **Do we agree the current `sample = TRUE` behavior is a calibration bug under the Rubin MI framework?** (My reading: yes, per the literature in §3. Szymek / Shinichi — does that match your intent when the current path was written?)
2. **Scope of the fix**:
   - [ ] (4.1) Gaussian — required. Add residual draw to `.predict_bace` gaussian branch.
   - [ ] (4.2) Poisson — required. Replace `round(exp(Liab))` with `rpois(n, exp(Liab))`.
   - [ ] (4.3) Categorical family — affects **binary, ordered/threshold, multinomial via OVR (current default), and multinomial via probit**. One conceptual change applied to four helpers (`.pred_threshold`, `.pred_threshold_forward`, `.pred_cat`, `.pred_cat_forward`): return per-iteration probability arrays and sample at the iteration level in `.impute_levels`. `.fit_predict_ovr` aggregates by iteration, not by posterior mean. Do we:
     - [ ] Fix all four helpers in one pass, or
     - [ ] Ship (4.1)+(4.2) first, measure Brier under the current categorical code with those fixes in place, then decide about (4.3)?
   - [ ] (4.3 sub-decision) For `prob_preds`: keep as posterior-mean summary (recomputed from the new per-iteration arrays; downstream-compatible) or migrate to per-iteration 3-array (adds memory, more expressive)?
3. **Safety clip**: keep the ±5σ clip after adding residual noise, or drop it entirely now that pathological draws are less of a concern? My vote: keep, belt-and-braces is cheap.
4. **Version / release handling**: this is a behavior change for `bace()` users. Bump the version and add a `NEWS.md` entry before/after the change?
5. **Documentation**: `?bace_final_imp` and `?bace` currently describe the imputations as "posterior predictive" already (see CLAUDE.md). We'd be making the docs newly accurate.

---

## 8. Appendix: observed metrics on the production-MCMC demo run

Trait_MAR mechanism, n_species = 60, n_cases = 150, rate = 0.35, n_final = 20, nitt = 50000, thin = 25, burnin = 10000.

| Variable | Type          | NRMSE | Correlation | Accuracy | Coverage95 | Brier |
|----------|---------------|-------|-------------|----------|------------|-------|
| y        | gaussian      | 0.157 | 0.393       | —        | 0.460      | —     |
| x1       | binary        | —     | —           | 0.710    | —          | 0.392 |
| x2       | multinomial3  | —     | —           | 0.404    | —          | 0.791 |
| x3       | poisson       | 0.177 | 0.983       | —        | 0.354      | —     |
| x4       | threshold3    | —     | —           | 0.911    | —          | 0.169 |

The striking pattern is **high point-estimate accuracy (x3 correlation = 0.98, x4 accuracy = 91%) paired with low calibration (y coverage = 0.46, x3 coverage = 0.35)**. That pattern is exactly what the diagnosis predicts: point estimates are fine; predictive intervals are too narrow.

---

## 9. References

- Dawid, A.P. (1982). The well-calibrated Bayesian. *Journal of the American Statistical Association* 77: 605–610.
- Gelman, A. & Hill, J. (2007). *Data Analysis Using Regression and Multilevel/Hierarchical Models*. Cambridge University Press.
- Gelman, A. & Rubin, D.B. (1992). Inference from iterative simulation using multiple sequences. *Statistical Science* 7: 457–472.
- Gneiting, T. & Raftery, A.E. (2007). Strictly proper scoring rules, prediction, and estimation. *JASA* 102: 359–378.
- Graham, J.W., Olchowski, A.E. & Gilreath, T.D. (2007). How many imputations are really needed? *Prevention Science* 8: 206–213.
- Hadfield, J.D. (2010). MCMC methods for multi-response generalized linear mixed models: the MCMCglmm R package. *JSS* 33(2).
- Little, R.J.A. & Rubin, D.B. (2019). *Statistical Analysis with Missing Data* (3rd ed.). Wiley.
- Rubin, D.B. (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley.
- Stekhoven, D.J. & Bühlmann, P. (2012). MissForest — non-parametric missing value imputation for mixed-type data. *Bioinformatics* 28: 112–118.
- van Buuren, S. (2018). *Flexible Imputation of Data* (2nd ed.). Chapman & Hall/CRC.
- van Buuren, S. & Groothuis-Oudshoorn, K. (2011). mice: multivariate imputation by chained equations in R. *JSS* 45(3).
