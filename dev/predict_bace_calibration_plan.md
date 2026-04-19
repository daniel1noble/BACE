# Plan: BACE posterior predictive calibration

**Author**: Daniel Noble
**Date**: 2026-04-19
**Status**: For collaborator review before any code change

## TL;DR

While building a simulation-based imputation benchmark for BACE (`dev/simulation_imputation_quality.R`), we found that BACE's **95% posterior predictive intervals under-cover by 40–50 percentage points** on gaussian and poisson traits. The cause is not MCMC budget or number of imputations — it is a structural choice in `.predict_bace()` where `sample = TRUE` returns a draw of the **posterior mean** `E[y | X, u]` rather than a draw of the **posterior predictive** `y | X, u`. Missing the observation-level variance (residual for gaussian; Poisson sampling noise; per-iteration class probabilities for categorical) shrinks every predictive interval below its nominal width.

This matters because multiple imputation (Rubin 1987) requires imputations to be draws from the posterior predictive — not point estimates — for downstream variance estimates to be valid. The fix is small (~10 lines per branch in `R/model_functions.R`) but it changes the behaviour of `sample = TRUE`, so it deserves discussion before we commit.

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

`R/model_functions.R:629-688` (`.pred_cat`, which `sample = TRUE` feeds):

```r
prop_results[[i]] <- colMeans(exp_liab_list[[i]] / exp_sum)
```

The probability matrix is the **posterior mean** of class probabilities across MCMC iterations. `.impute_levels(..., sample = TRUE)` then draws a class from those mean probabilities. This is one averaging step too many: posterior uncertainty in the probabilities is smoothed out before sampling, so the `n_final` imputations cluster more tightly than the true posterior predictive. Brier will be biased downward (optimistic) for the same reason.

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

### 4.3 Categorical / threshold branches — RECOMMENDED but with larger surface

Currently `.pred_cat()` / `.pred_threshold()` return `colMeans` of per-iteration class probabilities. `.impute_levels(sample = TRUE)` then draws from those mean probabilities.

Cleanest fix: have `.pred_cat()` / `.pred_threshold()` return the per-iteration probability array (dims: cell × class × iteration). `.impute_levels()` then picks a random iteration per draw and samples a class from that iteration's probability row.

This is more invasive than (4.1) and (4.2) because:
- `.pred_cat` / `.pred_threshold` return type changes.
- Downstream `prob_preds` storage and `pool_cat` in the benchmark would need updating.

A minimal-change alternative: keep the current return type but within `.impute_levels`, do a Monte Carlo draw that approximates per-iteration sampling by the formula

```r
probs_mean <- current posterior-mean probability row
# add uncertainty by jittering on the logit scale with sd inferred from MCMC var
```

This is hackier and I'd argue against. If we decide calibration matters for categorical at all, do the full fix.

---

## 5. Backwards-compatibility impact

This change modifies the behaviour of `bace_final_imp()` (the function that calls `.predict_bace(sample = TRUE)`) and therefore `bace()`. Specifically:

- Imputed point values will have slightly higher variance across the `n_final` datasets — this is the whole point.
- Point-estimate metrics (accuracy, NRMSE on the posterior mean, correlation) should change only marginally; they depend on the central tendency, not the spread.
- Any user analysis that treats individual `imputed_datasets[[i]]` as a deterministic best guess will see a small shift (now genuinely stochastic).
- Pooled posterior models (`pool_posteriors`) are unaffected — they already stack full MCMC chains.

I think the change is a **correctness fix that supersedes the current design**, rather than a breaking change. But it's worth asking whether Szymek and Shinichi agree with that framing before we commit.

---

## 6. Test plan

Before and after each fix:

1. **Unit-ish test**: on a single simulated dataset with known truth, compute coverage at (n_final = 20, production MCMC). Expected: 0.85–0.95 per continuous variable after fix; ~0.35–0.45 before fix (current behavior).
2. **AVONET benchmark regression**: re-run `avonet_mixed_benchmark.R`. Expected: point-estimate metrics within noise of current values; posterior pooling summaries unchanged.
3. **Point metrics on the simulation benchmark**: accuracy, NRMSE, correlation should not move more than ~2–3 pp. If they move more, something unintended happened.
4. **Convergence diagnostics** (`assess_convergence`): should be unaffected — the sampling change is only at the prediction-extraction step, not inside `.model_fit`.
5. **CRAN check**: `devtools::check()` passes with warnings = errors.

---

## 7. Decisions we need from collaborators

1. **Do we agree the current `sample = TRUE` behavior is a calibration bug under the Rubin MI framework?** (My reading: yes, per the literature in §3. Szymek / Shinichi — does that match your intent when the current path was written?)
2. **Scope of the fix**:
   - [ ] (4.1) Gaussian — required
   - [ ] (4.2) Poisson — required
   - [ ] (4.3) Categorical / threshold — recommended but more disruptive. Do we fix all three in one pass, or gate (4.3) on measuring Brier behavior after (4.1)+(4.2)?
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
