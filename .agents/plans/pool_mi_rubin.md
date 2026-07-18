# Plan: Rubin's-rules pooling pathway for BACE (`with_imputations()` + `pool_mi()`)

Status: proposed (Track D). Owner: —. Created 2026-07-11.

## 1. Motivation

BACE currently has one combiner, `pool_posteriors()`, which stacks per-imputation
MCMCglmm posterior chains into a mixture (Monte Carlo marginal posterior; **not**
Rubin's rules). It is wired to BACE's internal `bace_final` object and works only
when the per-dataset model is the BACE MCMCglmm imputation chain itself.

We want a **model-agnostic, fit-then-combine** MI layer with a **free choice of
combiner**:

1. `with_imputations()` fits ANY user model to each of the `M` completed datasets —
   frequentist (`phylolm`, `pgls`, `gls`, `lm`, `glm`) **or** Bayesian (`MCMCglmm`).
2. Combine those fits either way:
   - `pool_posteriors()` — stack posterior draws (exact marginal posterior). For
     Bayesian fits and **mandatory for variance components**.
   - `pool_mi()` — Rubin's rules. For frequentist fits, **and** as an option on
     Bayesian fixed effects when the user wants standard MI reporting + `fmi`.

So the design is a full matrix: **any model x either combiner** (with guardrails),
not the pigauto split of "frequentist->Rubin, Bayesian->stacking". BACE is the
Bayesian package, so it is the natural home for "Bayesian fits, combined either way."

### Bayesian fits: stacking vs Rubin — which, and why offer both

For MCMCglmm fits you can do either, and they answer slightly different questions:

- **Stacking** (`pool_posteriors`) is the *coherent* Bayesian combiner: concatenate
  the posterior draws across imputations -> the exact marginal posterior integrating
  over the imputation distribution. No normality assumption. Handles skewed posteriors
  (variance components!) correctly. **This is the principled default.**
- **Rubin on posterior summaries** (`pool_mi`) extracts `Q_m = colMeans(Sol)` and
  `U_m = diag(cov(Sol))` per fit and applies Rubin's rules. This is a *normal
  approximation* to the stacked posterior. It is worth offering for three concrete
  reasons: (a) it yields **`fmi`** (fraction of missing information) — a diagnostic
  stacking does not naturally give, and the single most useful "should I collect more
  data?" number for a comparative biologist; (b) it produces a standard Rubin table
  that matches reviewer expectations and is comparable across frequentist and Bayesian
  analyses; (c) it needs only scalar summaries, so it is cheap.

The two agree as `M -> infinity` and differ by `O(1/M)` at small `M` (stacking uses
between-var weight `(M-1)/M`, Rubin uses `(1+1/M)`); at the shipped `n_final = 50`
default the gap is ~4%. **Never Rubin-pool variance components / lambda / heritability**
(skewed, bounded) — route those to stacking always.

### Zero-refit path (BACE's own imputation model)

`bace_final_imp()` already fits an MCMCglmm to each imputed dataset (`$all_models`).
So when the analysis model IS the BACE imputation chain, both combiners run on the
*same* already-fitted models with no refit: `pool_posteriors()` stacks them (as today),
and `pool_mi()` can pool their posterior summaries for the `fmi` table. Users only call
`with_imputations()` when they want a *different* downstream model than the imputation
chain.

### Relationship to pigauto (companion package)

pigauto's `pool_mi()` **rejects MCMCglmm and points to `BACE::pool_posteriors()`**.
BACE's `pool_mi()` deliberately does the opposite — it *accepts* MCMCglmm (auto-
extracting posterior mean/var) with a NOTE recommending stacking as the principled
default. This is coherent: pigauto punts the Bayesian case to BACE; BACE then offers
both combiners on Bayesian fits. Column names / df formula still mirror pigauto so the
frequentist tables are identical, and BACE's `imputed_datasets` still feed
`pigauto::pool_mi()` directly for users who prefer it.

## 2. API (mirror pigauto verbs)

Two new exported functions. Names match pigauto (`with_imputations`, `pool_mi`) for
cross-package consistency; masking when both packages are attached is harmless because
the semantics are identical and `::` disambiguates.

```r
# 1. Fit a frequentist analysis model on each of the M completed datasets.
fits <- with_imputations(result, function(dat, tree) {
  dat$species <- rownames(dat)
  phylolm::phylolm(y ~ x1 + x2, data = dat, phy = tree, model = "lambda")
}, tree = my_tree)
# result may be bace_complete / bace_final / a bare list of data.frames.
# -> class "bace_mi_fits": list of M fits (or bace_mi_error for failed draws)

# 2. Pool via Rubin's rules.
pool_mi(fits)
# -> class "bace_pooled_mi": data.frame(term, estimate, std.error, df,
#      statistic, p.value, conf.low, conf.high, fmi, riv)
```

### `with_imputations(object, .f, ..., tree = NULL, .progress, .on_error)`
- Accepts `bace_complete` (use `$imputed_datasets`), `bace_final` (use
  `$all_datasets`), or a plain list of data.frames.
- `.f = function(dataset, ...)` returns a fitted model. **Model-agnostic**: any class
  works — frequentist (`coef()`+`vcov()` for `pool_mi`) or `MCMCglmm` (carries `$Sol`
  for either combiner).
- If `.f` declares a `tree` argument, pass the tree through (single-tree case). This
  keeps the door open for a future `multi_impute_trees()` analogue without an API break.
- `.on_error = c("continue","stop")`: capture per-imputation errors, drop later.
- Returns `bace_mi_fits` (list) with a `print` method. Records whether the fits are
  MCMCglmm (so the two combiners can dispatch).

### `pool_mi(fits, conf.level = 0.95, coef_fun = NULL, vcov_fun = NULL, df_fun = NULL)`
- **Accepts MCMCglmm** (unlike pigauto). When fits are MCMCglmm and no extractors are
  given, default `coef_fun = function(f) colMeans(f$Sol)` and
  `vcov_fun = function(f) cov(f$Sol)` (posterior mean / posterior covariance of the
  fixed effects; drop the `Species.` BLUP columns first, as in `dev/10`). Emit a
  one-time NOTE: "Rubin on posterior summaries is a normal approximation; for the exact
  Bayesian combiner use pool_posteriors(). Do not use pool_mi for variance components."
  For non-MCMCglmm fits, default to `stats::coef` / `stats::vcov`.
- Extract coef vector + `diag(vcov)` per fit; align terms **by name** (error on
  mismatch, as pigauto does — Rubin needs a common term set).
- Apply Rubin's rules (§3). `df_fun` supplied -> Barnard-Rubin small-sample df
  (for MCMCglmm there is no residual df; default `df_fun = NULL` -> Rubin 1987 df).
- Return `bace_pooled_mi` (data.frame) with print method reporting `fmi` and the
  "increase M until SEs stabilise" hint when `max(fmi) > 0.5`.

### `pool_posteriors()` — generalize to accept `with_imputations()` output
- Today it takes a `bace_final` object. Extend it to also accept a `bace_mi_fits` list
  of MCMCglmm fits (stack their `$Sol`/`$VCV` row-wise, same as now). This makes the
  two combiners symmetric: the same `fits` object goes to either `pool_posteriors()`
  (stacking, exact; use for variance components) or `pool_mi()` (Rubin, fmi table).
- Error clearly if handed non-MCMCglmm fits (no posterior draws to stack) — the mirror
  of `pool_mi`'s behaviour, pointing the user back to `pool_mi()`.

## 3. Rubin's-rules math (validated line-for-line against pigauto::pool_mi and mice::pool)

Per term, over M imputations with per-fit estimate `Q_m` and squared SE `U_m`:

```
theta_bar = mean(Q_m)
W         = mean(U_m)                         # within-imputation
B         = var(Q_m)                          # between, (M-1) denominator
T         = W + (1 + 1/M) * B                 # total variance
SE        = sqrt(T)
r         = (1 + 1/M) * B / W                 # relative increase in variance (riv)
lambda    = (1 + 1/M) * B / T                 # = r/(1+r)
v_old     = (M - 1) * (1 + 1/r)^2             # Rubin (1987) df
# Barnard-Rubin (1999) if df_fun gives complete-data residual df nu_com:
v_obs     = ((nu_com + 1)/(nu_com + 3)) * nu_com * (1 - lambda)
v_bar     = 1 / (1/v_old + 1/v_obs)
fmi       = (r + 2/(v_bar + 3)) / (r + 1)
t_stat    = theta_bar / SE ; p = 2*pt(-|t_stat|, v_bar)
```

Guards (match pigauto): `W == 0 -> r = Inf, lambda = 1, fmi = 1`; `r = Inf ->
v_old = M-1`; `SE == 0 -> p = NA`. Use the Barnard-Rubin df because comparative
datasets often have small complete-data df where the Rubin-1987 df is badly biased.

## 4. Critical-reviewer caveats to build in

1. **Congeniality (Meng 1994, Stat Sci 9:538).** Rubin's variance can be biased when
   the analysis model is uncongenial with the imputation model. Fine for regression
   coefficients under BACE's rich imputation model; document, don't block.
2. **Do not pool variance components / lambda / heritability on the raw scale** — their
   sampling distributions are bounded/skewed. Provide optional `transform` /
   `back_transform` hooks (Fisher-z for correlations, log for variances) so those pool
   on an approximately-normal scale. Default identity.
3. **For Bayesian fits, stacking is the principled default; Rubin is an approximation
   offered for `fmi` + reporting.** `pool_mi()` on MCMCglmm uses posterior mean/var as
   `(Q_m, U_m)` — a normal approximation to the exact stacked posterior. Document this
   in the NOTE and recommend `pool_posteriors()` when the posterior is skewed. **Never
   Rubin-pool variance components / lambda / heritability** — route those to stacking.

## 5. Numerical-correctness testing (priority)

- **Golden test (strongest anchor):** feed identical `(Q_m, U_m, nu_com)` tuples to
  `BACE::pool_mi`, `pigauto::pool_mi` (if installed), and `mice::pool`; assert
  estimate/SE/df/fmi agree to machine tolerance. mice is the reference implementation
  of Barnard-Rubin df.
- **Edge cases:** `M = 1` errors; term-name mismatch errors with a useful message;
  `nu_com = Inf` reduces exactly to Rubin (1987); `W = 0` gives `fmi = 1`.
- **Integration (frequentist):** `bace()` on a small sim ->
  `with_imputations(gls/lm/phylolm)` -> `pool_mi()`; assert column set, finiteness,
  positive SE.
- **Bayesian fits:** (a) `pool_mi()` on a `bace_mi_fits` of MCMCglmm reproduces a
  hand-computed Rubin table from `colMeans(Sol)` / `diag(cov(Sol))`; (b) the same fits
  fed to `pool_posteriors()` (stacking) vs `pool_mi()` (Rubin) give fixed-effect point
  estimates and SEs agreeing to O(1/M) at `n_final = 50`; (c) a guard test that
  Rubin-pooling a variance component is flagged/steered to stacking.
- **Two-combiner validation (manuscript figure):** on the reference datasets
  (`dev/09-11`), show `pool_mi` (Rubin) and `pool_posteriors` (stacking) coverage of
  known beta agree at `n_final = 50`, diverging by the expected O(1/M) at small M.

## 6. Scope

**In v1 (guaranteed + tested downstream fitters):** `lm`, `glm`, `nlme::gls`,
`phylolm::phylolm`, `phylolm::phyloglm`. `caper::pgls` best-effort.

**Deferred to v2:** a `multi_impute_trees()` analogue (pool across a posterior sample
of trees as well as imputations; Nakagawa & de Villemereuil 2019, Syst Biol 68:632).
BACE takes a single tree today; `with_imputations()`'s `tree` pass-through is designed
so this can be added without an API break.

## 7. Files & integration

- New: `R/with_imputations.R`, `R/pool_mi.R` (+ print methods).
- Reuse: `bace_complete$imputed_datasets` / `bace_final$all_datasets` (already produced).
- NAMESPACE exports; roxygen -> `man/`; `NEWS.md` entry; README "Currently implemented".
- References: Rubin (1987); Barnard & Rubin (1999, Biometrika 86:948); van Buuren
  (2018, FIMD 2nd ed. Sec 2.3-2.4); Meng (1994); Nakagawa & Freckleton (2008 TREE,
  2011 BES) for the ecology MI framing.

## 8. Sequencing

Track D. Land **after Track A** (correctness hardening: stale `K=3` comments, explicit
`sample=FALSE` at bace_imp.R:352, oracle-prior alignment) so the pipeline is hardened
before adding surface. Reuses the Track B reference benchmark for the two-combiner
figure. Effort: Rubin math + `pool_mi` + golden test ~1 day; the `with_imputations`
model-extractor generality across model classes is the bulk.
```
