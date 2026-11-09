# BACE 1.0.0 — defects found while building the 2026-11 simulation

Found 2026-11-09 by using BACE's own functions to build the manuscript's simulation. That is the
value of dogfooding: none of these would have shown up from the test suite, because every one of
them produces a plausible answer rather than an error.

Each entry says how it was established. **Verified** means reproduced directly here. **Reported**
means a survey agent found it and it is backed by a source reading, but it has not been
re-run independently — treat those as needing confirmation before acting.

---

## 1. A multinomial predictor loses its phylogenetic effect — VERIFIED — [issue #16](https://github.com/daniel1noble/BACE/issues/16)

**Where:** `R/simulate_simBACE.R:502-509`

By line 450, `linear_pred` holds the intercept **plus** `u_species`, `u_phylo` and the residual.
The multinomial-predictor branch then does:

```r
if (length(which_deps) == 0) {
  liab_linear <- liab_intercept                                    # throws linear_pred away
} else {
  liab_linear <- liab_intercept + linear_pred - intercepts$predictors[vari]
}
```

The no-dependency branch discards `linear_pred` entirely, so an independent multinomial predictor
gets **no random effects of any kind** — no phylogeny, no species effect, no residual.

**Evidence** (n = 150, 3-level multinomial predictor, `phylo_signal = 0.9`):

| case | var(u_phylo) generated | realised Pagel's lambda |
|---|---|---|
| no dependencies | 1.392 | **0.08** |
| depends on a gaussian trait | 1.392 | 0.34 |

The effect is generated and then not used. (In the dependent case the 0.34 largely comes *through*
the parent, not from its own effect: at `phylo_signal = 0` the same setup still gives 0.21.)

**Fix:** delete the `if` and always use the `else` expression. When there are no dependencies,
`linear_pred - intercepts$predictors[vari]` *is* the random effects, so the one expression is
correct in both cases.

**Scope:** multinomial *predictors* only. Threshold (`:518`), binary (`:487`) and poisson (`:490`)
all use `linear_pred` directly and are correct. A multinomial *response* is correct.

---

## 2. `pool_mi(df_fun = )` silently does nothing for every phylogenetic model — VERIFIED — [issue #17](https://github.com/daniel1noble/BACE/issues/17)

**Where:** `pool_mi()`, the `df_fun` block.

```r
if (!is.null(df_fun)) {
  nu_com <- tryCatch(df_fun(fits[[1]]), error = function(e) NA_real_)
  if (is.numeric(nu_com) && length(nu_com) == 1L && is.finite(nu_com) && nu_com > 0) {
    ...Barnard-Rubin...
  }                       # <- no else, no warning: falls back to uncapped Rubin 1987 df
}
```

`stats::df.residual()` returns **NULL** for `gls` (verified), and reportedly for `phylolm` and
`MCMCglmm` too. So the check fails and the small-sample correction is silently skipped.

**This is in the documentation:** `?pool_mi`'s own example passes `df_fun = stats::df.residual`.
A user following it gets byte-identical output to passing nothing, with no indication.

**Fix:** warn when `df_fun` is supplied but returns something unusable, naming the class it was
called on. Optionally ship a helper that knows the df of `gls`/`phylolm` fits, and correct the
example.

---

## 3. `with_imputations()` silently drops the tree unless the argument is literally named `tree` — VERIFIED (source) — [issue #18](https://github.com/daniel1noble/BACE/issues/18)

**Where:** `with_imputations()`:

```r
f_formals <- names(formals(.f))
if (!is.null(tree) && "tree" %in% f_formals && !("tree" %in% ...)) call_args[["tree"]] <- ...
```

If `.f`'s phylogeny argument is named anything else (`phy`, `phylo`) **and has a default**, the
tree the user passed is dropped and the default is used — no warning. A phylogenetic model then
silently conditions on the wrong tree. Reported coefficient: 0.4407 (wrong tree) vs 0.9785
(correct). Without a default you at least get an error.

**Fix:** if `tree` is supplied and `.f` has no formal named `tree`, warn.

---

## 4. Barnard-Rubin df diverges from `mice` at zero between-imputation variance — REPORTED

`v_old <- ifelse(is.finite(r) & r > 0, (M-1)*(1+1/r)^2, M-1)`. The correct limit as `r -> 0` is
`Inf`, not `M-1`; `mice::barnard.rubin()` floors lambda at 1e-4 so df tends to `dfcom`. Reported
effect: df 3.73 where mice gives 56.09, a CI 1.43x too wide, and `fmi` reported as 0.297 when it
should be ~0. Matters whenever an imputation contributes no between-imputation variance.

---

## 5. Smaller things — REPORTED, worth confirming

- **`beta_resp` as a named numeric vector is matched positionally**, names discarded, so
  `beta_resp = c(x2 = 2.0, x1 = 0.0)` assigns 2.0 to `x1`. A named *list* is matched correctly.
  Silent mis-assignment. Either honour the names or reject a named vector.
- **Poisson columns come back as `double`** from `get_imputed_data()`, so `.get_type()` reclassifies
  them as gaussian on any second pass — a silent family switch.
- **Column order changes**: the cluster column comes back last, so completed data is not
  `identical()` to the input. Never index completed data by position.
- **`sigmas` cannot be scalars** — `list(sigma_species = 1, sigma_residual = 1)` fails with
  `'names' attribute [3] must be the same length as the vector [1]`. Needs a clear error.
- **`n_species` is an upper bound, not the tip count** (species are sampled with replacement), so
  `Ntip(tree)` can be less than asked for. Worth documenting.
- **Fitted objects are ~35 MB** at modest size (`$pooled_models` 14.9 MB, `$final_results` 18 MB).
  Fine for one analysis; a campaign must extract statistics in the worker and never save the fit.

---

## 6. Softmax overflow in `mnom_liab2cat()` — VERIFIED — [issue #19](https://github.com/daniel1noble/BACE/issues/19)

`probs <- exp(liab_full) / rowSums(exp(liab_full))` in `R/simulate_auxiliary.r:196`. A liability
over ~710 makes `exp()` return `Inf`, and `Inf/Inf` is `NaN`, which `sample(prob = )` rejects:
`"NA in probability vector"`. Killed 3 of 54 scenarios outright, all at n = 500 with high signal
and related traits. One-line fix: subtract the row maximum before exponentiating (softmax is
shift-invariant, so this changes nothing mathematically).

## 7. Unguarded `coda::geweke.diag()` aborts `bace()` — VERIFIED — [issue #20](https://github.com/daniel1noble/BACE/issues/20)

`R/model_functions.R:1635` calls `coda::geweke.diag(param)` with no `tryCatch`, so a parameter
whose draws have no variance takes down the whole run with `"zero-variance series"` — after the
expensive fitting is already done. The package wraps the identical call at
`R/phylo_signal_summary.R:850`, so the hazard is known; it is just missing here. Killed 2 of 54
scenarios. Deterministic, so resuming does not help.

---

## Not a bug, but it must be said in the manuscript

`phylo_signal` is phylogenetic heritability **on the latent scale**:

```
                          sigma2_phylo
    phylo_signal  =  ------------------------------------------------
                     sigma2_phylo + sigma2_species + sigma2_residual
```

For a **gaussian** trait that equals Pagel's lambda, and `sim_bace` delivers it accurately —
measured 0.29 / 0.70 / 0.91 against nominal 0.3 / 0.7 / 0.9, unchanged by whether the traits are
related to each other.

For **threshold and count** traits the observed-scale lambda is always lower, because thresholding
and Poisson sampling add variance that was never in the linear predictor. This is the classic
liability-scale vs observed-scale heritability distinction (Dempster & Lerner 1950, *Genetics*
35:212-236; Nakagawa & Schielzeth 2010, *Biol Rev* 85:935-956). Measured with BACE's default
sigmas, at nominal 0.3 / 0.7 / 0.9:

| trait | realised observed-scale lambda |
|---|---|
| gaussian | 0.23 / 0.61 / 0.86 |
| poisson | 0.06 / 0.29 / 0.62 |
| ordinal (4 level) | 0.08 / 0.45 / 0.79 |
| binary | **0.00 / 0.04 / 0.18** |

A binary trait at nominal 0.3 has effectively no phylogenetic signal at all. Two things follow:

1. The simulation must **specify signal on the latent scale and report realised observed-scale
   lambda per dataset** as a measured result, not assert it.
2. The default `sigmas` make it worse than it has to be: 1.0 for gaussian variables but
   sqrt(0.1) = 0.316 for every non-gaussian one, a ~10x smaller latent variance. Raising them
   helps the discrete traits, but raising them carelessly hurts the response, because the
   fixed-effect variance contributed by upstream traits is **not** in the denominator above.
   Any choice here needs calibrating and checking, not guessing.
