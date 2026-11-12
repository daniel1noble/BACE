# BACE 1.0.3

* **The trait simulator produced different data on different machines, from the
  same seed.** `sample_random_effects()` drew the phylogenetic effects with
  `MASS::mvrnorm()`, which builds its matrix square root from
  `eigen(Sigma, symmetric = TRUE)`. Eigenvectors are **not unique**: LAPACK
  fixes no sign convention, so any implementation may return a column negated,
  and for repeated eigenvalues the whole subspace is arbitrary. Every such
  choice is a valid square root, so each machine drew a different -- though
  correctly distributed -- set of traits from an identical random stream.

  Measured across three machines running the same seed, tree and RNG settings:

  | | `x1[1:3]` | realised lambda |
  |---|---|---|
  | Accelerate (macOS) | 7.169, -1.485, 0.273 | 0.827 |
  | OpenBLAS (Linux) | 2.047, -1.537, 2.636 | 0.875 |
  | AOCL (Linux) | 4.074, -1.417, 4.599 | 0.850 |

  The trees, tip labels, `runif()`, `sample()` and the masked cells were
  **identical** on all three, which is what isolated the fault to the matrix
  square root rather than to the random-number stream. Negating half the
  eigenvector columns on one machine reproduces it exactly.

  The draw now goes through the **Cholesky factor**, which is unique for a
  positive-definite matrix, so implementations agree to rounding error instead
  of choosing between valid alternatives. It still consumes exactly `n` normals,
  so trees, missingness patterns and everything else downstream are unchanged.
  A near-degenerate matrix gets a deterministic ridge, never an eigen fallback.

  Recovery is unaffected: imputation accuracy against known truth is
  statistically identical before and after (mean r 0.510 vs 0.528 over 30
  seeds, p = 0.85). The simulated VALUES all change, because the old ones were
  machine-dependent and there was nothing to preserve.

  `tests/testthat/test-simulator-reproducibility.R` pins this, including a
  golden test with a fixed tree and seed at a tolerance of 1e-8 -- loose enough
  for rounding, tight enough to catch a sign flip. **Running the test suite on
  any machine is therefore a cross-machine check.**

* `test-numerical_correctness.R`'s imputation-recovery test ran on 40 species,
  leaving 12 hidden cells. A correlation on 12 points is far too noisy to
  assert anything: over 20 seeds it ranged from -0.29 to +0.93 and fell below
  its own threshold in 7 of them, so it was passing on a lucky seed rather than
  on a property of the code. At 120 species the same measurement is stable
  (mean 0.724, sd 0.153, minimum 0.406, none below threshold).

# BACE 1.0.2

* **`plot_steps_convergence()` drew bars on a different yardstick from the
  one in its own caption.** The panel divided each step by
  `mean(sd_df[[var]])`. For a numeric variable that is the imputation spread
  and is correct; for a PROPORTION (binary, ordinal, categorical) `sd_df`
  holds `sqrt(p(1-p))`, a function of the mean, which the check explicitly
  refuses to treat as a spread -- it measures proportions in proportion units.
  A settled 40-cell binary was drawn at 0.050 where the panel header printed
  0.025, and a near-certain one (p = 0.999) at 0.092 against a true 0.005, an
  **18.4-fold** overstatement that put a settled variable above the red alarm
  line. Its fallback, `stats::sd(m[keep])`, was the circular yardstick removed
  from the check itself earlier in this version.

  `assess_convergence()` now reports the yardstick each verdict was reached
  with, as `spread` in the per-variable statistics, and the panel uses it
  instead of deriving its own. A convergence object from an older version
  records none; the panel falls back as before, or says so rather than
  drawing a misleading bar.

* **A third audit, of the 2026-11-10 fixes themselves, found six defects that
  those fixes had introduced.** All are corrected here. Recorded because the
  pattern matters: each round of fixes has created defects of its own, and the
  only thing that caught them was auditing the fix rather than trusting it.

  - **A dotted response COLUMN NAME made every categorical imputation collapse
    onto the baseline level.** Earlier on 2026-11-10 the level-label strip was
    changed from greedy (`sub(".*\\.", "")`, which loses a LEVEL containing a
    dot) to non-greedy (`sub("^trait.*?\\.", "")`). But MCMCglmm names traits
    `trait<response>.<level>`, so the non-greedy rule loses a RESPONSE NAME
    containing a dot -- and those are common, because `read.csv()` and
    `data.frame()` turn every space in a column heading into a dot. Renaming a
    column from `diettype` to `diet.type`, nothing else changed, moved the
    imputed spread from 11/18/7 to **36/0/0**: the two non-baseline levels
    became un-imputable, accuracy against the hidden truth fell from 65 per
    cent to 2.5 per cent, and the between-imputation variance became exactly
    zero, so any downstream interval ignored missing-data uncertainty. Neither
    loose rule is right in general; the strip is now the exact, regex-escaped
    `^trait<response>.` taken from the formula's LHS, which is correct for a
    dotted response, a dotted level, and both at once.

  - **A colon in a factor level made heritability collapse to zero.** The
    multinomial diagonal matcher split the column name on `":"` and required
    exactly two pieces, which assumes no colon inside a level name. Renaming
    one level from `b` to `b:1` moved `H2_mean` from 0.988 ("high") to
    **0.000** ("low"), with no warning. Diagonals are now identified by shape
    -- the same trait identifier either side of one colon -- however many
    colons the level names contain.

  - **The new pooling unit-correction never fired for a categorical response,
    nor for any interaction term.** It matched a Sol column only when the whole
    name was exactly a z-scored variable, so it saw `x1` but not `x1:x2` and
    not `trait<resp>.<level>:x1` -- i.e. it did nothing at all for multinomial
    slopes. Measured at 8.6 to 23.6 per cent on multinomial interval widths,
    against the ~4 per cent that justified the correction in the first place.
    A column's units now come from every standardised variable named in it.

  - **`get_scaling()`'s own worked example used the wrong divisor**, and
    contradicted the `@param model` text three lines above it: it took the
    predictor's constants from its own model (the observed-only SD) rather
    than from the model the coefficient belongs to. Measured at up to 16 per
    cent on the recovered slope, systematically upward. The `@return` field
    also described a two-level structure where the object is three-level
    (run -> model -> variable).

  - **Known limitation, deliberately not approximated:** the pooling correction
    harmonises the SCALE but not the CENTRE, so the pooled intercept still
    mixes runs centred on different means (about 4.6 per cent on its interval
    width). The exact correction is easy with main effects only but an
    interaction changes the main effects too, so a partial application would be
    wrong rather than incomplete. The slopes, which are what is normally
    interpreted, are exact.

* **`pool_posteriors()` stacked runs that were in different units.** Each run
  standardises its PREDICTORS by the spread of its own completed data, so run
  m's slope is "per sd_m of x" and run m+1's is "per sd_m+1 of x"; stacking
  them combines numbers measured with different rulers. (The response's divisor
  is stable -- its missing rows are re-masked before scaling, so it is the
  observed-only sd every run -- so the drift is entirely in the predictors.)
  Each run's coefficients are now rescaled to a common reference divisor before
  stacking, which leaves the overall scale of the returned object unchanged and
  only touches columns whose name exactly matches a z-scored variable; factor
  levels and random-effect BLUPs are not standardised and are left alone.

  Measured at 50 per cent missingness, where `sd(x1)` ranged 16.5 to 23.5 over
  15 runs: against the reference of converting each run with its own constants
  and then stacking, the pooled estimate was out by 0.37 per cent and the
  interval width by 4.2 per cent; both are now 0.00 per cent. The November
  audit put this at 8.6 per cent and 10 per cent, which did not reproduce --
  the effect is real but small, and it does not grow with the drift, because
  the per-run divisors scatter around their mean and largely average out.
  A `bace_final` from an earlier version carries no constants and pools
  uncorrected, as before.

* **The z-scoring constants are now kept, and there is `get_scaling()` to read
  them.** BACE standardises continuous variables before fitting, so the fitted
  models -- `$pooled_models` and the per-run models -- carry coefficients in SD
  units, not data units. The mean and SD used were computed, applied and then
  **discarded**, so no exact back-transform was possible from BACE's output at
  all: a true slope of 2.0 came back as 0.785 with nothing recording the
  divisor.

  They are recorded per imputation run **and per model**, because that is how
  they are computed. A variable's divisor is not one number: as the response of
  its own model it is the observed-only SD (stable), but as a *predictor* in
  another variable's model it is recomputed from that run's imputed values and
  moves run to run -- measured at 21.6 to 24.9 against an observed-only 22.9, a
  spread of 8.5 per cent. A coefficient is in units of both sides, so
  converting the slope of `y` on `x1` needs
  `get_scaling(fit, "x1", model = "y")`, not `x1`'s own. With the recorded
  constants the example above back-transforms to 2.0048 against a truth of 2.0;
  the nearest guess available before was out by 3.3 per cent.

  Checkpoints written by earlier versions carry no constants, so they are now
  rejected by the shape check and recomputed rather than silently resumed with
  the constants missing.

* **`phylo_signal_summary()` reported heritability of 1.0000 where the truth
  was 0.50.** The probit latent scale has no size of its own, so the residual
  is FIXED to identify it, and the H2 formula wrote that `1` in by hand:
  `V_A / (V_A + V_S + 1)`. A user prior that pins the scale elsewhere scales
  everything on it, `V_A` included, so the numerator scaled and the denominator
  did not. Measured on one 100-tip binary fit with a true H2 of 0.50: pinning
  the residual at 1 / 2 / 4 gave `V_A` 12.0 / 21.4 / 41.2 and a reported H2 of
  0.810 / 0.884 / 0.929 against the model's own 0.810 / 0.805 / 0.795. Worse,
  the standard MCMCglmm idiom `R = list(V = 1, nu = 0.002)` leaves the residual
  **free**: it ran to 2e19 and H2 was reported as **1.0000**. `bace()` already
  refuses an R override for these types; this function exposed the same
  documented `prior` argument with no such guard. It now reads the residual the
  model actually used (identical results under the default prior) and refuses a
  free residual, naming the fix. `R_species` had the same hard-coded `1`.

* **`phylo_signal_summary()` could fit a different model family from
  `bace()` for the same column.** It built its own type summary with a looser,
  NA-sensitive integer rule — `is.integer(x) || (is.numeric(x) && all(!is.na(x)
  & x == floor(x)))` against `bace()`'s plain `is.integer(x)`. A whole-number
  **double**, which is what `read.csv()` gives for a count column, was gaussian
  to `bace()` and **poisson** here, so the same data got a different variance
  decomposition: H2 0.492 against 0.766-0.787 on a column whose true value was
  0.8, flipping the printed band from "high" to "moderate". And because
  `all(!is.na(x) & ...)` is FALSE as soon as one value is missing, a single
  `NA` silently changed the family back. Both paths now use
  `.summarise_var_types()`, so there is one rule.

* **Random slopes in `ran_phylo_form` are now an error, not a silent
  substitution.** BACE fits random-INTERCEPT models only. The formula parser
  could read a slope, but the prior builder cannot construct a correct one:
  it widens `V` from 1x1 to k x k and widens `alpha.mu`/`alpha.V` to match
  while leaving `nu` at 0.002, which is improper for k > 1 (an inverse-Wishart
  on a k x k matrix needs `nu > k - 1`; Hadfield's course notes use `nu = k`
  for a k-dimensional parameter-expanded prior). The dimension itself is
  sniffed by a regex that only recognises `us(...)`, so `idh(1 + x1 + x2)`
  fell through to a hardcoded 2. Worse, on the **default** `species = FALSE`
  path the slope was discarded outright: `~ 1 + x1 | Species` and
  `~ 1 | Species` produced byte-identical models, so the user asked for one
  thing and silently got another. Supplying a slope on either path now errors,
  naming what was given and the supported form. Proper support is a feature
  request ([#22](https://github.com/daniel1noble/BACE/issues/22)), not a fix.

* **A declared-but-unobserved factor level shifted every imputation onto the
  wrong category.** When a level of an ordered or categorical response is never
  seen in training, MCMCglmm re-indexes the response to `1:K_eff` in
  declaration order, so the prediction comes back with fewer columns than the
  factor has levels. The right vocabulary is then the **observed** levels in
  declaration order. `.impute_levels()` passed no observed set, so
  `.align_levels_to_probs()` fell back to *the first `n_cols` declared levels*
  — correct only when the unobserved levels happen to sit at the end. With
  levels `2 | 3 | 4` and only `3` and `4` observed, it returned `2 | 3`, moving
  every imputed value one category down, reported only as a warning naming no
  variable (`'(unknown)'`). The observed set and the response name are now
  passed through, and the other call site already did this correctly. Reached
  in practice: the 2026-11 campaign's ordinal trait realises three declared
  levels at low phylogenetic signal, with the top one observed once or not at
  all.

* `assess_convergence()`'s reported `thresholds` are derived from the code's
  own defaults rather than written out by hand, so they cannot drift from the
  values actually applied again.

* This is the version the 2026-11 simulation campaign requires. The guards in
  `run_one.R` and `run_all.sh` refuse to start below it, because an older build
  reporting the same version number produces results that look entirely normal
  and mean nothing.

# BACE 1.0.1

* **The chained-equations convergence check was defective and has been
  replaced.** The old check compared each variable's movement between runs to
  the variable's own *level* (`%` change) and its own mean (the energy and
  Wasserstein variants), so its verdict tracked where a variable sat rather
  than whether it had settled: a mean-centred trait could *never* pass
  (division by ~0), while a trait far from zero passed even while steadily
  drifting. Fed chains that are converged by construction, the old check
  passed ~3% of them at default settings; fed steady drifts at level 20, it
  passed essentially all of them. In practice this meant `bace()` silently
  retried (5 → 8 → 11 runs, ~4.8x the compute) on nearly every real dataset
  and then warned that convergence had failed when it had not. Separately,
  a chain whose imputation variance collapsed to zero — the textbook
  non-convergence of van Buuren & Groothuis-Oudshoorn (2011, p. 37) — was
  special-cased as *converged*.

  The new check (`assess_convergence()`, method `"stability"`) watches the
  mean *and* the SD of each variable's imputed values across chained runs,
  treats the first half as warm-up, measures all movement in units of the
  variable's imputation spread, and flags three named pathologies: **drift**
  (still travelling; the one thing more runs can cure — `bace()` now retries
  only for this), **oscillation** (variables flip-flopping in step; retrying
  is pointless and `bace()` now says so and proceeds), and **stuck**
  (collapsing imputation variance) and **diverging** (the spread growing
  without bound). Categorical and ordinal variables are monitored on every
  level and take the worst level's verdict — watching only the first level
  is blind to a chain that flips its cells between the others.

  Calibration on constructed chains (pinned in the test suite): settled
  chains pass >= 99% at `runs = 10` wherever their mean sits and however
  many cells are imputed; drift, loud oscillation, collapse and divergence
  are each caught >= 99%; false oscillation stays under 5% even at 1000
  imputed cells and the most pessimistic between-run noise we have measured.
  Short chains say what they could not do: below `runs = 5` the drift test
  cannot run and below `runs = 7` the flip-flop test cannot, and the output
  carries a note saying so. A pass means *no pathology detected*, not a
  certificate — no clear-cut convergence test exists for chained equations
  (van Buuren 2018, ch. 6.5.2; Oberman, van Buuren & Vink 2021,
  arXiv:2110.11951) — so the trace plots remain the primary diagnostic.

  Breaking changes, all loud: the per-criterion methods
  (`"summary.acf"`, `"summary.percentage"`, `"summary.trend"`,
  `"summary.geweke"`) and the distribution methods (`"energy"`,
  `"wasserstein"`, `"all"`) are defunct with explanatory errors
  (`"summary"` itself is a deprecated alias for `"stability"`); the
  `pct_change_threshold`, `lag` and `alpha` arguments are gone (replaced by
  `drift_threshold` and `stuck_ratio`); `plot()` type `"pct_change"` is
  replaced by `"steps"` (movement against the Monte-Carlo noise floor) and
  the energy/Wasserstein plot helpers are defunct. Fits made with earlier
  versions re-assess cleanly with the new check.

* **`with_imputations()` no longer drops the phylogeny in silence** (#18). The
  tree was forwarded only when `.f` declared a formal literally named `tree`.
  A user who wrote `function(d, phy = my_tree)` — the natural thing to write
  with a tree in scope — got their *default* tree in every fit, so every
  per-imputation model conditioned on the wrong phylogeny with no error and no
  warning. On one such pair the reported coefficient was 0.4407 with the tree
  dropped against 0.9785 with it passed, and nothing downstream looked wrong.
  Supplying a tree that `.f` cannot receive is now an error that names `.f`'s
  actual arguments and how to fix the call; a function taking `...` now
  receives the tree too.

* **The stability check called settled variables "drift" most of the time**
  (#20). When the imputed cells had no within-run spread, `.conv_assess_series()`
  fell back to `sd()` of the very series it was judging. That makes
  `step_mean = mean|diff(m)| / sd(m)` scale-free by construction: for any iid
  wobble, however small, it converges to `2/sqrt(pi) = 1.128`, against a
  `wander` gate of 1.00. A perfectly settled variable was therefore reported
  as drift **81-89 per cent of the time**, and more runs made it worse (0.80 at
  `runs = 10`, 0.88 at `runs = 50`), because the statistic is scale-free. Two
  ordinary situations land on that branch: a variable with exactly **one**
  missing cell (the SD of a length-1 vector is `NA`, recorded as 0), and a
  **poisson** variable whose rounded point imputations all land on the same
  integer. Each false verdict costs `bace()` a full extra chained phase and
  still ends `converged = FALSE`. The yardstick is now the spread of the
  variable's **observed** values — a scale independent of the series under
  test — and when no such scale exists the verdict is `"unassessable"` rather
  than an invented one. Settled series now return `"ok"` 100 per cent of the
  time at every `runs`, a genuine drift of 3 data-SD is still caught 100 per
  cent of the time, and series with a real within-run spread are unchanged.
  This was the same defect the 1.0.1 rewrite set out to fix — a derived
  quantity used as if it were an independent measurement — reintroduced in the
  fix itself.

* **`assess_convergence()` reported thresholds it was not using** (#20). The
  returned `thresholds` vector was written out by hand in two places and had
  drifted from the code: `osc_floor` was reported as 0.50 while 0.60 was
  applied. It is now derived from `.conv_assess_series()`'s own defaults, so
  the two cannot diverge again, and the numeric-only entries are named as such
  (`drift_threshold_numeric`, `osc_floor_numeric`) — proportion series are
  judged on a tightened, per-variable bar reported in each assessment's `ref`.
  The `wander` gate is also now documented as unreachable for proportion
  series (a per-run step cannot exceed 1), which is a detection gap, not a
  wrong number: drift and alternation are still covered by the other tests.

* **The poisson rate guard no longer truncates in silence** (#20). The cap
  (`3 * max(observed) + 5`) exists to stop log-link chain excursions — an
  `exp(20.4)` liability with a maximum observed count of 2 — becoming
  astronomical imputations, and it stays. But it cannot distinguish an
  excursion from genuine extrapolation, so under MAR-on-high (the missing
  counts really being larger than the observed ones, which is the regime
  multiple imputation exists for) it truncates legitimate posterior-predictive
  draws. Measured on an 80-species pipeline with the 20 largest-covariate
  counts hidden: the summed cell means came out 8.3 per cent below truth and
  the between-imputation SD shrank to 0.33x its uncapped value on the worst
  cells, with nothing reported. BACE now warns when the cap actually binds,
  naming the variable, the cap and how many rates were truncated, so the
  trade-off is visible rather than discovered.

* **A categorical level whose name contains a dot could never be imputed**
  (#20). `.pred_cat_forward_core()` labelled the predicted-probability columns
  by taking everything after the *last* dot of MCMCglmm's trait name, so with
  levels `forest | forest.edge | open` the trait `traithabitat.forest.edge`
  was labelled `edge`. `.predict_bace(sample = TRUE)` builds its sampling
  vocabulary with `intersect(levels_var, colnames(probs))`, so that column was
  dropped and the remaining mass renormalised: across 210 posterior-predictive
  draws `forest.edge` was imputed **0 times** although it was 35 per cent of
  the observed data. The probabilities themselves were always correct — only
  the label was wrong, which is why nothing failed. Labels are now stripped by
  the same non-greedy `trait<response>.` prefix rule `.pred_cat()` already
  used, and a column name that is not a declared level now raises a warning
  instead of silently deleting a category. Dotted and plain level names now
  give identical imputations (68/28/44 in both). The `sample = FALSE` path was
  unaffected, since `.impute_levels()` falls back to positional alignment.

* **Multinomial `H2` and `R_species` were inflated** (#20).
  `.get_multinom_diag_cols()` matched the residual diagonal with a pattern
  requiring a *numeric* trait index, but MCMCglmm labels traits by factor
  level (`traity.B`), never by number. The match therefore always failed on
  real fits and a silent fallback took the first J columns — which in
  MCMCglmm's column-major order is column 1 of R: one diagonal entry plus
  J−1 **off-diagonal** ones. Since R is fixed at `(I + 11')/K`, that sum is
  exactly 1 for every K instead of the true `tr(R) = 2(K-1)/K`, so the
  residual trace was understated by 25 per cent at K=3 and 40 per cent at K=6,
  and reported heritability grew too large with the number of categories. The
  package's own test fabricated digit-labelled names, so it passed over a
  broken path. Diagonals are now matched by splitting on `:` and comparing the
  two trait identifiers, which also copes with dotted level names, and the
  no-columns fallback is `2J/K` rather than `J²/K` (equal only at K=3).

* **The reported gaussian posterior SD was not a standard deviation** (#20).
  The back-transform multiplied the whole `.pred_cont()` data frame by `sd`
  and added `mean`, so the location shift reached `post_sd` as well:
  `sd(a + bX) = |b| sd(X)` has no `a` term. On a variable with mean 1000 and
  sd 0.01 the reported SD was 1000.003 instead of 0.0031 — out by a factor of
  3e5 — and it was **negative** for any variable with a negative mean. Imputed
  values were never affected (only column 1 is consumed), but `post_sd` is
  returned to the user. Columns are now transformed individually.

* **`pool_mi()` reported missing information where there was none** (#20).
  `B = 0` (every imputation agreeing exactly, which happens whenever the
  analysis model uses only fully observed variables) and `W = 0` share a
  fallback, but their limits differ: as `r → ∞` the df tend to `M-1`, while as
  `r → 0` they tend to **infinity**. Sending both to `M-1` made
  `fmi = (r + 2/(v+3))/(r+1)` report `2/(M+2)` of the information missing when
  none was. BACE now returns `df = Inf`, `fmi = 0`, `riv = 0`, matching the
  limit `mice` reaches by clamping lambda (df ≈ 9e8, fmi ≈ 2e-9).

* **A character vector of formulas silently lost every response but the
  first** (#20). `fixformula = c("y ~ x1", "z ~ x1")` is not a list, so it took
  the single-formula branch, where `stats::as.formula()` keeps only the first
  element (R warns only that `formula()` is deprecated for length > 1) and
  auto-expands from it. `z` disappeared entirely, and the guard that every
  variable with missing data has a formula was skipped on that same path, so
  nothing complained. A character vector of several formulas is now treated as
  the list the user meant, in `bace()`, `bace_imp()` and `bace_final_imp()`.

* **`with_imputations()` now pairs each fit with the tree that actually
  generated its imputation** (#18). The pairing was re-derived here as plain
  modular recycling, but `bace_final_imp()` assigns trees with
  `.select_tree_indices()`, which *spreads* across the supplied trees when
  there are more trees than imputations. With 100 posterior trees and 50
  imputations the imputations were generated under trees 1, 3, 5, ..., 99
  while the downstream fits used trees 1, 2, 3, ..., 50 — so 49 of 50 models
  conditioned on a tree that had not generated their data. The old guard
  compared only the tree *counts*, so passing the same 100 trees (the correct
  thing to do) matched on count and warned about nothing. The pairing is now
  read from the object's recorded `$tree_index`, with the supplied trees
  checked against it by content hash rather than by count; a set that does not
  match falls back to `bace_final_imp()`'s own rule and says so. Runs with no
  more trees than imputations are unaffected — the two rules always agreed
  there.

* **`pool_mi(df_fun = )` no longer skips the small-sample correction in
  silence** (#17). `stats::df.residual()` returns `NULL` for `gls` (and for
  `phylolm` / `MCMCglmm`), and an unusable value fell through without a
  signal, so `pool_mi(fits, df_fun = stats::df.residual)` on the models
  phylogenetic users actually fit produced output byte-identical to applying
  no correction at all. The reported degrees of freedom were the uncapped
  Rubin (1987) values, which can exceed the complete-data df absurdly: on a
  30-tip `gls` the intercept was reported with df = 1.2e8 where the correct
  Barnard & Rubin (1999) value is 26.2, and the slope with df = 8.5 against
  4.2 — making nominal 95% intervals about 94% and 92%. BACE now derives the
  complete-data df itself (`nobs - n_coef`, verified equal to
  `df.residual()` where both are defined) and says so, or warns loudly when
  no value can be derived.

* **`sim_bace()`: a multinomial predictor with no dependencies kept none of
  its random effects** (#16). The no-dependency branch discarded the linear
  predictor and with it the phylogenetic, species and residual effects, so
  such a trait had essentially no phylogenetic signal however `phylo_signal`
  was set, silently. Both branches now use the same expression: liability SD
  across species goes from 0.000 to 1.000, and simulated traits respond to
  `phylo_signal` again (observed-scale lambda 0.002 at 0.10 vs 0.317 at 0.95).
  Multinomial *responses* were never affected.

* Fixed a crash in `sim_bace()` when a multinomial predictor's liabilities
  were large (softmax overflow; GitHub issue #19).

* Fixed a crash in the final-model diagnostics when a chain was too short or
  degenerate for `effectiveSize()`/`geweke.diag()`/`acf()`; each is now
  reported as NA instead of aborting the fit (GitHub issue #20).

# BACE 1.0.0

First CRAN release. Earlier versions were distributed from GitHub only. The
1.0.0 number marks a stable user-facing API, not a change in method. Every
change below is additive, with one small exception noted under the per-phase
MCMC settings (`fixformula` can no longer be abbreviated as far as `f`/`fi`).

* **Long runs can now be resumed instead of restarted.** Give `bace()` (or
  `bace_imp()` / `bace_final_imp()`) a `checkpoint_dir` and the convergence
  chain is saved when it finishes, and each final imputation run as soon as it
  completes. Calling again with the same arguments and the same directory
  loads that work and computes only what is missing — so a session that dies
  during the final runs picks up where it left off rather than starting over.
  (The convergence chain is saved once it finishes, so a crash during that
  phase still restarts it.) A resumed run gives
  exactly the same answer as one that was never interrupted, which is why
  `seed` is required. If any setting that would change the results differs
  from the saved work, the resume is refused and the differing settings are
  named, so results from two different analyses are never mixed. Default
  `NULL`, i.e. off, and nothing else changes.

* **The final imputation runs can now have their own MCMC settings.**
  `bace()` gains `final_nitt`, `final_thin` and `final_burnin` (default
  `NULL` = reuse the chained-phase settings, exactly as before). The two
  phases do different jobs — the chained phase only needs point imputations
  to settle, while the final runs supply the draws that are pooled for
  inference — so they often want very different chains, for example
  `nitt = 50000, thin = 40` with `final_nitt = 5e5, final_thin = 450` (scale
  the thinning with the chain: every retained draw is stored for all
  `n_final` runs). Previously this required bypassing `bace()` and calling
  `bace_imp()` and `bace_final_imp()` separately, which is what the bird case
  study had to do. (The amphibian study splits the phases for a different
  reason — memory on an 18,000-row dataset — which this does not change.)
  The new arguments sit at the end of `bace()`'s argument list, so existing
  calls are unaffected; the only change is that `fixformula` can no longer be
  abbreviated as far as `f` or `fi` (`fix` and longer still work).

* **New `bace_validate()`: check how well BACE imputes YOUR data.** The
  natural call is on a finished run — `bace_validate(fit)` reuses the run's
  data, formulas, tree(s) and settings (now stored on the result as
  `$settings` and `$trees_used`) and re-runs the pipeline on freshly masked
  observed cells, scoring the imputations against the hidden truth; the
  direct data-and-formulas form works without a prior run. zRMSE (with its
  mean-filling baseline), 95% interval coverage (with a simulated fair
  expectation for counts) and accuracy (with fractional credit for tied
  votes) match the manuscript's definitions exactly; balanced accuracy and
  the Brier score are reported in addition. Repetitions (`reps`) give an
  honest between-repetition uncertainty. Inspired by Amelia's `overimpute()`
  (Honaker, King & Blackwell 2011, *J Stat Softw* 45(7)), but deliberately
  harsher (a share of cells is hidden at once). The docs state plainly what
  the check cannot do: it validates prediction where data were observed and
  cannot detect missing-not-at-random bias.

* **New `get_imputation_intervals()`: calibrated per-cell intervals for
  imputed values.** The standard multiple-imputation t-interval (Rubin 1987),
  which holds its 95% level for any `n_final` of about 20 or more — unlike
  the raw 2.5/97.5% quantiles of the imputed values (≈92% at `n_final = 50`),
  which the documentation previously suggested. The `n_final` documentation
  now points here.

* **Phylogenetic (tree) uncertainty can now be propagated through the
  imputations.** `bace()`, `bace_imp()` and `bace_final_imp()` accept a
  `multiPhylo` (or list of trees) as `phylo`: the convergence chain conditions
  on the first tree (its stabilisation diagnostic needs a fixed target), and
  each of the `n_final` final imputation runs is fit under its own tree —
  tree `m` for run `m`, selected deterministically from the run index (an
  evenly spaced subsequence when more trees than runs are supplied, balanced
  recycling otherwise), so `seed` reproducibility and serial/parallel equality
  are untouched. Between-imputation variance then includes tree uncertainty:
  proper multiple imputation with the phylogeny as a nuisance parameter drawn
  from its posterior (Rubin 1987; Nakagawa & de Villemereuil 2019, *Systematic
  Biology* 68:632–641). `bace_final_imp()` records `$tree_index` and
  `$tree_hashes`; `pool_posteriors()` detects differing trees, aligns
  species-level BLUPs by name across runs and drops internal-node BLUPs (node
  identities differ between trees), with `BACE_pooling$trees_vary` set;
  `with_imputations(tree = <multiPhylo>)` passes tree `m` to the downstream
  fit on imputation `m`. Runtime cost is negligible (`inverseA` is ~0.05 s at
  10,000 tips; the MCMC fits are indifferent to which tree they receive).

* **Trees are validated up front, with the offending tree named.** Tip
  coverage of the data, negative or zero-length branch lengths (with the
  standard `1e-8` remedy in the message), and every condition
  `MCMCglmm::inverseA()` enforces — most importantly ultrametricity, which
  previously surfaced as an opaque error deep inside a fit (and for a
  `multiPhylo`, only after hours) — are now checked at entry. Node labels
  (e.g. support values) are stripped internally: `inverseA()` rejects
  duplicated ones, and they fragment internal bookkeeping. A warning flags
  trees whose root depths suggest mixed branch-length scalings.

* **`inverseA()` is computed once per distinct tree per phase** instead of
  once per model fit (it was recomputed for every variable × sweep × run),
  via an internal per-tree cache cleared at each phase entry. Results are
  bitwise unchanged; the saving is seconds on long runs, but the cache is
  what makes per-run trees cheap.

# BACE 0.1.0

First public release.

### Fixes before first release (2026-10)

Found while verifying a co-author's review of the pre-release builds (four
reported issues, all confirmed, plus one found during verification); earlier
development builds are affected. None of these crash - they silently degrade
results, which is why they are listed prominently.

* **Warm-up sweeps are robust to transient fit failures.** With several
  variables missing, a warm-up sweep can resample a rare factor level out of
  existence, making the next conditional model's equations singular (most
  likely at strong phylogenetic signal). Such a failure now skips that
  variable's update for that sweep (the previous sweep's values carry
  forward - a skipped Gibbs update) with a warning, instead of aborting the
  run; a failure in the final, kept sweep still errors, so returned results
  are never silently degraded.

* **The default categorical prior no longer warns on every model build.**
  With the default `BACE.gelman = 2`, every categorical prior construction
  emitted "Not enough complete cases to compute Gelman prior or gelman=2"
  even though the pseudo-Gelman prior *is* the documented default - 60+
  spurious warnings per `bace()` call on a dataset with one categorical
  variable, enough to bury real warnings behind R's 50-warning cap. The
  warning now fires only on a genuine unrequested fallback (`gelman = 1`
  with fewer complete cases than J). `bace_final_imp()` also no longer
  builds (and then discards) a categorical prior when
  `ovr_categorical = TRUE` delegates prior construction to the one-vs-rest
  binary sub-models - matching `bace_imp()`, which never built one there.
  Prior values, imputations and RNG streams are unchanged.

* **Categorical imputation silently ignored the phylogeny.** The default
  (multinomial) prediction path looked species effects up under the wrong
  MCMCglmm column name, never matched, and used 0 - predictions equalled a
  non-phylogenetic model's. Fixed (exact prefix matching, with a tripwire
  warning should naming ever drift again); with `species = TRUE`, the second
  species effect is now also included for categorical and ordinal predictions.
* **Final imputations understated uncertainty when two or more variables were
  missing.** Each final run filled the variables once, in formula order, from
  the same point-imputed start, so the first-listed variable conditioned on
  identical point estimates of the others in every run (between-imputation
  variance collapsed; formula order mattered). The final phase now performs
  `sweeps` stochastic passes per run (default 3; warm-up passes use a short
  post-burn-in tail, so the default costs about 1.4x a single pass; collapses
  to 1 pass when at most one variable is missing, where a single pass is
  exact). New optional `seed` gives per-run RNG streams: results reproduce
  and are identical whether runs execute serially or in parallel; the
  parallel fallback now redoes only failed runs.
* **Reported category probabilities were overconfident.** `prob_preds` (and
  the chained phase's argmax fills) evaluated probabilities at posterior-mean
  coefficients; they now average the per-draw probabilities over the
  posterior (imputed values themselves were always proper draws and are
  unchanged in distribution). `prob_preds` is now documented.
* **MCMC diagnostics now cover ordinal cutpoints and the final-phase
  models.** Previously only chained-phase fixed effects and variances were
  checked; the models actually pooled for inference were never diagnosed.
  `bace_final_imp()` results carry a `diagnostics` summary (worst-case ESS /
  Geweke across runs, including one-vs-rest sub-models); printing shows
  flagged parameters only.
* **Priors are now user-settable.** New `priors` argument (named list, one
  element per variable) on `bace()`, `bace_imp()` and `bace_final_imp()`;
  blocks `G`/`R`/`B` replace the defaults wholesale, with validation before
  any MCMC runs. Overriding `R` is refused for threshold/ordinal/categorical
  variables (fixed latent scale). Defaults are unchanged. `bace()` now warns
  about unknown arguments instead of silently ignoring them, and both phases
  read the same `BACE.gelman` option (the final phase used to hard-code a
  different categorical fixed-effect prior than the chained phase).

* Core imputation pipeline: `bace()` orchestrates chained-equations imputation
  (`bace_imp()`), convergence assessment of the chained loop
  (`assess_convergence()`, with a family of `plot_*_convergence()`
  visualisations), multiple posterior-predictive final imputations
  (`bace_final_imp()`, parallelisable via `n_cores`), and posterior pooling
  (`pool_posteriors()`), with automatic retry of the chained phase on
  non-convergence.
* Every per-variable imputation model is a phylogenetic mixed model fitted
  with MCMCglmm; supported response types are gaussian, poisson (count),
  binary and ordinal (threshold scale), and unordered categorical
  (multinomial probit by default; one-vs-rest binary threshold models
  available via `ovr_categorical = TRUE`).
* Deterministic variable-type detection from column storage: doubles are
  gaussian, non-negative integers are poisson counts, factors map to
  threshold/multinomial models by level count and orderedness.
* Rubin's-rules pathway for arbitrary downstream analysis models:
  `with_imputations()` fits any model to each imputed dataset (with optional
  phylogeny pass-through) and `pool_mi()` combines estimates with Rubin's
  rules (Barnard-Rubin small-sample degrees of freedom, `fmi`/`riv`
  diagnostics; matches `mice::pool` to machine precision on identical inputs).
* Accessors for fitted objects: `get_pooled_model()`, `get_imputed_data()`.
* Pre-flight phylogenetic signal screening with `phylo_signal_summary()`
  (latent-scale H2 for all variable types, plus lambda, K and D where
  applicable, with MCMC reliability flags).
* Simulation engine `sim_bace()` for phylogenetically structured mixed-type
  datasets with known ground truth, missingness mechanisms, and optional
  within-species replication (`sim_tree()` simulates the phylogeny alone).
* The public API is deliberately small: imputation pipeline, pooling,
  accessors, diagnostics/plots, signal screening, and the two simulation
  entry points. All other helpers are internal and unexported.
