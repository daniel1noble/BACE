# Plan: four post-submission features (2026-11)

Dan picked these four (2026-11-07) from the "what would make BACE better" discussion:

1. **Tree-per-imputation** — accept a `multiPhylo` and use tree *m* for final run *m*,
   so between-imputation variance includes phylogenetic (tree) uncertainty.
2. **`bace_validate()`** — mask a fraction of *observed* cells, re-impute, score
   recovery on the user's own data (zRMSE, interval coverage, category accuracy).
3. **Checkpoint/resume** — long runs (22–28 h in our case studies) survive crashes;
   per-final-run caching with a settings fingerprint.
4. **Runtime** — attack where the time actually goes (final phase dominates), without
   touching `n_final` (it is load-bearing for per-cell interval calibration,
   `R/bace.R:24-33`).

Provenance: anchored to a full code map (2026-11-07) and then adversarially reviewed
by three independent passes (compatibility, statistics, tests) against the actual
code. The reviews confirmed every mechanical claim (per-run RNG streams, positional
pooling, the Jacobi convergence structure, dots handling, the retry ladder) and found
seven must-fix defects in the original feature specs, all folded in below — the most
consequential: multiPhylo today dies with an opaque `unused argument (2)` inside
`inverseA`; pooling's real failure mode is *silent tip-BLUP misalignment* (not a hard
error); my original zRMSE spec disagreed with the manuscript's definition by ~√2; and
the original feature-1 gate could not detect tree-1 initialization leakage.

## Ground rules that shaped every design decision

These came out of the code map and reviews; they are *constraints*, not preferences:

- **New arguments must be explicit.** `bace()` warns on any `...` and the workers
  silently drop dots (`R/bace.R:166-178`, `R/bace_imp.R:83`). Every new argument is
  added to the signature of `bace()` *and* the relevant worker, plus roxygen.
- **Features must land at the `bace_imp()`/`bace_final_imp()` level, with `bace()`
  forwarding.** Both case studies bypass `bace()` entirely and call the two phases
  directly (`vignettes/amphibian_ctmax/00_setup.R:33-48`,
  `vignettes/bird_social/02_run_bace.R:44-57`).
- **`bace()` calls `bace_final_imp` at TWO sites** — converged (`R/bace.R:299-317`)
  and non-converged/skip_conv (`R/bace.R:350-368`). Every new forwarded argument is
  plumbed through *both*, and each feature adds a forwarding test that runs with
  `skip_conv = TRUE` (the configuration the examples and case studies use).
- **`.model_fit`'s signature is frozen** at its current 9 arguments
  (`R/model_functions.R:16`): three existing test files mock it with the full
  signature and would break on any change (`test-final-sweeps.R:52-53, 210-211,
  242-243`). Per-run trees arrive through the *existing* `tree` argument.
- **Default paths stay byte-identical.** The bird vignette pins results to commit
  `dc79be1` with recorded seeds; the serial==parallel and fallback-stream tests
  (`tests/testthat/test-final-sweeps.R:131-191`) define the RNG contract. No new RNG
  consumption on any default path; no default value changes.
- **The final-phase RNG design is our friend.** With `seed` set, run *m*'s entire
  randomness comes from a stream computed from `(seed, m)` alone
  (`R/bace_final_imp.R:180-207, 221-224`); runs re-seed at entry, so any subset can
  execute in any order with bitwise-identical results (the failed-worker retry
  already relies on this, `R/bace_final_imp.R:52-56, 448-461`). Resume (3) and
  per-run trees (1) both build on it.
- **The convergence phase is unseeded by design** (`seed` reaches only the final
  phase, `R/bace.R:80-82`) — it is persisted as an object, never replayed.

## Speed is a first-class constraint (Dan, 2026-11-07)

Accuracy and calibration come first, but BACE's runtime (22–28 h in the case
studies) means **no feature may slow the default path, and every PR reports
wall-time** on a fixed reference fixture (hook: the existing `benchmark.yml` /
`sim-reference-benchmark.yml` GHA workflows and `dev/benchmark_engine.R`). Where the
time actually goes, measured: the final phase dominates (bird: 20 × 55–75 min final
runs vs 2.5 h chained), and within it, full MCMC fits — `n_final × sweeps ×
n_missing_vars` of them. The levers that move hours are 4b/4d (chain sizing), 4c
(parallelising the chained phase), and the already-shipped `n_cores`.

**The multiPhylo speed hazard was measured and is negligible** (2026-11-07, this
machine; independently reproduced by the review): `MCMCglmm::inverseA(tree,
nodes = "ALL")` costs 0.02 s at 100 tips, 0.05 s at 10,000 tips (~2 MB per Ainv).
Fifty posterior trees add ~2.5 s of matrix work and ~110 MB (computed lazily per
worker) to a multi-hour run; the MCMCglmm fits themselves are indifferent to which
tree they receive. Tree-per-imputation buys phylogenetic-uncertainty propagation at
essentially zero runtime cost. The same measurement honestly demotes 4a: caching
inverseA saves ~30 s on a case-study run — hygiene plus the keying mechanism
feature 1 needs, not a performance claim.

**The one unmeasured speed cost is feature 3's checkpoint WRITE time** (the plan
documents file *size* but gzip `saveRDS` on ~GB objects can cost minutes per run —
possibly hours across 50 runs). Feature 3's PR benchmarks the write path and exposes
a compression flag (`compress = FALSE` is usually the right trade on a 22-h run).

## Order of work

```
4a  inverseA cache            (tiny, standalone, zero statistical change)
1   multiPhylo                (touches the final-run loop)
3   checkpoint/resume         (same loop; after 1 so the fingerprint includes
                               trees from day one)
4b  final-phase MCMC args     (plumbing; what the case studies already hand-roll)
2   bace_validate()           (a consumer of everything above)
4c  convergence-phase         (Jacobi confirmed; needs its redesigned step-zero
    parallelism                test first; separate PR)
4d  ESS floor per run         (reframed by review; benchmark-gated; separate PR,
                               may not land)
```

Each lands as its own PR with its own tests; `devtools::test()` green and
`R CMD check --as-cran` clean are gates for every one. NEWS.md bullet per feature
(house style: bold one-sentence headline, plain-English mechanism, `NEWS.md:1-32`).

---

## Feature 1 — Tree-per-imputation (`multiPhylo`)

**What.** `phylo` accepts a `multiPhylo` (ape class; ape already in Imports). The
convergence phase uses `trees[[1]]`; final run *m* uses a tree chosen
**deterministically from the run index, never drawn**, so no RNG is consumed and the
stream contract survives:

- `n_trees ≥ n_final`: an evenly spaced subsequence
  `trees[round(seq(1, n_trees, length.out = n_final))]` (posterior tree files are
  usually chain-ordered or MCC-first, so taking the first 50 of 1,000 would not be
  exchangeable); a message names which trees were used.
- `n_trees < n_final`: recycle `((m-1) %% n_trees) + 1`. `K ≥ n_final` is NOT
  required — at `n_final = 50, K = 20` balanced reuse deflates the expected
  between-tree component of B by only ~3% (review's calculation), and the formally
  exact framework for reuse is nested/two-stage MI, which plain Rubin approximates
  adequately here. Documented, with the advice to supply a representative random
  sample of posterior trees.

Statistical rationale (goes in the roxygen and, later, the ms): this is proper MI
with the phylogeny treated as a nuisance parameter drawn from its (empirical)
posterior — Rubin (1987, ch. 3–4); tree-sample averaging precedent: Huelsenbeck,
Rannala & Masly 2000 *Science* 288:2349; Housworth & Martins 2001 *Syst Biol*
50:628; de Villemereuil et al. 2012 *BMC Evol Biol* 12:102; Rangel et al. 2015
*Evolution* 69:1301; software precedent *mulTree* (Guillerme & Healy). **The
convergence phase conditions on a single tree by necessity, not convenience**: the
stabilisation diagnostic reads the settling of deterministic point imputations
(`sample = FALSE` invariant, `R/bace_imp.R:331-342`), and cycling trees across
convergence runs would inject tree-to-tree jitter into exactly the trajectories the
diagnostic needs stable. Tree uncertainty is propagated where it belongs — the
between-imputation variance of the final draws.

**Code changes.**
- **Shared validator helper** (new, e.g. `.validate_trees()`), called at the top of
  `bace_imp()` AND `bace_final_imp()` (so the two-phase workflow is covered; `bace()`
  inherits it). Today a multiPhylo dies inside `inverseA` with
  `unused argument (2)` (reproduced by review: inverseA treats it as a pedigree).
  Checks per tree: class; tip labels match the data's cluster values; zero/negative
  branch lengths; **tip+node label uniqueness** (inverseA hard-stops on duplicated
  node labels — e.g. support values stored as node labels — which would otherwise
  surface hours into the final phase; strip or error at entry); and a
  **scale-consistency warning** when across-tree root-depth CV exceeds ~10–25%
  (σ²_phylo is defined relative to branch-length scale; a BEAST posterior's genuine
  root-height variation should pass, an accidentally rescaled tree should not).
- `R/bace.R` — accept phylo|multiPhylo; pass `trees[[1]]` to `bace_imp()` and the
  selected list to `bace_final_imp()` at **both** call sites; `phylo_signal = TRUE`
  preview runs on `trees[[1]]` with a message (`phylo_signal_summary` hard-requires
  a single phylo, `R/phylo_signal_summary.R:165-166`).
- `R/bace_imp.R` — accept multiPhylo by taking `trees[[1]]` (with the validator and
  a message), so the documented two-phase workflow works unchanged.
- `R/bace_final_imp.R` — thread the per-run tree through `run_one`/`.one_final_run`
  to the existing `tree`/`phylo` arguments of `.model_fit` (:350) and
  `.fit_predict_ovr` (:301). Record `$tree_index` (integer vector, length
  `n_final`) and a hash of each tree on the `bace_final` output — today no result
  object stores the tree at all.
- `R/pool_posteriors.R` — **the one real hazard, corrected by review.** Stacking
  takes dims/colnames from run 1 and inserts later runs BY POSITION
  (`R/pool_posteriors.R:165-235`), and `Sol` carries tip + internal-node BLUP
  columns (`pr = TRUE`, `inverseA(nodes = "ALL")`). The review *reproduced* the real
  failure mode: same-tip trees give **identical column-name sets in different
  orders**, so nothing ever hard-errors — tip BLUPs land under the wrong tip names
  and node columns silently refer to different clades. Fix: when `$tree_index` says
  the runs' trees differ (detection must come from metadata — column names cannot
  reveal it), stack fixed effects, VCV, CP and Liab as now (meaning is
  tree-invariant), align tip-BLUP columns **by name**, and drop node columns
  identified as the **complement of the tip-label set** (never a `Node\d+` regex —
  labelled trees name node columns after the labels). Set
  `BACE_pooling$trees_vary = TRUE`; document that the pooled skeleton's `$X`/`$Z`
  (copied wholesale from run 1, `:267`) describe tree 1 only.
- `R/with_imputations.R` — `tree =` already anticipates this (docstring "leaves room
  for a future posterior-tree workflow"): accept a list of trees, inject tree *m*
  into fit *m*, validate list length against the number of datasets.

**Does it break anything?** No existing test constrains `class(phylo)`; "multiPhylo"
appears nowhere in R/, tests/, or vignettes. Single-tree calls take the
`n_trees == 1` path, byte-identical. Tree choice consumes no RNG, so the
reproducibility contract and published case-study seeds survive. The OVR path flows
automatically (its J sub-fits go through `.model_fit`). The existing tip-mismatch
test (`test-validation.R:29-45`) is `expect_error` with no regexp, so earlier,
friendlier messages cannot break it.

**Workflow impact.** None for single trees. For a `multiPhylo`: swap one argument,
in either the one-call or the two-phase workflow. Pooled-BLUP caveat in roxygen and
vignette.

**Tests.**
- CRAN-safe (no MCMC): class acceptance; per-tree tip-mismatch errors naming the
  offending tree, asserted against BOTH entry points (`bace_imp`, `bace_final_imp`);
  node-label-uniqueness and scale-consistency warnings; subsequence/recycling
  arithmetic via an exported-internal helper; `with_imputations(tree = <list>)` with
  plain `lm` fits (list length validation, recycling, `.f` without a tree formal).
- Tiny-MCMC, `skip_on_cran`: `$tree_index` recorded; single-tree results unchanged;
  the `phylo_signal` preview message fires.
- **The identity test (review's redesign — the original "two tree groups differ"
  was vacuous** since per-run streams make any two runs differ): with seed *s* and
  trees (A, B), run multiPhylo `n_final = 4` plus single-tree-A and single-tree-B
  runs under the same seed; assert multi runs 1,3 ≡ single-A runs 1,3 and multi
  runs 2,4 ≡ single-B runs 2,4, and single-A ≢ single-B. This catches both
  "tree silently ignored" and stream contamination at no extra fit cost.
- Serial == parallel under multiPhylo (`skip_on_os("windows")` — GHA runs
  windows-latest, so the omission is a red CI, not a style point).
- Pooling: CRAN-safe via a handcrafted fake `bace_final` whose per-run `Sol`
  matrices have identical column-name sets in permuted order — assert by-name tip
  alignment, node-column dropping, `trees_vary` driven by `$tree_index`.
- Statistical gate (dev/ script; redesigned by review because the original
  "variance must not shrink" gate was blind to **tree-1 initialization leakage** —
  every final run starts from the tree-1 converged state, and residual dependence
  after 3 sweeps pulls run-m draws toward the tree-1 predictive centre, shrinking
  the between-tree variance component in exactly the anti-conservative direction):
  on a clade-uncertainty sim, (a) a gold-standard arm runs the *chained phase
  itself* on tree m (cheap at sim scale) — the multiPhylo shortcut's B and
  per-tree-group mean separation must match it within ~10%/MC error; (b) a sweeps
  dose-response (sweeps ∈ {1, 2, 3, 5, 10}, paired seeds) — tree-group separation
  and B must plateau by sweeps = 3, or the default is raised for multiPhylo runs.
  Pre-registered thresholds in the script header, template
  `dev/15_sweeps_calibration.R`.

---

## Feature 2 — `bace_validate()` (overimputation diagnostics)

**What.** One exported call that answers "does the imputation model predict well
where data were observed": mask a fraction of **observed** cells, run the full
pipeline on the masked data, score recovery against the held-out truth. Replaces
the manual recipe the vignette already teaches (`vignettes/bace.qmd:708-789`) and
both case studies hand-rolled.

**The honesty caveat goes in roxygen AND the print method** (review must-fix):
MCAR masking of observed cells validates predictive accuracy *on the observed-data
region*. It is necessary, not sufficient, for trusting imputations of the
actually-missing cells, and it cannot detect MNAR bias (in comparative data the
genuinely missing species are often small/rare — an extrapolation region). Design
precedent and citation: Amelia's `overimpute()` (Honaker, King & Blackwell 2011,
*J Stat Softw* 45(7)). MAR-style stress masking (mask probability depending on
observed covariates; pattern already in `dev/15_sweeps_calibration.R:14-16`) is the
documented v2 extension.

**Design decisions (all review-corrected).**
- **A fresh run on masked data, never mutation of an existing object** (pinned by
  `test-performance-fixes.R:28-77`). Masking changes the input data upstream of
  `missing_matrix` (`R/bace_imp.R:177-181`), which picks the masked cells up
  automatically.
- Masking: MCAR v1, `prop` per modelled variable, stratified for factors so no
  observed level empties (guard + **warning** — the dev masker silently `next`s on
  <10 observed cells, `dev/benchmark_engine.R:153`; the port adds the warning).
  Never the cluster column; only variables with a formula. Per-rep mask seeds are
  derived deterministically from `(seed, rep)` **outside** the pipeline call, so
  they can never interleave with the final-phase stream derivation.
- **zRMSE — the manuscript's definition, exactly** (review must-fix; my original
  spec was wrong twice): numerator = RMSE of the **across-`n_final` mean** at
  masked cells; denominator = SD of the variable's **remaining observed cells**
  (the ms's `s_obs`, `ms/ms.qmd:308-313`; also what the rubin-study scorer does,
  `dev/rubin_study_2026-09-26/code/rubin_cell.R:110-114`). Do **not** port
  `dev/benchmark_engine.R`'s `.metrics_continuous` as-is — it averages
  per-single-draw RMSEs, a different estimand that is ~√2 larger for a calibrated
  imputer and incomparable to the published figures.
- **Interval coverage — per-cell Rubin t-interval, matching the ms**
  (`ybar ± t_{M-1,0.975} · sqrt((1 + 1/M)·s²)`, `ms/ms.qmd:324`,
  `dev/rubin_study_2026-09-26/code/rubin_lib.R:80-95`), which is ≈nominal at any
  M ≥ ~20 — not empirical quantiles, whose ~92% ceiling at `n_final = 50`
  (`R/bace.R:24-33`) would be misread against the ms's ~95%. The print shows
  observed coverage next to its expectation and a binomial MC band given
  `n_masked`; zero-width intervals (degenerate low counts) are flagged.
- **Poisson noise floor** (review): a count trait with weak between-species signal
  shows zRMSE near 1 under a *perfect* imputer (irreducible Poisson noise — the ms
  needed oracle floor lines for exactly this). Print a plug-in floor
  (`sqrt(mean imputed rate)/s_obs`) next to the poisson zRMSE, with roxygen text
  "zRMSE near the floor means the imputer extracts all available signal".
- Categorical/threshold: accuracy and balanced accuracy vs the modal-class
  baseline; **Brier from the MI-averaged `prob_preds`** (the documented pooling,
  `R/bace.R:113-124`) — not from imputed-class indicator proportions (the dev
  helper's different definition; pick one, this one).
- **A named scorer seam**: internal
  `.score_validation(truth, masked_cells, imputed_datasets, prob_preds, types)` —
  so metrics are testable against a handcrafted fake without MCMC. Types come from
  the pipeline's own `$types`, never re-detected (the signal-path detector
  classifies double-stored whole numbers differently,
  `R/phylo_signal_summary.R:288-321`).
- Output: class `bace_validation` — `$table` (variable × metric × baseline ×
  expectation/floor), `$masked_cells`, `$settings`, `$n_final`; print with
  interpretation bands; `reps` (default 1) with per-rep and aggregated tables and
  per-metric MC error.

**Signature.** Mirrors `bace()` explicitly (dots are a dead end), *including*
`skip_conv` and `max_attempts` — the auto-retry ladder multiplies per-rep cost (up
to `max_attempts` chained phases), and the cost note says so.

**Does it break anything?** No — a new export consuming the public pipeline.

**Workflow impact.** Additive; the vignette's benchmarking section becomes one
call; recommended "validate, then impute".

**Tests.**
- CRAN-safe: masked cells ⊆ observed; cluster untouched; level-preservation and
  <10-cells warnings; mask reproducibility (same seed ⇒ same mask; reps differ);
  `.score_validation` arithmetic against a handcrafted fake `bace_complete`
  (feasible: `get_imputed_data` needs only class + `$imputed_datasets`,
  `R/accessors.R:108-124`; `prob_preds` are plain matrices) with hand-computed
  answers — including **zRMSE ≈ 1 for a mean-filling imputer and 0 for a perfect
  one**, and a replicated-species case (make.unique rownames make name-joins
  non-trivial, `R/bace_final_imp.R:383-389`).
- Tiny-MCMC end-to-end (`skip_on_cran`): structure-only assertions (class, table
  columns, finite metrics) on the small fixture. Metric-level gates live in a dev/
  script on a ≥60-species fixture with `prop ≥ 0.3` — the conditions under which
  the borrowed bounds (`test-imputation_accuracy.R:36-47`) are stable; an 8-tip
  fixture at `prop = 0.1` masks 1–3 cells and makes r/coverage degenerate.

---

## Feature 3 — Checkpoint/resume

**What.** `checkpoint_dir` on `bace()`, `bace_imp()` and `bace_final_imp()`. Final
phase: one rds per completed run. Convergence phase: one rds for the finished
`bace` object. Resume = same arguments + same directory; completed work loads, only
missing runs execute. Replaces the hand-built scaffolding in both case studies.

**Why the design is safe.** The per-run payload already exists
(`list(models, dataset, prob_preds, ovr_models)`, `R/bace_final_imp.R:424-425`);
post-loop assembly is deterministic and positional (`:469-498`) so cached runs
splice by index. With `seed`, resume is bitwise-identical to an uninterrupted run
(stream-per-run property). **`checkpoint_dir` therefore requires `seed`** — hard
error otherwise, message explaining why. A cached run is final — warm-up sweep
failures legitimately alter a run's internal RNG consumption
(`R/bace_final_imp.R:398-406`), so cached runs are never re-run-and-compared.
Migration note (review): the bird case study already passes `seed`; the amphibian
helper seeds externally and calls `bace_final_imp` *without* the seed argument
(`vignettes/amphibian_ctmax/00_setup.R:77-81`) — migrating it to checkpointing
changes its seeding semantics, and its archived draws are not bitwise-recoverable
under the stream contract. Say so in its README when migrating.

**Fingerprint (manifest.rds).** Stores the settings *themselves* so mismatches can
be named: data, tree(s) + `tree_index`, built formulas, post-transform
`ran_phylo_form`, nitt/thin/burnin lists AFTER `nitt_cat_mult` (and after any 4d
resize), `eff_sweeps`, priors, `ovr_categorical`, `species`, `seed`, `n_final`, the
`BACE.gelman` option, package version, **and R version as its own named field**.
Large objects are compared by hash — **with the serialization header stripped**
(`serialize(version = 2)` minus its fixed 14-byte header, then `tools::md5sum` on a
temp file): the review reproduced that raw `serialize()` output embeds the running
R version and native encoding, so after an R upgrade identical data would hash
differently and the refusal would name the *wrong* field ("data" instead of
"R version"). No new dependency (digest stays out; base `::` use is the precedent).
On mismatch: refuse, naming the differing fields; never silently mix.

**Code changes.**
- `R/bace_final_imp.R` — in `run_one` (`:429-432`): check cache → load or execute →
  write. **Fork-safe atomic write: `tempfile(tmpdir = checkpoint_dir)` then
  `file.rename()` with its return value checked** (review: `tempfile()` defaults to
  node-local `tempdir()`, and `file.rename` across filesystems — NFS scratch is the
  motivating case — returns FALSE rather than falling back to copy; unchecked, the
  checkpoint is silently lost). Document that rename atomicity on NFS is
  best-effort. Expose `checkpoint_compress` (default TRUE; FALSE recommended for
  big runs) and **benchmark the write path in this PR** — the one unmeasured speed
  cost.
- `R/bace_imp.R` / `R/bace.R` — persist the finished `bace` object (the phase is
  unseeded; persisted, never replayed). The retry ladder re-runs `bace_imp` with
  `runs × 1.5^attempt` (`R/bace.R:258-286`); each attempt's object overwrites the
  chain checkpoint with its `runs` and attempt number recorded, so a resumed
  `bace()` continues the ladder instead of restarting it. Plumb `checkpoint_dir`
  through **both** `bace_final_imp` call sites.
- `R/bace.R` — on resume with `final_results` incomplete, re-enter the final phase
  before pooling (pooling is called at `R/bace.R:328` and `:374`).

**Does it break anything?** Default `checkpoint_dir = NULL` is a no-op; no RNG or
behaviour change. Known documented cost: per-run files are heavy (models saved with
`pr/pl/saveX/saveZ = TRUE`, `R/model_functions.R:129-130`; ~GB per imputation on
18k rows in the amphibian study). A `checkpoint_slim` option is future work (slim
reload is proven viable, `vignettes/bird_social/02_run_bace.R:52-74`).

**Workflow impact.** One new argument (+ seed if not already set). Case-study
scripts' manual `saveRDS` chunking disappears.

**Tests.**
- **The resume oracle, redesigned** (review must-fix — the original "fail run 3 on
  its first call" is swallowed by the warm-up tolerance and can never crash a run;
  and call-counting mocks cannot identify "run 3" under mclapply): fixture with
  **one missing variable so `eff_sweeps` collapses to 1** (`:177-178`) and every
  fit is fatal, `n_cores = 1`; run A uninterrupted with seed; run B with
  `checkpoint_dir` and a mock that fails at the full-nitt fit of run 3
  (house pattern keys on nitt, `test-final-sweeps.R:241-248`), wrapped in
  `expect_error`; assert runs 1–2's cache files exist; resume (mock gone) →
  `all_datasets` identical to A.
- Parallel resume as a **real-fit** test (`skip_on_cran`, `skip_on_os("windows")`):
  resume under `n_cores = 2` equals serial resume. (Mock state is serial-only;
  bindings do reach forks but counters don't return.)
- CRAN-safe: `checkpoint_dir` without `seed` errors; fingerprint mismatch on
  handcrafted manifests (changed nitt / data / tree / **package version** /
  R version) refuses and names the field; atomic-write naming; cached files loaded
  not recomputed (spy counts fits, serial).
- **Retry-ladder resume** (tiny-MCMC, mock `assess_convergence` to force one
  retry): interrupt after attempt 2's chain checkpoint, resume, assert the ladder
  continues from attempt 2 rather than refusing on the changed `runs`.

---

## Feature 4 — Runtime

Where the time goes, measured: fits per final run = `eff_sweeps × n_missing_vars`
(+J per OVR categorical) MCMCglmm fits; bird study: each final run 55–75 min at
nitt = 500k vs a 2.5 h chained phase.

**4a. Cache `inverseA` per tree (hygiene, not a performance claim).**
`.model_fit` recomputes `inverseA(tree, nodes = "ALL")` on every fit
(`R/model_functions.R:19`). Measured: 0.05 s at 10k tips → ~30 s saved per
case-study run; the real value is the per-tree keying feature 1 needs, plus removing
the same redundancy in `phylo_signal_summary` (`:229, 415-427`). Structure (review):
**`.get_ainv(tree)` delegating to an inner `.compute_ainv(tree)`** — the
once-per-tree assertion spies on the inner function (the outer is called per fit by
design); lazy per-worker computation; `.model_fit`'s signature unchanged. Key must
distinguish trees with identical tip labels but different topologies (hash of the
tree's edge matrix + lengths + labels, not object identity — R has no stable
identity across forks).
Tests: inner-function spy (once per tree per phase); **collision test** — two trees,
same tips, different topologies → different Ainv retrieved per run; byte-identical
results with/without cache under the same seed.

**4b. Final-phase MCMC arguments.**
`bace()` currently reuses one `nitt/thin/burnin` for both phases; the case studies
hand-tune the final phase an order of magnitude longer (bird: 50k vs 500k). Add
`final_nitt/final_thin/final_burnin = NULL` to `bace()`. **Resolution lives in
`bace()`**: `final_nitt <- final_nitt %||% nitt` before the two `bace_final_imp`
call sites — `.standardize_mcmc_params(NULL)` errors loudly (`R/bace_imp.R:24-29`),
so a forwarded NULL cannot slip through silently, but the substitution must still
be in `bace()`, not the worker.
Breaks nothing: NULL defaults inherit byte-identically (verified pathway:
both workers standardize and apply `nitt_cat_mult` internally and identically).
Tests: forwarding spy incl. the `skip_conv = TRUE` branch; list-valued `final_nitt`
through `bace()`; `final_nitt` × `nitt_cat_mult` interaction (the multiplier must
still apply, `R/bace_final_imp.R:144-153`); defaults unchanged (the genuine pin is
`test-posterior_sampling.R:355-363`; note `test-bace_wrapper.R:196-203` self-adjusts
to chain length and guards nothing).

**4c. Parallelise the convergence phase across variables.**
**The Jacobi claim is now confirmed by two independent reviews**: within a run every
variable is fit against the previous run's dataset with only its own response
re-NA'd (`R/bace_imp.R:233, 244-248`); updates land in a write-only copy (`:356-358`)
published at run end (`:366`). The convergence diagnostic consumes only the dataset
sequence (`R/convergence_diagnostics.R:59`), whose Markov kernel is order-invariant
given Jacobi, and the phase is unseeded, so parallelising changes no distributional
property and violates no reproducibility contract. The final phase is intentionally
Gauss-Seidel (`R/bace_final_imp.R:374-376`) — the asymmetry is correct.
**Step-zero verification test, redesigned** (review: the naive permutation test
fails spuriously because `.data_prep` consumes global RNG — empirical draws for
factor warm-starts at r == 2, `R/prep_functions.R:211-214` — in visit order): use an
all-gaussian-predictor fixture (warm-start is then deterministic means,
`:203-209`) with mocked deterministic fits, and assert `pred_missing_run` is
byte-identical under permuted formula order. New `n_cores_conv = 1` default; the
PR's wall-time report states the realistic ceiling (Σtᵢ / max tᵢ over per-variable
fit times). Separate PR.

**4d. Per-run ESS floor (reframed by review; benchmark-gated; may not land).**
My original target — pooled ESS across stacked chains — was the wrong quantity: the
MC error of a pooled functional = within-run chain noise (per-run ESS) + B/M
between-imputation noise, and B/M is reducible only via `n_final`. Once per-run ESS
is adequate, extra nitt buys nothing; sizing *down* to a pooled target would drive
per-run ESS to target/M (≈20 at target 1000, M 50) — territory the package's own
diagnostics flag as unhealthy (ESS ≤ 200, `R/model_functions.R:1638-1640, 1731`) and
where per-run `prob_preds` degrade. Reframe: **opt-in `target_ess_per_run`
(default NULL; floor 200, matching `.check_final_diagnostics`; cf. Vehtari et al.
2021, *Bayesian Analysis* 16:667)** — estimate per-variable ESS/iteration from the
convergence phase's last-run models, resize `nitt_list` up or down to the floor,
report pooled draws = M × per-run ESS as information. Explicit `nitt` always
honoured. Resolved lists go into the checkpoint manifest.
Gate (dev/ study, `dev/15_sweeps_calibration.R` template): pooled estimates and CI
widths within MC error of the fixed-length reference; per-cell interval coverage
unchanged; **zero new diagnostic flags**; categorical `prob_preds` stability; wall
time reported. If the gate fails, drop the feature and keep the study as the
documented reason.

**Not a lever: `n_final`** (empirical per-cell intervals need it, `R/bace.R:24-33`).

---

## Cross-feature interactions

- **1 + 3:** the fingerprint includes the tree list and `tree_index`; each cached
  run records its tree index — a cache can never silently pair a run with a
  different tree.
- **1 + 4a:** the Ainv cache is keyed per tree hash; multiPhylo gets one entry per
  distinct tree, computed lazily per worker.
- **1 + pooling consumers:** `pool_posteriors` consumes only `$all_models`;
  `with_imputations`/`pool_mi` consume only `$all_datasets` — the dataset-pooling
  route is tree-agnostic. OVR categoricals stay invisible to `pool_posteriors` by
  design (`R/bace_final_imp.R:312-318`).
- **4b/4d + 3:** chain lengths enter the fingerprint after resolution, so resumes
  compare what actually ran.
- **2 + everything:** `bace_validate()` is a pure consumer and inherits multiPhylo,
  checkpointing and the final-phase arguments for free.

## What this does NOT change

- The two-phase contract: convergence predicts with `sample = FALSE`
  (`R/bace_imp.R:305, 331-342`), final draws with `sample = TRUE`.
- Default MCMC settings, priors, `n_final = 50`, `sweeps = 3`, `.model_fit`'s
  signature.
- Object shapes — only additive fields (`$tree_index`, checkpoint metadata).
- CRAN v0.1.0 submission content (on hold awaiting Dan's go). All four features are
  0.2.0 material unless Dan decides otherwise.

## Effort estimate (coding + tests, excluding gate studies)

| Feature | Size | Risk |
|---|---|---|
| 4a inverseA cache | ~½ day | low |
| 1 multiPhylo | ~3–4 days (pooling alignment + gate study are most of it) | medium |
| 3 checkpoint/resume | ~2–3 days (fingerprint + retry-ladder semantics) | medium |
| 4b final_* args | ~½ day | low |
| 2 bace_validate | ~2–3 days (correct metrics + seam + vignette rewrite) | low |
| 4c conv parallelism | ~1–2 days after the redesigned step-zero test | medium |
| 4d per-run ESS floor | ~2 days + calibration study | high (gated) |
