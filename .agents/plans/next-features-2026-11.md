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

Everything below is anchored to the current code (mapped 2026-11-07; all `file:line`
refer to that state). The plan was adversarially reviewed before being adopted.

## Ground rules that shaped every design decision

These came out of the code map and are *constraints*, not preferences:

- **New arguments must be explicit.** `bace()` warns on any `...` and the workers
  silently drop dots (`R/bace.R:166-178`, `R/bace_imp.R:83`). Every new argument is
  added to the signature of `bace()` *and* the relevant worker, plus roxygen.
- **Features must land at the `bace_imp()`/`bace_final_imp()` level, with `bace()`
  forwarding.** Both case studies bypass `bace()` entirely (memory/chunking) and call
  the two phases directly (`vignettes/amphibian_ctmax/00_setup.R:33-42`,
  `vignettes/bird_social/_bird_social.qmd:95-117`). A `bace()`-only feature misses
  the documented big-data workflow.
- **Default paths stay byte-identical.** The bird vignette pins results to commit
  `dc79be1` with recorded seeds; the serial==parallel and fallback-stream tests
  (`tests/testthat/test-final-sweeps.R:131-191`) define the RNG contract. Therefore:
  no new RNG consumption on any default path, and no default value changes.
- **The final-phase RNG design is our friend.** With `seed` set, run *m*'s entire
  randomness comes from a stream computed from `(seed, m)` alone
  (`R/bace_final_imp.R:180-207, 221-224`); runs re-seed themselves at entry, so any
  subset can be executed in any order with bitwise-identical results. The existing
  failed-worker retry already relies on this (`R/bace_final_imp.R:52-56, 448-461`).
  Resume (feature 3) and per-run trees (feature 1) both build on it.
- **The convergence phase is unseeded by design** (`seed` goes only to the final
  phase, `R/bace.R:80-82`), so it can never be re-derived from settings — it must be
  persisted as an object, not replayed.

## Speed is a first-class constraint (Dan, 2026-11-07)

Accuracy and calibration come first, but BACE's runtime (22–28 h in the case
studies) means **no feature may slow the default path, and every PR reports
wall-time** on a fixed reference fixture (hook: the existing `benchmark.yml` /
`sim-reference-benchmark.yml` GHA workflows and `dev/benchmark_engine.R`). Where
the time actually goes, measured: the final phase dominates (bird: 20 × 55–75 min
final runs vs 2.5 h chained), and within it, full MCMC fits — `n_final × sweeps ×
n_missing_vars` of them. So the levers that move hours are 4b/4d (chain sizing),
4c (parallelising the chained phase), and the already-shipped `n_cores`.

**The multiPhylo speed hazard was measured and is negligible** (2026-11-07,
this machine): `MCMCglmm::inverseA(tree, nodes = "ALL")` costs 0.02 s at 100 tips,
0.05 s at 10,000 tips (2.2 MB per Ainv). Fifty posterior trees therefore add ~2.5 s
of matrix work and ~110 MB (computed lazily per worker, only for the trees that
worker's runs use) to a multi-hour run; the MCMCglmm fits themselves are
indifferent to which tree they receive (same dimensions, same sparsity class).
Tree-per-imputation buys phylogenetic-uncertainty propagation at essentially zero
runtime cost. The same measurement honestly demotes 4a: caching inverseA saves
~0.05 s × ~600 fits ≈ 30 s on a case-study run — it stays in the plan as hygiene
and as the keying mechanism feature 1 needs, not as a performance claim.

## Order of work

```
4a  inverseA cache            (tiny, standalone, zero statistical change)
1   multiPhylo                (touches the final-run loop)
3   checkpoint/resume         (touches the same loop; after 1 so the fingerprint
                               includes trees from day one)
4b  final-phase MCMC args     (plumbing; what the case studies already hand-roll)
2   bace_validate()           (a consumer of everything above)
4c  convergence-phase         (needs its own verification first; separate PR)
    parallelism
4d  ESS-targeted auto-sizing  (benchmark-gated; separate PR, may not land)
```

Each lands as its own PR with its own tests; `devtools::test()` green and
`R CMD check --as-cran` clean are gates for every one. NEWS.md bullet per feature
(house style: bold one-sentence headline, plain-English mechanism, `NEWS.md:1-32`).

---

## Feature 1 — Tree-per-imputation (`multiPhylo`)

**What.** `phylo` accepts a `multiPhylo` (ape class; ape already in Imports,
`DESCRIPTION:29-44`). The convergence phase uses `trees[[1]]`; final run *m* uses
tree `((m-1) %% n_trees) + 1` — **deterministic from the run index, never drawn**,
so no RNG is consumed and the stream contract survives. Between-imputation variance
then automatically includes tree uncertainty (standard MI reasoning; one tree per
imputation is the Rubin-style approach to tree uncertainty).

**Code changes.**
- `R/bace.R` — accept `phylo` of class `phylo` or `multiPhylo`; validate every tree
  once at entry (class, tip labels match the data's cluster column, zero/negative
  branch lengths — today there is NO tree validation anywhere,
  `R/bace_imp.R:98-129`); pass `trees[[1]]` to `bace_imp()` and the whole list to
  `bace_final_imp()`. `phylo_signal = TRUE` preview runs on `trees[[1]]` with a
  message (`phylo_signal_summary` hard-requires a single `phylo`,
  `R/phylo_signal_summary.R:165-166`).
- `R/bace_final_imp.R` — thread a per-run tree: `.one_final_run(run, ...)` closes
  over a single `phylo` today and hands it to `.model_fit`/`.fit_predict_ovr` at
  lines 301 and 350; replace with `trees[[tree_index(run)]]`. Record
  `$tree_index` (integer vector, length `n_final`) and the tree list (or its hash)
  on the `bace_final` output — today *no result object stores the tree at all*
  (`R/bace_imp.R:374-382`, `R/bace_final_imp.R:475-499`, `R/bace.R:385-397`), which
  also matters for feature 3's fingerprint.
- `R/pool_posteriors.R` — **the one real hazard.** Stacking takes dimensions and
  column names from the FIRST run's model and inserts later runs BY POSITION
  (`R/pool_posteriors.R:165-235`); `Sol` contains tip AND internal-node BLUP columns
  (fits use `pr=TRUE` and `inverseA(nodes="ALL")`). Different topologies give node
  columns that either differ in count (hard error) or silently refer to different
  clades. Fix: when the runs' trees differ, stack fixed effects, VCV, CP and Liab as
  now (their meaning is tree-invariant), match tip-BLUP columns *by name*, drop
  internal-node columns with a one-time message, and set
  `BACE_pooling$trees_vary = TRUE`. Document that the pooled skeleton's `$X`/`$Z`
  (copied wholesale from run 1, `R/pool_posteriors.R:267`) describe tree 1 only.
- `R/with_imputations.R` — `tree =` already anticipates this ("leaves room for a
  future posterior-tree workflow", `R/with_imputations.R:19-22`): accept a list of
  trees, inject tree *m* into fit *m*.

**Does it break anything?** No existing test constrains `class(phylo)`;
"multiPhylo" appears nowhere in R/, tests/, or vignettes. Single-tree calls take the
`n_trees == 1` path and are byte-identical. The reproducibility tests
(`test-final-sweeps.R:131-191`) survive because tree choice consumes no RNG. The
OVR-categorical path flows automatically (its J sub-fits call `.model_fit` with the
same `phylo` argument, `R/model_functions.R:521-539`).

**Workflow impact.** None for single trees. For a `multiPhylo`, the user swaps one
argument: `bace(..., phylo = my_posterior_trees)`. The two-phase workflow gets the
same argument on `bace_final_imp()` directly. Pooled-BLUP caveat goes in the
roxygen and the vignette.

**Tests.**
- CRAN-safe (no MCMC): multiPhylo accepted; per-tree tip-mismatch errors name the
  offending tree; recycling rule when `n_final > n_trees`; `$tree_index` recorded;
  single-tree results unchanged (object-shape assertions only).
- Mocked (`local_mocked_bindings(.model_fit)` spy, the house pattern,
  `test-final-sweeps.R:50-62`): run *m* receives tree `tree_index(m)`.
- Real tiny-MCMC (`skip_on_cran`): 2 very different trees × `n_final = 4`, same
  seed → serial == parallel; datasets from the two tree groups differ.
- Pooling: `pool_posteriors` under differing trees returns, with node columns
  dropped and metadata set; under identical trees, output unchanged.
- Statistical gate (dev/ script, not a unit test): on simulated data with 50
  posterior trees (topological uncertainty concentrated on a clade), total
  imputation variance ≥ single-tree variance, and single-tree runs on each of the
  50 trees bracket the multiPhylo pooled estimate. Pre-registered thresholds in the
  script header, template `dev/15_sweeps_calibration.R:20-43`.

---

## Feature 2 — `bace_validate()` (overimputation diagnostics)

**What.** One exported call that answers "can I trust the imputations on *my*
data": mask a fraction of **observed** cells, run the full pipeline on the masked
data, score recovery against the held-out truth. This replaces a documented manual
recipe — the main vignette teaches exactly this loop by hand
(`vignettes/bace.qmd:708-789`), and both case studies hand-rolled it
(amphibian CV: r/RMSE/coverage95; bird: accuracy/balanced accuracy/RPS/calibration).

**Design decisions.**
- **A fresh run on masked data, never mutation of an existing object.** The
  proper-MI regression test pins that final fits see exactly the originally-missing
  cells as NA (`tests/testthat/test-performance-fixes.R:28-77`); masking must change
  the input data, upstream of `missing_matrix` (`R/bace_imp.R:177-181`), which then
  picks the masked cells up automatically.
- Masking: MCAR v1 (`prop` per modelled variable), stratified for factors so no
  observed level is emptied (guard + warn; a disappearing level changes the model
  space). Never mask the cluster column (hard stop exists, `R/bace_imp.R:106-108`).
  Only variables that appear as a response in some formula can be masked.
- Scoring (per variable, against held-out truth):
  - gaussian/poisson: zRMSE (RMSE / SD of the held-out values — 0 = exact,
    1 = no better than mean-filling, the ms definition), Pearson r, and 95%
    interval coverage from across-run empirical quantiles of the `n_final` values.
  - categorical/threshold: accuracy and balanced accuracy vs the modal-class
    baseline, Brier score from the MI-averaged `prob_preds` (averaging a variable's
    matrices across runs is the documented pooling, `R/bace.R:113-124`).
  - Baselines always reported next to the metric (mean-fill / modal class), so the
    table is interpretable without reading the docs.
- Output: class `bace_validation` — `$table` (variable × metric × baseline),
  `$masked_cells`, `$settings`, `$n_final`; print method with interpretation bands
  (house pattern: `phylo_signal_summary`); `reps` argument for repeated random masks
  (default 1).
- Implementation reuses, moved from dev/ into R/ with tests:
  `dev/benchmark_engine.R` already implements `.apply_mask` (:145), continuous
  metrics (:175-192), majority vote + categorical/ordinal metrics (:197-226). The
  unexported `apply_missingness()` (`R/simulate_auxiliary.r:258-272`) is the MCAR
  primitive.
- Types come from the pipeline's own detection (`.get_type`,
  `R/prep_functions.R:123-168`), never re-detected — the signal-path detector
  classifies double-stored whole numbers differently
  (`R/phylo_signal_summary.R:288-321`), a divergence we must not import.

**Signature.** Mirrors `bace()` explicitly (dots are a dead end):
`bace_validate(fixformula, ran_phylo_form, phylo, data, prop = 0.1, reps = 1,
seed = NULL, ...all bace() args...)`. It forwards to `bace()` internally; `seed`
drives both the masks and the final-phase streams.

**Does it break anything?** No — it is a new export that only consumes the public
pipeline. No existing test asserts the export list.

**Workflow impact.** Additive. The vignette's benchmarking section is rewritten to
one call; recommended before any production run ("validate, then impute").
Honest-cost note in the docs: one full pipeline run per rep; users can pass lighter
MCMC settings for a quick read, with the caveat that coverage reflects those
settings (and ~92% is the n_final = 50 ceiling, `R/bace.R:24-33`).

**Tests.**
- CRAN-safe: masked cells ⊆ observed cells; cluster column untouched; factor-level
  preservation guard fires; metric arithmetic verified on synthetic imputation
  lists with hand-computed answers (no MCMC — feed a fake `bace_complete`).
- Real tiny-MCMC (`skip_on_cran`): end-to-end on an 8–30 tip fixture; output shape;
  metrics within the sanity bounds already used by
  `tests/testthat/test-imputation_accuracy.R:1-30` (r > 0.30, accuracy ≥ modal
  baseline, etc.).
- Edge cases: a variable with too few observed cells to mask → skipped with a
  warning; `reps > 1` returns per-rep and aggregated tables.

---

## Feature 3 — Checkpoint/resume

**What.** `checkpoint_dir` on `bace()`, `bace_imp()` and `bace_final_imp()`.
Final phase: one rds per completed run. Convergence phase: one rds for the finished
`bace` object. Resume = call again with the same arguments and directory; completed
work is loaded, only missing runs execute. This is exactly what both case studies
hand-built (`vignettes/amphibian_ctmax/03_run_full.R:7-57`,
`vignettes/bird_social/02_run_bace.R:44-98`) — the feature deletes that scaffolding.

**Why the design is safe (from the code map).**
- The per-run payload already exists: each run returns
  `list(models, dataset, prob_preds, ovr_models)` (`R/bace_final_imp.R:424-425`),
  and everything after the loop is deterministic positional post-processing
  (`:469-498`) — cached runs splice in by index.
- With `seed`, resume is bitwise-identical to an uninterrupted run (stream-per-run
  property above). **`checkpoint_dir` therefore requires `seed`** — hard error
  otherwise, with a message explaining why (without a seed the serial path consumes
  the global RNG sequentially, so a resume could not reproduce the original draws,
  `R/bace_final_imp.R:222, 462-467`).
- A cached run is final. Warm-up sweep failures legitimately change a run's internal
  RNG consumption (`R/bace_final_imp.R:398-406`), so two executions of the "same"
  run can differ — never re-run-and-compare a cached run.

**Fingerprint.** `manifest.rds` in the checkpoint dir stores the *settings
themselves* (not just a hash, so mismatches can be named): data, tree(s) (feature 1
means trees are per-run inputs; the result object cannot reproduce them), built
formulas, post-transform `ran_phylo_form`, nitt/thin/burnin lists *after*
`nitt_cat_mult` scaling (`R/bace_final_imp.R:139-153`), `eff_sweeps`, priors,
`ovr_categorical`, `species`, `seed`, `n_final`, the `BACE.gelman` option (read per
fit via `getOption`, `R/bace_final_imp.R:327-334`), and the package version. Large
components (data, trees) are compared via base-R hashing: serialize to a temp file,
`tools::md5sum()` — **no new dependency** (digest is absent from DESCRIPTION; the
package already uses undeclared base packages via `::`, `DESCRIPTION:36-44`). On
mismatch: refuse with the differing field names; never silently mix.

**Code changes.**
- `R/bace_final_imp.R` — inside `run_one` (`:429-432`): check cache → load or
  execute → write. Fork-safe write: unique temp file + `file.rename()` (atomic; the
  mclapply workers all write to the same dir, `R/bace_final_imp.R:436`).
- `R/bace_imp.R` / `R/bace.R` — persist the finished `bace` object (the convergence
  phase is unseeded, so it is persisted, never replayed). The auto-retry loop
  re-runs `bace_imp` with `runs × 1.5^attempt` from scratch (`R/bace.R:267-273`);
  each attempt's object overwrites the chain checkpoint with its `runs` recorded,
  so a resumed `bace()` continues the retry ladder instead of restarting it.
- `R/bace.R` — on resume, if `final_results` is incomplete, re-enter the final
  phase before pooling (pooling assumes a complete `final_results`,
  `R/bace.R:385-397`).

**Does it break anything?** Default `checkpoint_dir = NULL` is a no-op; no RNG or
behaviour change. Verbose messaging slots into the existing four step banners
(`R/bace.R:219-326`). Known cost to document, not solve, in v1: per-run files are
heavy because models are saved with `pr/pl/saveX/saveZ = TRUE`
(`R/model_functions.R:129-130`) — the case studies measured ~GB per imputation on
18k rows. v1 documents sizes; a `checkpoint_slim` option is future work (the case
studies prove slim chain objects reload fine,
`vignettes/bird_social/02_run_bace.R:52-74`).

**Workflow impact.** One new argument. The case-study scripts' manual
`saveRDS`/chunking disappear; the bird study's one-run-per-process pattern becomes
`bace_final_imp(..., checkpoint_dir = d, seed = s)` run as many times as needed.

**Tests.**
- The resume oracle (`skip_on_cran`): run A uninterrupted with seed; run B with
  `checkpoint_dir` and a mocked `.model_fit` that fails at run 3 on the first call;
  call B again (mock removed) → `all_datasets` identical to A. This extends the
  existing serial==parallel contract (`test-final-sweeps.R:131-191`).
- CRAN-safe: `checkpoint_dir` without `seed` errors; fingerprint mismatch (changed
  nitt, changed data, changed tree) refuses and names the field; cached files are
  loaded not recomputed (spy counts fits); atomic-write naming; manifest records
  attempt number across the retry ladder.
- Parallel: resume under `n_cores = 2` equals serial resume (`skip_on_os("windows")`
  per house convention).

---

## Feature 4 — Runtime

Four levers, separated because their risk profiles differ. Evidence for where time
goes: fits per final run = `eff_sweeps × n_missing_vars` (+J per OVR categorical)
MCMCglmm fits (`R/bace_final_imp.R:230-245`); bird study: each final run 55–75 min
at nitt = 500k vs a 2.5 h chained phase — the final phase dominates at production
sizes (`vignettes/bird_social/02_run_bace.R:14-20`).

**4a. Cache `inverseA` per tree (do first; hygiene, not a performance claim).**
`.model_fit` recomputes `MCMCglmm::inverseA(tree, nodes = "ALL")` on EVERY fit
(`R/model_functions.R:19`) — `n_final × sweeps × n_vars` (+ OVR levels) identical
computations on an unchanging tree. Measured 2026-11-07: 0.05 s per call at 10,000
tips, so the saving is ~30 s on a case-study run — small, but the hoist is what
gives feature 1 its natural per-tree keying, and it removes the same redundancy in
`phylo_signal_summary` (computes A once, then `.fit_signal_one` recomputes per
variable, `R/phylo_signal_summary.R:229, 415-427`). Hoist into `.get_ainv(tree)`
with a per-call-chain cache (environment keyed by tree identity), computed lazily so
mclapply workers only pay for the trees their runs use. Zero statistical change —
same matrix.
Tests: mocked — `.get_ainv` called once per tree per phase; real tiny-MCMC —
results byte-identical with/without cache under the same seed.

**4b. Final-phase MCMC arguments (explicit, simple).**
`bace()` currently reuses one `nitt/thin/burnin` for both phases
(`R/bace.R:227-231, 299-317`); the case studies hand-tune the final phase an order
of magnitude longer than the chained phase (bird: 50k vs 500k). Add
`final_nitt = NULL, final_thin = NULL, final_burnin = NULL` to `bace()` (NULL =
inherit, exact current behaviour). Pure plumbing; per-formula lists already work
(`.standardize_mcmc_params`, `R/bace_imp.R:9-30`).
Breaks nothing: defaults inherit. Tests: forwarding asserted via spy; defaults
unchanged (reuse existing arithmetic assertions).

**4c. Parallelise the convergence phase across variables (verify first).**
The map found the convergence update is **Jacobi, not Gauss-Seidel**: within a run,
every variable is fit against the *previous* run's data with only its own response
re-NA'd; updates become visible only at the next run
(`R/bace_imp.R:244-248, 356-358`). The per-variable fits within a run are therefore
mutually independent and could run on `n_cores` — the convergence phase is entirely
serial today. Bird mixed: 2.5 h chained → ~/min(n_vars, cores).
*Step zero is a verification test of the Jacobi invariant* (fit order permutation
leaves results unchanged given fixed fits), because this contradicts the usual MICE
mental model and must be true for parallelisation to be valid. Note the final phase
is intentionally Gauss-Seidel (sweeps update `data_current` in place,
`R/bace_final_imp.R:374-376`) — the two phases differ and that difference is
correct (convergence wants settling diagnostics; final wants a joint draw).
The convergence phase is unseeded, so there is no reproducibility contract to
preserve; the macOS-BLAS fork caveat gets the same fallback the final phase has.
Separate PR; off by default (`n_cores_conv = 1`).

**4d. ESS-targeted auto-sizing (benchmark-gated; may not land).**
Opt-in `target_pooled_ess = NULL` on `bace_final_imp()`. When set: use the
convergence phase's last-run models to estimate per-variable ESS-per-iteration
(diagnostics machinery exists: `.diagnose_one_model`,
`R/model_functions.R:1629-1650`; pooled draws = `n_final × (nitt−burnin)/thin`,
`R/pool_posteriors.R:167`), then resize `nitt_list` so pooled ESS ≈ target. Resolved
lists go into the checkpoint manifest (interaction with feature 3).
**Must be opt-in**: five test files pin the `(nitt−burnin)/thin` arithmetic
(`test-mcmc_param_lists.R:249-282`, `test-final_imp_extras.R:128-147`,
`test-posterior_sampling.R:355-363`, `test-bace_wrapper.R:196-203`,
`test-final-sweeps.R:64-95`) and explicit `nitt` must always be honoured.
Gate: a dev/ calibration study on the `dev/15_sweeps_calibration.R` template —
pooled estimates and CI widths within Monte-Carlo error of the fixed-length
reference, per-cell interval coverage unchanged, wall-time saving reported.
If the gate fails, the feature is dropped and the study still documents why.

**Not a lever: `n_final`.** Per-cell prediction intervals are empirical quantiles
across the `n_final` values (50 → ~92% effective coverage, `R/bace.R:24-33`);
trading `n_final` for chain length changes interval calibration.

---

## Cross-feature interactions (so nobody trips later)

- **1 + 3:** the fingerprint includes the tree list; each cached run records its
  `tree_index`. A cache can never silently pair a run with a different tree.
- **1 + 4a:** the Ainv cache is keyed per tree — multiPhylo gets one entry per tree.
- **1 + pooling consumers:** `pool_posteriors` consumes only `$all_models`
  (`R/pool_posteriors.R:108-315`); `with_imputations`/`pool_mi` consume only
  `$all_datasets` — the dataset-pooling route is tree-agnostic and needs no change.
  OVR categoricals are invisible to `pool_posteriors` by design (per-variable slot
  NULL, `R/bace_final_imp.R:312-318`) — unchanged.
- **4b/4d + 3:** chain lengths enter the fingerprint *after* resolution (post
  `nitt_cat_mult`, post auto-size), so resumes compare what actually ran.
- **2 + everything:** `bace_validate()` is a pure consumer; it inherits multiPhylo,
  checkpointing (useful for expensive validations) and the final-phase arguments
  for free once they exist.

## What this does NOT change

- The two-phase contract: convergence predicts with `sample = FALSE` (point
  estimates settle; `R/bace_imp.R:305, 331-342`), final draws with `sample = TRUE`.
  Untouched by all four features.
- Default MCMC settings, priors, `n_final = 50`, `sweeps = 3`.
- The `bace_complete`/`bace_final`/`bace_pooled` object shapes — only *additive*
  fields (`$tree_index`, checkpoint/provenance metadata); no existing field moves.
- CRAN v0.1.0 submission content (on hold awaiting Dan's go). All four features are
  post-submission (0.2.0 material) unless Dan decides otherwise.

## Effort estimate (coding + tests, excluding gate studies)

| Feature | Size | Risk |
|---|---|---|
| 4a inverseA cache | ~½ day | low |
| 1 multiPhylo | ~2–3 days (pooling alignment is most of it) | medium |
| 3 checkpoint/resume | ~2–3 days (fingerprint + retry-ladder semantics) | medium |
| 4b final_* args | ~½ day | low |
| 2 bace_validate | ~2 days (port + metrics tests + vignette rewrite) | low |
| 4c conv parallelism | ~1–2 days after the Jacobi verification | medium |
| 4d ESS auto-sizing | ~2 days + calibration study | high (gated) |
