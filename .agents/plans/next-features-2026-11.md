# Plan: four post-submission features (2026-11)

Dan picked these four (2026-11-07) from the "what would make BACE better" discussion:

1. **Tree-per-imputation** — accept a `multiPhylo` and use tree *m* for final run *m*,
   so between-imputation variance includes phylogenetic (tree) uncertainty.
2. **`bace_validate()`** — mask a fraction of *observed* cells, re-impute, score
   recovery on the user's own data (zRMSE, interval coverage, category accuracy).
3. **Checkpoint/resume** — long runs (22–28 h in our case studies) survive crashes;
   per-final-run caching with a settings fingerprint.
4. **Runtime** — attack where the time actually goes (final phase dominates), without
   touching `n_final` (load-bearing for per-cell interval calibration, `R/bace.R:24-33`).

## Provenance

Built on a full code map (2026-11-07, four parallel readers; every claim anchored to
`file:line`), then adversarially reviewed twice:

- **Round 1** (compatibility, statistics, tests): confirmed the mechanical core;
  found 7 must-fixes (multiPhylo's opaque death in `inverseA`; pooling's silent
  tip-misalignment mode; a zRMSE spec off by ~√2 from the ms; a gate blind to tree-1
  initialization leakage; two test designs that could not detect their targets; a
  missing MNAR caveat). All folded in.
- **Round 2** (internal consistency, statistics-by-simulation, implementation-by-
  experiment, citations-by-web): confirmed ~40 anchors and reproduced every measured
  claim; found 4 further must-fixes (fingerprint omitted the converged dataset;
  validator omitted ultrametricity; formula serialization captures environments;
  the gate's acceptance rule was statistically underpowered) plus citation
  corrections. All folded in.

**Verified-OK roster** (round 2 ran these; round 3 need not re-litigate):
`serialize(version = 2)` has an exactly-14-byte header carrying the writer R
version (strip-and-hash is version-stable); `%||%` already exists package-wide
(`R/with_imputations.R:158`) and in base R ≥ 4.4; mclapply forks see parent
environments copy-on-write and worker writes never propagate back;
`file.rename` returns FALSE across filesystems (verified against a ram-disk);
`inverseA` costs 0.016–0.05 s at 100–10,000 tips (~2 MB per Ainv); the per-cell
Rubin t-interval is numerically nominal for gaussian at any M ≥ 5 (simulated
0.9498 at M = 50); the recycling deflation arithmetic checks out; the zRMSE and
t-interval specs match `ms/ms.qmd:308-330` and `rubin_lib.R:80-95` exactly.

## Ground rules that shaped every design decision

- **New arguments must be explicit.** `bace()` warns on any `...` and the workers
  silently drop dots (`R/bace.R:166-178`, `R/bace_imp.R:83`). Every new argument is
  added to the signature of `bace()` *and* the relevant worker, plus roxygen.
- **Features must land at the `bace_imp()`/`bace_final_imp()` level, with `bace()`
  forwarding.** Both case studies bypass `bace()` and call the phases directly
  (`vignettes/amphibian_ctmax/00_setup.R:33-48`,
  `vignettes/bird_social/02_run_bace.R:44-57`).
- **`bace()` calls `bace_final_imp` at TWO sites** — converged (`R/bace.R:299-317`)
  and non-converged/skip_conv (`R/bace.R:350-368`). Every new forwarded argument is
  plumbed through *both*, and **every feature that adds one carries a forwarding
  test that runs with `skip_conv = TRUE`** (features 1, 3 and 4b each list one).
- **`.model_fit`'s signature is frozen** at its current 9 arguments
  (`R/model_functions.R:16`): one test file mocks it three times with the full
  signature (`test-final-sweeps.R:52, 210, 242`; two other files use partial
  mocks). Per-run trees arrive through the *existing* `tree` argument.
- **Default paths stay byte-identical.** The bird vignette pins results to commit
  `dc79be1` with recorded seeds; the serial==parallel and fallback-stream tests
  (`tests/testthat/test-final-sweeps.R:131-191`) define the RNG contract. No new
  RNG consumption on any default path; no default value changes.
- **The final-phase RNG design is our friend.** With `seed` set, run *m*'s entire
  randomness comes from a stream computed from `(seed, m)` alone
  (`R/bace_final_imp.R:180-207, 221-224`); runs re-seed at entry, so any subset can
  execute in any order bitwise-identically (the failed-worker retry already relies
  on this, `:52-56, 448-461`).
- **The convergence phase is unseeded by design** (`seed` reaches only the final
  phase, `R/bace.R:80-82`) — it is persisted as an object, never replayed. This has
  a hard consequence used repeatedly below: *two invocations of the chained phase
  under identical arguments produce different converged datasets*, so anything that
  must be comparable across invocations must share ONE converged `bace` object.

## Speed is a first-class constraint (Dan, 2026-11-07)

Accuracy and calibration come first, but BACE's runtime (22–28 h in the case
studies) means **no feature may slow the default path, and every PR reports
wall-time** on a fixed reference fixture (hook: the existing `benchmark.yml` /
`sim-reference-benchmark.yml` GHA workflows and `dev/benchmark_engine.R`). Where
the time goes, measured: the final phase dominates (bird: 20 × 55–75 min final runs
vs 2.5 h chained) — `n_final × sweeps × n_missing_vars` full MCMC fits. The levers
that move hours are 4b/4d (chain sizing), 4c (parallelising the chained phase), and
the already-shipped `n_cores`.

**The multiPhylo speed hazard was measured twice and is negligible**:
`inverseA(tree, nodes = "ALL")` costs 0.05 s at 10,000 tips (~2 MB per Ainv), so 50
posterior trees add ~2.5 s and ~110 MB to a multi-hour run; the MCMC fits are
indifferent to which tree they receive. The same measurement demotes 4a to hygiene
(~30 s saved per case-study run) — its real value is the per-tree keying feature 1
needs. **Checkpoint writes are the other measured cost**: gzip `saveRDS` on
MCMC-like payloads ran 28× slower for a 4% size saving (1.65 s vs 0.06 s on 49 MB),
so per-run checkpoint files default to `compress = FALSE` (flag provided), and
feature 3's PR benchmarks the write path at case-study scale.

## Order of work

```
4a  inverseA cache            (tiny, standalone, zero statistical change)
1   multiPhylo                (touches the final-run loop)
3   checkpoint/resume         (same loop; after 1 so the fingerprint includes
                               trees from day one)
4b  final-phase MCMC args     (plumbing; what the case studies already hand-roll)
2   bace_validate()           (a consumer of everything above)
4c  convergence-phase         (Jacobi confirmed; redesigned step-zero test first;
    parallelism                separate PR)
4d  per-run ESS floor         (benchmark-gated; separate PR, may not land)
```

Each lands as its own PR with its own tests; `devtools::test()` green and
`R CMD check --as-cran` clean are gates for every one. NEWS.md bullet per feature
(house style, `NEWS.md:1-32`). When feature 1's rationale moves into roxygen/ms,
the same commit adds the verified citations to `bib/BACE.bib` (none are in it yet
except `nakagawa2018` — whose year field must be corrected to 2019 when cited).

---

## Feature 1 — Tree-per-imputation (`multiPhylo`)

**What.** `phylo` accepts a `multiPhylo` (ape class; ape already in Imports). The
convergence phase uses `trees[[1]]`; final run *m* uses a tree chosen
**deterministically from the run index, never drawn** (no RNG consumed; the stream
contract survives):

- `n_trees ≥ n_final`: an evenly spaced subsequence
  `trees[round(seq(1, n_trees, length.out = n_final))]` — verified duplicate-free
  for all (n_final, n_trees) pairs up to 120×400. Posterior tree files are
  chain-ordered, so this is thinning, which suffices as rationale. Docs require
  post-burn-in posterior samples (LogCombiner/TreeAnnotator output) and tell users
  to drop any leading summary/consensus tree before passing the set (tree 1 also
  anchors the convergence phase).
- `n_trees < n_final`: recycle `((m-1) %% n_trees) + 1`. `K ≥ n_final` is NOT
  required: full accounting (round-2 simulation) — balanced reuse at
  `n_final = 50, K = 20` deflates E[B]'s tree component by ~3.3% AND under-inflates
  the (1+1/M) correction by ~5.2%, together ~6% of the *tree-variance component*
  (≈2% of total variance when trees contribute 30% of it). Documented with the
  advice to prefer `K ≥ n_final` when trees are plentiful. The formally exact
  framework for reuse is nested/two-stage MI — Shen (2000, Harvard PhD thesis);
  Rubin (2003, *Statistica Neerlandica* 57:3-18, single-authored); cf. Harel &
  Schafer (2003) — which plain Rubin approximates adequately at this deflation.

**Selection is single-sourced in `bace_final_imp()`**: `bace()` forwards the *full*
multiPhylo; the worker selects and records `$tree_index` as indices into the
*user's original object* (if `bace()` pre-selected, the worker's indices would
collapse to 1:n_final and provenance would be lost).

Statistical rationale (for roxygen and, later, the ms): proper MI with the
phylogeny as a nuisance parameter drawn from its (empirical) posterior — Rubin
(1987, ch. 3–4). **Primary precedent: Nakagawa & de Villemereuil (2019, *Syst
Biol* 68(4):632-641, doi 10.1093/sysbio/syy089)** — treats true trees as missing
data, one analysis per tree, pooled by Rubin's rules, and finds ~50 trees suffice;
already in `bib/` as `nakagawa2018` (fix year to 2019). Tree-sample averaging
precedent more broadly: Huelsenbeck, Rannala & Masly 2000 *Science* 288:2349;
Housworth & Martins 2001 *Syst Biol* 50:628; de Villemereuil et al. 2012 *BMC Evol
Biol* 12:102; Rangel et al. 2015 *Evolution* 69:1301; software precedent *mulTree*
(Guillerme & Healy). **The convergence phase conditions on a single tree by
necessity, not convenience**: the stabilisation diagnostic reads the settling of
deterministic point imputations (`sample = FALSE` invariant,
`R/bace_imp.R:331-342`); cycling trees would inject tree-to-tree jitter into
exactly the trajectories it needs stable. Tree uncertainty is propagated where it
belongs — the between-imputation variance of the final draws.

**Code changes.**
- **Shared validator `.validate_trees()`**, called at the top of `bace_imp()` AND
  `bace_final_imp()` (`bace()` inherits it). Today a multiPhylo dies inside
  `inverseA` with `unused argument (2)` (reproduced twice). Per tree, it:
  1. checks class and that tip labels match the data's cluster values;
  2. **strips node labels unconditionally** (with a message). This kills three
     birds: `inverseA` hard-stops on duplicated node labels (e.g. BEAST support
     values — reproduced); differing node labels would otherwise fragment the 4a
     cache key (collision reproduced: identical trees with different node labels
     give identical Ainv values but different rownames); and pooling's
     node-column identification simplifies to the `Node<k>` complement;
  3. **runs `.get_ainv()` (feature 4a) on the tree up front** — this is the
     airtight check: it surfaces *every* condition `inverseA` enforces, most
     importantly **ultrametricity** (`scale = TRUE` default at both production
     call sites; round-2 reproduced that even a `write.tree(digits = 6)`
     round-trip of a genuinely ultrametric tree fails it — i.e. exactly the
     posterior files this feature ingests), names the offending tree index
     *before* hours of fitting, guarantees tolerance-consistency with the fit
     path by construction, and pre-warms the Ainv cache in the parent before
     workers fork (copy-on-write). Cost: ~2.5 s for 50 large trees. Also checks
     zero/negative branch lengths with a clear message (cheaper to say than to
     let inverseA's error surface);
  4. warns on **scale inconsistency** (across-tree root-depth CV > ~10–25%):
     genuine BEAST root-height variation passes; an accidentally rescaled tree
     does not. Trees must come from one analysis in common units.
- `R/bace.R` — accept phylo|multiPhylo; `trees[[1]]` to `bace_imp()`; the full
  multiPhylo to `bace_final_imp()` at **both** call sites; `phylo_signal = TRUE`
  preview runs on `trees[[1]]` with a message.
- `R/bace_imp.R` — accept multiPhylo by taking `trees[[1]]` (validator + message).
- `R/bace_final_imp.R` — select per-run trees; thread through
  `run_one`/`.one_final_run` to the existing `tree`/`phylo` arguments of
  `.model_fit` (:350) and `.fit_predict_ovr` (:301). Record `$tree_index` and a
  per-tree hash (via the stripped-serialize hash, feature 3's helper) on the
  `bace_final` output.
- `R/pool_posteriors.R` — the real hazard (round-1 reproduced): same-tip trees
  give **identical Sol column-name sets in different orders**, so positional
  stacking silently misaligns tip BLUPs; nothing ever hard-errors. Fix, with the
  mechanics round 2 demanded: when `$tree_index` says trees differ (column names
  cannot reveal it), stack fixed effects, VCV, CP, Liab as now; identify each
  run's random columns via its own `colnames(model$Z)` (the order-aligned tail of
  Sol under `pr = TRUE`); classify tip vs node by **exact match of the level part
  against the known tip set, per random term** (prefixes differ: `animal.`,
  `Species2.` under `species = TRUE`, trait prefixes for categorical — and tip
  labels may contain dots, so never `strsplit` on "."); align tip columns by name;
  drop node columns; set `BACE_pooling$trees_vary = TRUE`; document that the
  pooled skeleton's `$X`/`$Z` (copied from run 1, `:267`) describe tree 1 only.
- `R/with_imputations.R` — accept a multiPhylo or bare list of phylo: tree *m*
  injected into fit *m*, **recycled by the same `((m-1) %% K) + 1` rule as
  bace_final_imp** (decided: recycling, not strict length validation). A single
  phylo keeps today's same-tree-for-all behaviour. When `object` is a
  bace_final/bace_complete carrying `$tree_index`, verify the supplied list is
  consistent with it and message otherwise. Escape hatch documented for an `.f`
  whose own `tree` argument expects a list: wrap it
  (`.f = function(d, tree) g(d, my_trees)`).

**Does it break anything?** No test constrains `class(phylo)`; "multiPhylo" appears
nowhere in R/, tests/, vignettes. Single-tree calls take the `n_trees == 1` path,
byte-identical. Tree choice consumes no RNG. OVR flows automatically. The existing
tip-mismatch test is `expect_error` with no regexp — new messages can't break it.

**Workflow impact.** None for single trees; one argument swap for a multiPhylo, in
either workflow. Pooled-BLUP caveat in roxygen and vignette.

**Tests.**
- CRAN-safe (no MCMC): class acceptance; per-tree tip-mismatch naming the tree,
  against BOTH entry points; node-label stripping message; **non-ultrametric tree
  errors at entry naming the tree** (the `write.tree` round-trip case); negative
  branch length errors; scale-consistency warning; subsequence/recycling
  arithmetic via an exported-internal helper; `with_imputations` tree-list tests
  with plain `lm` (recycling rule, single-tree unchanged, `.f` without a tree
  formal, `$tree_index` consistency message).
- Tiny-MCMC (`skip_on_cran`): `$tree_index` and per-tree hashes recorded;
  single-tree results unchanged; `phylo_signal` preview message; **forwarding spy
  under `skip_conv = TRUE`** (trees reach `bace_final_imp` at the second call
  site).
- **The identity test** (round-1 redesign, round-2 correction): **all three arms
  run `bace_final_imp` on ONE shared converged `bace` object** (chained phase run
  once, on tree A — three separate pipeline invocations would each produce a
  different unseeded converged state and the equalities could never hold). With
  seed *s*: multiPhylo(A,B) `n_final = 4` vs single-tree-A vs single-tree-B on the
  same object; assert multi runs 1,3 ≡ single-A runs 1,3; multi runs 2,4 ≡
  single-B runs 2,4; single-A ≢ single-B. Catches "tree silently ignored" and
  stream contamination.
- Serial == parallel under multiPhylo (`skip_on_os("windows")`).
- Pooling: CRAN-safe via a handcrafted fake `bace_final` whose per-run Sol
  matrices have identical name sets in permuted order, **including one dotted tip
  label (e.g. "G.fortis") and a `species = TRUE`-style second random term**;
  assert by-name tip alignment, node dropping, `trees_vary` driven by
  `$tree_index`.
- **Statistical gate** (dev/ script; redesigned twice). Purpose: detect tree-1
  initialization leakage (every final run starts from the tree-1 converged state;
  residual dependence after 3 sweeps would shrink the between-tree variance
  component — anti-conservative in exactly the direction the feature exists to
  fix). Round-2 power analysis showed my "within ~10%/MC error" rule was
  ill-defined and underpowered (a TOST at a 10% margin fails a PERFECT
  implementation ~98% of the time at 25 seeds; a naive mean-ratio rule passes a
  10%-deflated one 43% of the time). The redesigned gate:
  1. **Gold-standard arm with replication**: ≥2 independent chained phases per
     tree (the phase is unseeded, so single chained states would confound tree
     separation with chain Monte-Carlo noise); ANOVA decomposition estimates the
     within-tree between-replicate variance and subtracts it before comparing
     between-tree components. The within-tree replicate spread is itself the
     direct answer to "is initialization forgotten?"
  2. **One-sided, pre-registered decision rule** on the quantity that matters
     (deflation of B), sized to the rule: SE of the mean log-ratio ≤ ~0.029 →
     ~100 seeds at `n_final = 50` for 90% pass-when-true at a 10% margin (or ~35
     seeds at a 15% margin); per-seed noise reduced by estimating B from per-cell
     components where possible. Thresholds in the script header
     (`dev/15_sweeps_calibration.R` template).
  3. **Sweeps dose-response** (sweeps ∈ {1, 2, 3, 5, 10}, paired seeds):
     tree-group separation and B must plateau by sweeps = 3, else the default is
     raised for multiPhylo runs.

---

## Feature 2 — `bace_validate()` (overimputation diagnostics)

**What.** One exported call answering "does the imputation model predict well where
data were observed": mask a fraction of **observed** cells, run the full pipeline
on the masked data, score recovery against the held-out truth. Replaces the manual
recipe the vignette teaches (`vignettes/bace.qmd:708-789`) and both case studies
hand-rolled.

**The honesty caveats go in roxygen AND the print method.** (1) MCAR masking
validates predictive accuracy *on the observed-data region*: necessary, not
sufficient, for trusting imputations of the actually-missing cells; it cannot
detect MNAR bias (the genuinely missing species are often small/rare — an
extrapolation region). (2) Phrase the precedent as **"inspired by / cf. Amelia's
`overimpute()`"** (Honaker, King & Blackwell 2011, *J Stat Softw* 45(7)) and say
explicitly that Amelia overimputes *sequentially, one cell at a time*, while
`bace_validate()` masks a fraction *simultaneously* — deliberately harsher, so the
two diagnostics will not produce identical numbers. (The Amelia paper itself: "it
cannot speak to the untestable assumption of MAR" — same framing.) MAR-style
stress masking (mask probability depending on observed covariates; pattern in
`dev/15_sweeps_calibration.R:14-16`) is the documented v2 extension.

**Design decisions.**
- **A fresh run on masked data, never mutation of an existing object** (pinned by
  `test-performance-fixes.R:28-77`); masking changes the input data upstream of
  `missing_matrix` (`R/bace_imp.R:177-181`).
- Masking: MCAR v1, `prop` per modelled variable, stratified for factors so no
  observed level empties (guard + **warning** — the dev masker silently skips
  <10-cell variables, `dev/benchmark_engine.R:153`; the port warns). Never the
  cluster column; only variables with a formula. Per-rep mask seeds derived
  deterministically from `(seed, rep)` **outside** the pipeline call. If
  `seed = NULL` (and no checkpointing), mask seeds are drawn from the session RNG
  and **recorded in `$settings`** so the masks remain reportable.
- **zRMSE — the manuscript's definition exactly** (verified against
  `ms/ms.qmd:308-313` and `rubin_cell.R:110-114`): numerator = RMSE of the
  across-`n_final` **mean** at masked cells; denominator = SD of the variable's
  **remaining observed cells**. Do not port `.metrics_continuous` as-is (per-draw
  RMSE ≈ √2 larger; different estimand).
- **Interval coverage — per-cell Rubin t-interval**
  (`ybar ± t_{M-1,0.975} · sqrt((1 + 1/M)·s²)`, matching `ms/ms.qmd:324` and
  `rubin_lib.R:80-95`; simulated nominal: 0.9498 at M = 50 for gaussian).
  Round-2 refinements:
  - **Ship an exported helper** computing per-cell Rubin t-intervals from
    `$all_datasets` (≈5 lines; exists as dev code `rubin_lib.R`) and point the
    `n_final` roxygen at it — otherwise validate certifies an interval the
    package's own docs don't teach (the documented empirical-quantile recipe
    yields ~91–92% at M = 50, and users would wrongly read a 3-point gap as
    miscalibration).
  - **Rate-aware expectation for counts**: discreteness pulls true calibration of
    the t-interval to ~0.93–0.94 at mean rates below ~0.6 (simulated; consistent
    with the campaign's own 0.943–0.964). When a poisson variable's mean imputed
    rate < ~2, print a simulated-under-calibration expectation next to the
    observed coverage (one vectorised rpois simulation, milliseconds).
  - **The MC band**: per-cell indicators are positively correlated (all cells
    share the same M fits), so a binomial band is anti-conservative (simulated:
    true SD 1.6–3.7× binomial; 14–50% false alarms on a perfect imputer). At
    `reps = 1`, label the binomial band an *optimistic lower bound* (or use a
    jackknife over imputation runs, which captures the shared-fit component); at
    `reps ≥ 2`, the usable band comes from between-rep variation.
  - Zero-width intervals (degenerate low counts) flagged.
- **Poisson noise floor — two reference lines, precisely worded** (round 2: my
  single plug-in floor sits at ~half the achievable floor in weak-signal regimes
  and would make optimal imputers look bad): (i) `sqrt(mean imputed rate)/s_obs`
  as the **absolute sampling-noise bound** ("no imputer can beat this"); (ii) the
  **model-implied calibration value** `sqrt((1 + 1/M) · mean per-cell
  variance)/s_obs` ("a calibrated imputer should score near this") — both
  computable from the imputations; (ii) is what "extracting all available signal"
  actually looks like.
- **Proportion-scale note**: the ms clips proportion traits to [1e-4, 1−1e-4] and
  scores on the logit scale (`ms/ms.qmd:313`, `rubin_lib.R::score_scale`); BACE
  sees 0–1 variables as plain gaussian. v1 ships an optional per-variable
  `scale = c("raw", "logit")` argument (default raw) with a documented note that
  raw-scored proportions are not comparable to the published figures.
- Categorical/threshold: accuracy and balanced accuracy vs modal baseline;
  **Brier from the MI-averaged `prob_preds`** (the documented pooling,
  `R/bace.R:113-124`) — not from imputed-class indicator proportions.
- **A named scorer seam**
  `.score_validation(truth, masked_cells, imputed_datasets, prob_preds, types)` —
  testable against a handcrafted fake without MCMC. Types come from the pipeline's
  own `$types`, never re-detected.
- Output: class `bace_validation` — `$table` (variable × metric × baseline ×
  expectation/floor), `$masked_cells`, `$settings` (incl. recorded mask seeds),
  `$n_final`; print with interpretation bands; `reps` (default 1) with per-rep and
  aggregated tables and per-metric MC error.

**Signature.** Mirrors `bace()` explicitly — *including* `skip_conv` and
`max_attempts` (the retry ladder multiplies per-rep cost; the cost note says so).
**Checkpointing composes via per-rep subdirectories**
(`file.path(checkpoint_dir, sprintf("rep%03d", rep))`): each rep's masked data
differs by design, so a shared directory would be *refused* by feature 3's own
fingerprint rule — "inherits checkpointing" means explicit mirrored arguments plus
this subdirectory rule, not literally for free.

**Tests.**
- CRAN-safe: masked ⊆ observed; cluster untouched; level-preservation and
  <10-cells warnings; mask reproducibility (same seed ⇒ same mask; reps differ;
  `seed = NULL` records drawn seeds); `.score_validation` against a handcrafted
  fake with **hand-computed zRMSE (mean-filler ≈ 1, perfect = 0), t-interval
  coverage, both poisson reference lines, and Brier** — including a
  replicated-species case (make.unique rownames).
- Tiny-MCMC (`skip_on_cran`): structure-only assertions on the small fixture.
  Metric-level gates live in a dev/ script on a ≥60-species fixture with
  `prop ≥ 0.3` (the conditions where the borrowed bounds are stable; an 8-tip
  fixture at `prop = 0.1` masks 1–3 cells and degenerates).

---

## Feature 3 — Checkpoint/resume

**What.** `checkpoint_dir` on `bace()`, `bace_imp()` and `bace_final_imp()`. Final
phase: one rds per completed run. Convergence phase: one rds for the finished
`bace` object. Resume = same arguments + same directory; completed work loads, only
missing runs execute.

**Who writes what (round-2 decision).**
- `bace_imp(checkpoint_dir)` writes its finished `bace` object. **No seed
  requirement here** — the phase is unseeded by design and is persisted, never
  replayed. (`runs` is recoverable as `length(out$data) - 1`; the object itself is
  the checkpoint.)
- `bace_final_imp(checkpoint_dir, seed)`: **seed required at THIS level** (hard
  error without it; with `seed = NULL` the serial path consumes global RNG and a
  resume could not reproduce the original draws, `R/bace_final_imp.R:222,
  462-467`). Per-run rds written in `run_one` (`:429-432`): check cache → load or
  execute → write.
- `bace()` owns the **retry-ladder metadata**: after each `assess_convergence` it
  records `{attempt, runs_actually_used}` in the manifest (bace_imp cannot know
  the attempt; `match.call()` records the symbol `runs_new`, not its value). The
  manifest treats runs/attempt as *ladder-mutable* fields — a resumed `bace()`
  continues the ladder rather than refusing or restarting. On resume with
  `final_results` incomplete, `bace()` re-enters the final phase before pooling
  (pooling is called at `R/bace.R:328` and `:374`).

A cached run is final — warm-up sweep failures legitimately alter a run's internal
RNG consumption (`:398-406`), so cached runs are never re-run-and-compared.
Migration note: the bird case study already passes `seed`; the amphibian helper
seeds externally *without* the seed argument
(`vignettes/amphibian_ctmax/00_setup.R:77-81`) — migrating it changes its seeding
semantics, and its archived draws are not bitwise-recoverable under the stream
contract. Say so in its README when migrating.

**Fingerprint (manifest.rds).** Stores the settings *themselves* so mismatches can
be named. Fields: **the supplied `bace_object`'s last dataset under its own named
field** (round-2 must-fix: every final run conditions on it,
`R/bace_final_imp.R:120`; the chained phase is unseeded, so re-running a crashed
two-phase script regenerates a *different* converged state under identical
arguments — without this field a resume would silently mix runs from two
realizations; the two-phase docs direct users to checkpoint `bace_imp` as well so
the converged state itself is pinned); user data; tree(s) + `tree_index` +
per-tree hashes; **formulas and post-transform `ran_phylo_form` stored and
compared as `deparse1()` strings, never serialized** (round-2 must-fix,
reproduced: `serialize()` follows a formula's environment — a formula built inside
a function dragged 4 MB of enclosure into the hash and made content-identical
formulas compare unequal); nitt/thin/burnin lists AFTER `nitt_cat_mult` (and after
any 4d resize); `eff_sweeps`; priors; `ovr_categorical`; `species`; `seed`;
`n_final`; the `BACE.gelman` option; package version; **R version as its own named
field**. Large objects (data, trees, priors, numeric lists) are hashed as
`serialize(version = 2)` **minus its fixed 14-byte header** (verified: the header
carries the writer's R version; stripping it makes hashes version-stable) via
`tools::md5sum` on a temp file in node-local `tempdir()`. On mismatch: refuse,
naming the differing fields; never silently mix. **Documented recovery for the
version-refusal case** (the most likely crash cause is a BACE bug): reinstall the
pinned version (`remotes::install_github("...@<sha>")` — the pattern the bird
vignette already teaches) to finish the run under the version that started it;
the strict default stands.

**Write mechanics.** Fork-safe atomic write:
`tempfile(tmpdir = checkpoint_dir)` then `file.rename()` **with its return value
checked** (FALSE across filesystems — verified; unchecked, the checkpoint is
silently lost); rename atomicity on NFS is best-effort, documented.
`checkpoint_compress = FALSE` default for per-run files (measured: 28× time for 4%
size on MCMC payloads). Resume sweeps stale temp files older than the manifest,
with a message. Known documented cost: per-run files are heavy (models carry
`pr/pl/saveX/saveZ = TRUE`; ~GB per imputation on 18k rows); `checkpoint_slim` is
future work (slim reload proven viable, `vignettes/bird_social/02_run_bace.R:52-74`).

**Does it break anything?** Default `checkpoint_dir = NULL` is a no-op; no RNG or
behaviour change.

**Tests.**
- **The resume oracle** (round-1 redesign): fixture with ONE missing variable so
  `eff_sweeps` collapses to 1 (`:177-178`) and every fit is fatal; `n_cores = 1`
  (mock counters are serial-only); run A uninterrupted with seed; run B with
  `checkpoint_dir` and a mock failing the full-nitt fit of run 3 (house pattern
  keys on nitt), wrapped in `expect_error`; assert runs 1–2's cache files exist;
  resume → `all_datasets` identical to A.
- Parallel resume as a real-fit test (`skip_on_cran`, `skip_on_os("windows")`):
  resume under `n_cores = 2` equals serial resume.
- CRAN-safe: `checkpoint_dir` without `seed` errors **at bace_final_imp/bace()**
  (and is accepted by `bace_imp` without seed); fingerprint mismatch on
  handcrafted manifests (changed nitt / data / tree / **converged dataset** /
  package version / R version) refuses naming the field; **formula environments do
  not affect the fingerprint** (two content-identical formulas built in different
  environments match); cached files loaded not recomputed (serial spy);
  atomic-write naming; **a mocked `file.rename` returning FALSE errors rather than
  losing the run**; `compress = FALSE` round-trip.
- **Retry-ladder resume** (tiny-MCMC, mocked `assess_convergence`): interrupt
  after attempt 2's chain checkpoint, resume, ladder continues from attempt 2.
- **bace()-level end-to-end** (`skip_on_cran`): `bace(seed, checkpoint_dir)`
  interrupted mid-final-phase resumes and **pools identically** to uninterrupted.
- **Forwarding spy under `skip_conv = TRUE`**: `checkpoint_dir` reaches
  `bace_final_imp` at the second call site.

---

## Feature 4 — Runtime

**4a. Cache `inverseA` per tree (hygiene + the keying feature 1 needs).**
`.model_fit` recomputes `inverseA(tree, nodes = "ALL")` every fit
(`R/model_functions.R:19`); measured saving ~30 s per case-study run. Structure:
`.get_ainv(tree)` delegating to `.compute_ainv(tree)`; cache in a package-local
environment, **cleared at entry of `bace_imp()`, `bace_final_imp()` and
`phylo_signal_summary()`** (round 2: an uncleared namespace cache leaks ~2 MB per
distinct tree per session; clearing at phase entry bounds memory, makes the spy
assertion true in serial, and still pre-warms workers — they fork *after* entry
and inherit the warm cache copy-on-write; worker writes never propagate back, so
under parallel the true behaviour is **once per tree per worker** — documented,
and cheap at 0.05 s). **Cache key: hash of edge matrix + edge lengths + tip labels
+ node labels + Nnode** (round 2 reproduced a collision on node labels alone —
same values, different Sol rownames; moot after the validator strips node labels,
but `phylo_signal_summary` can receive unstripped trees, so the key stays
complete). Also removes the same redundancy in `phylo_signal_summary`
(`:229, 415-427`). Zero statistical change.
Tests: inner-function (`.compute_ainv`) spy, `n_cores = 1`, asserting **at most
once per tree per phase (serial)**; collision test — same tips, different
topologies AND same topology, different node labels → distinct cache entries,
correct Ainv per run; byte-identical results with/without cache under one seed.

**4b. Final-phase MCMC arguments.**
`bace()` reuses one `nitt/thin/burnin` for both phases; the case studies hand-tune
the final phase an order of magnitude longer (bird: 50k chained vs 500k final).
Add `final_nitt/final_thin/final_burnin = NULL` to `bace()`; **resolution lives in
`bace()`** (`final_nitt <- final_nitt %||% nitt` — `%||%` already exists
package-wide at `R/with_imputations.R:158` — before BOTH `bace_final_imp` call
sites; `.standardize_mcmc_params(NULL)` errors loudly, so a slipped NULL cannot be
silent, but the substitution still belongs in `bace()`).
Breaks nothing: NULL defaults inherit byte-identically (both workers standardise
and apply `nitt_cat_mult` internally and identically). Tests: forwarding spy incl.
the `skip_conv = TRUE` branch; list-valued `final_nitt` through `bace()`;
`final_nitt` × `nitt_cat_mult` interaction (`R/bace_final_imp.R:144-153`);
defaults unchanged (the genuine pin is `test-posterior_sampling.R:355-363`;
`test-bace_wrapper.R:196-203` self-adjusts and guards nothing).

**4c. Parallelise the convergence phase across variables.**
**Jacobi confirmed twice**: within a run, every variable is fit against the
previous run's dataset with only its own response re-NA'd
(`R/bace_imp.R:233, 244-248`); updates land in a write-only copy (`:356-358`)
published at run end (`:366`). The diagnostic consumes only the dataset sequence
(`R/convergence_diagnostics.R:59`), order-invariant given Jacobi; the phase is
unseeded, so no reproducibility contract exists. The final phase is intentionally
Gauss-Seidel — the asymmetry is correct.
**Step-zero verification test** (round-1 redesign): all-gaussian-predictor fixture
(the r == 2 warm-start is then deterministic means, `R/prep_functions.R:203-209`;
factor warm-starts consume RNG in visit order and would confound the permutation)
with mocked deterministic fits; assert `pred_missing_run` byte-identical under
permuted formula order. New `n_cores_conv = 1` default; PR wall-time report states
the realistic ceiling (Σtᵢ / max tᵢ). Separate PR.

**4d. Per-run ESS floor (benchmark-gated; may not land).**
Round-1 reframe stands: pooled-stacked ESS was the wrong target (MC error =
within-run chain noise + B/M, and B/M is reducible only via `n_final`); round 2
adds that per-cell interval calibration is *insensitive* to per-run ESS (each
dataset contributes one draw per cell) — the floor's real consumers are per-run
diagnostics and `prob_preds` stability. Opt-in **`target_ess_per_run`** (default
NULL): estimate per-variable ESS/iteration from the convergence phase's last-run
models and resize `nitt_list` to the floor. Round-2 safeguards: **up-sizing only
by default** (never shrink below the user's explicit/default nitt; an
`allow_shrink` escape may come later with its own gate); variables **without a
last-run model keep their nitt with a message** (OVR path stores NULL,
`R/bace_imp.R:308, 346-348`; a fully slim object — the pattern the plan itself
endorses — warns and makes the option a no-op). **Floor default 200, correctly
anchored** (round-2 citation check): the single-chain ESS ≥ 200 convention is
Tracer's (Rambaut et al. 2018, *Syst Biol* 67:901-904) and matches the package's
own `.check_final_diagnostics` threshold (`R/model_functions.R:1638-1640, 1731`);
Vehtari et al. 2021 (*Bayesian Analysis* 16:667) is cited only for the modern
multi-chain recommendation (rank-normalised R-hat < 1.01, ESS > 400 across ≥4
chains — it never states 200), with a note that MCMCglmm runs one chain per fit.
Resolved lists go into the checkpoint manifest. Gate (dev/ study): pooled
estimates and CI widths within MC error of fixed-length reference; per-cell
coverage unchanged; zero new diagnostic flags; `prob_preds` stability; **realized
final-phase ESS reported against the convergence-phase prediction** (the
transferability assumption is itself checked); wall time reported. If the gate
fails, drop the feature and keep the study as the documented reason.

**Not a lever: `n_final`** (`R/bace.R:24-33`).

---

## Cross-feature interactions

- **1 + 3:** the fingerprint includes the tree list, `tree_index` and per-tree
  hashes; each cached run records its tree index.
- **1 + 4a:** Ainv cache keyed per tree hash; validator pre-warms it pre-fork.
- **1 + pooling consumers:** `pool_posteriors` consumes only `$all_models`;
  `with_imputations`/`pool_mi` only `$all_datasets` (tree-agnostic). OVR stays
  invisible to `pool_posteriors` by design.
- **2 + 3:** `bace_validate` checkpoints into **per-rep subdirectories** (a shared
  directory would be refused — each rep's data legitimately differs).
- **4b/4d + 3:** chain lengths enter the fingerprint after resolution.
- **2 + everything else:** inherits multiPhylo and the final-phase arguments via
  explicit mirrored arguments.

## Code-hygiene fixes to make in passing (round-2 citation audit)

- `R/phylo_signal_summary.R:689` miscites the rescaling constant: "Amemiya (1981)
  Econometrica 49:1483" → the paper is *Journal of Economic Literature*
  19(4):1483-1536, and the constant `c = 16·sqrt(3)/(15·π)` is better attributed
  to the logistic-probit scaling literature (Johnson & Kotz 1970; Zeger, Liang &
  Albert 1988 *Biometrics* 44:1049-1060; cf. Hadfield's CourseNotes), keeping
  Amemiya as a cf. Fix before feature 2's roxygen copies it.
- When feature 1 lands in roxygen/ms, add the verified references to
  `bib/BACE.bib` (Nakagawa & de Villemereuil 2019 — correct `nakagawa2018`'s year;
  Huelsenbeck 2000; Housworth & Martins 2001; de Villemereuil 2012; Rangel 2015;
  mulTree; Honaker 2011; Vehtari 2021; Rambaut 2018; Shen 2000; Rubin 2003).

## What this does NOT change

- The two-phase contract (`sample = FALSE` convergence / `sample = TRUE` final).
- Default MCMC settings, priors, `n_final = 50`, `sweeps = 3`, `.model_fit`'s
  signature.
- Object shapes — only additive fields.
- CRAN v0.1.0 submission content (on hold awaiting Dan's go). All four features
  are 0.2.0 material unless Dan decides otherwise.

## Effort estimate (coding + tests only; gate studies listed separately)

| Feature | Coding + tests | Gate study | Risk |
|---|---|---|---|
| 4a inverseA cache | ~½–1 day | — | low |
| 1 multiPhylo | ~3 days (pooling alignment is most of it) | ~1–2 days (replicated gold arm, ~100 seeds) | medium |
| 3 checkpoint/resume | ~2–3 days | — (write-path benchmark in PR) | medium |
| 4b final_* args | ~½ day | — | low |
| 2 bace_validate | ~2–3 days | ~½ day (dev metric gates) | low |
| 4c conv parallelism | ~1–2 days | — (step-zero test is the gate) | medium |
| 4d per-run ESS floor | ~2 days | ~1 day (calibration study) | high (gated) |
