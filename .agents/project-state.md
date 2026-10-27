# Project state (updated 2026-10-26)

Keep this current — it's the shared "where are we" so nobody re-derives it.

## Health
- `devtools::test()` = **1418 pass / 0 fail / 0 skip** (~1.4 min; R 4.4.2,
  MCMCglmm 2.36). (One benign, stochastic warning from the intentional
  non-convergence fixture in test-bace_wrapper.R.)

## 2026-10: Ayumi's review → five pre-release fixes (all landed on main)
Ayumi's PR #15 review raised four implementation issues; verification
confirmed **all four** and found a fifth. Full write-up:
[investigation-2026-10.md](investigation-2026-10.md); plan + status:
[plans/fixes-2026-10-ayumi.md](plans/fixes-2026-10-ayumi.md). Headlines:
1. (worst, new) Default categorical predictions silently dropped ALL
   phylogenetic effects (BLUP name mismatch) — the 2026-08 categorical sim
   numbers were stale by construction. Regenerated 2026-10-27 (1000
   reps/cell, same seeds): accuracy 0.799/0.815 (MCAR/MAR), balanced
   accuracy 0.651/0.642 — ~90% of the 0.722/0.709 Bayes ceilings (from
   ~88%) — prob-set coverage 0.992/0.993, runtime halved (22→12.5 s/rep).
2. Final phase understated imputation uncertainty when ≥2 variables missing
   (one pass from a fixed start; first-listed variable saw identical
   predictor fills in every run) → new `sweeps` (default 3) + `seed`
   arguments on `bace_final_imp()`/`bace()`.
3. `prob_preds` were overconfident plug-ins → posterior-averaged.
4. Diagnostics never covered cutpoints or final-phase models → they do now
   (`bace_final$diagnostics`).
5. Priors now user-settable per variable (`priors=`); defaults unchanged.
Sweeps calibration (dev/15) DONE 2026-10-27: **default sweeps = 3 confirmed**
(matches the 8-pass reference on every metric at ~half the cost); the 1-pass
arm quantified the old bug — slope biased up ~15% of truth via incoherent
y–x pairs, between-imputation spread understated ~13–20% (full table in
investigation-2026-10.md). Study A re-run DONE 2026-10-27: numerically
identical to the 2026-08 run (max |Δ slope| ≈ 1e-8; 240/240 ok) — no
regression, auto-collapse confirmed on the y-only path. R CMD check --as-cran
PASSED 2026-10-27 (only the "New submission" NOTE + coveralls-badge 403 false
positive). dev/13 categorical cells DONE 2026-10-27 (numbers above).
**Whole 2026-10 validation campaign COMPLETE (2026-10-27):** calibration,
Study A, categorical cells, SIMULATION_REPORT, R CMD check, tutorial
re-render + publish (render now ~7 min — vectorised prediction removed the
old 1.5 h bottleneck; bib restored after reference-manager overwrites, see
95997ef), tarball rebuilt + uploaded to win-builder R-devel (results email
pending). CRAN submission remains ON HOLD until Dan posts the PR #15 reply
and gives the go.

**2026-10-27, Shinichi's Rubin study imported:** `dev/rubin_study_2026-09-26/`
— BACE vs frequentist MI, 18 scenarios (n 100/300/1000 × λ 0.3/0.7/1 ×
ρ 0/0.5, 200 datasets/cell, M=20), continuous + discrete (threshold) traits;
entry point `report/study.html`; code/records/pages included; raw 1.4 GB rds
pools in a separate archive (not imported). Three load-bearing facts:
(1) the study **independently found problem 1** (2026-09, a month before
PR #15): "BACE as shipped was not proper" — all 20 final imputations shared
one posterior-mean-filled start, correlation coverage collapsed at n=1000;
(2) its patched arm ("bace_chain", sequential final imputations) restored
coverage — and its discussion proposes exactly the independent-runs/sweeps
design BACE now ships; (3) **the BACE arms predate the 2026-10 fixes**, so
before these numbers support the ms, the BACE arm should be re-run on the
fixed package (HPC launchers are in `code/`). NOTE: bundle README says
"Private; not for circulation" — import committed locally but NOT pushed
until Dan confirms with Shinichi (BACE repo is public).

**2026-10-27, Rubin-study harness converted for the re-run (in place, per
Dan):** arm `bace` = **"Current BACE"** (corrected package, plain
`BACE::bace()`); "as shipped" + `bace_resid` retired; `bace_chain` kept,
pinned to `sweeps = 1`/link; display labels match the ms software table
(Rphylopars / castor). New local launcher
`code/rubin_campaign_local.sh` (tiered, resumable); `script -> code` symlink
restores his path convention; Rphylopars 0.3.10 + castor 1.8.7 installed
(his exact versions). Package hardened (a80586d): warm-up sweep fit failures
skip that variable's update instead of killing the cell (surfaced by his
λ=1 fixed-thresholds fixture; lambda-1 cells would otherwise fail at
elevated rates). Validated: smoke cell green (all arms, zero errors),
BACE suite 1422/0, HIS suite 398/0 on the corrected package. Morning's
win-builder upload predates a80586d — rebuild tarball + re-upload once the
day's package changes settle.

**2026-10-27, campaign scoped to BACE-only and LAUNCHED (Dan's call):** goal
is accuracy / low bias / nominal coverage for Current BACE across the full
grid, not a competitor race (freq arms remain one `ARMS=` env away; the
castor Table-1-row question is moot). Pilot on this Mac: 72 cells = all 18
scenarios × 5/5/2 seeds (n=100/300/1000), CORES=7, output
`dev/rubin_study_2026-09-26/results_2026-10/` (untracked until the results
policy is decided). **n=100 tier early read (30/30 cells, zero errors, zero
fit failures):** per-value 95% coverage 0.93–0.97 for every trait × λ;
zRMSE falls with λ exactly as it should (c1: 0.88 → 0.66 → 0.19);
`frac_B0 = 0` everywhere (the improper-MI signature is gone); discrete
accuracy bin/cat3/ord = 0.77/0.71/0.66 with credible-set coverage ≈ nominal;
pooled slope/cor 95% intervals covered 30/30. **Watch-item (verified
pre-existing, NOT from our fixes):** at ρ=0.5 point estimates attenuate
toward zero as λ grows (slope bias ≈ −0.05/−0.05/−0.13 at λ=0.3/0.7/1,
5 seeds). Shinichi's 2026-09 run shows the same pattern with ~40× the seeds
(−0.016/−0.019/−0.074), his report names it ("chained BACE pulled the slope
toward zero at n ≤ 300", fading at n=1000), and joint-normal Rphylopars sits
at −0.02 at λ=1 — i.e. a known cost of chained univariate phylo models with
30% missing in *all* traits at complete signal; conservative direction
(attenuation + over-wide intervals), needs an honest ms sentence, not a code
fix. **Bonus finding:** his 2026-09 campaign silently dropped ~16% of λ=1
cells (198/200 recorded failures at λ=1, dominated by "Mixed model equations
singular") — the hardened package completed all 11 λ=1 pilot cells (10 + a
bonus seed 6), removing a selection bias his λ=1 numbers carried.

**2026-10-27, "50+ warnings per cell" census → benign, and fixed (65e6fa7):**
a warn=1 rerun of a warner cell showed every warning is ONE message — the
pseudo-Gelman prior announcing itself — fired on every categorical prior
build because the default `BACE.gelman = 2` satisfied `gelman == 2` in the
warning condition: 20 final runs × 3 sweeps = exactly 60 per cell. Clean
cells are the λ=1 cells where fixed thresholds leave cat3's third class
unrealised, the harness's droplevels makes cat3 *binary*, and no categorical
model exists (also why λ=1 "cat3" results are really a binary trait — his
documented DGP degeneracy, not a bug). Fix: warning now fires only for the
unrequested fallback (gelman = 1, complete cases < J); `bace_final_imp` no
longer builds-and-discards a categorical prior under `ovr_categorical =
TRUE` (now matches `bace_imp`; prior values/RNG unchanged); suite 1424/0.
**Operational: campaign workers run the INSTALLED package — do NOT
reinstall BACE while cells are running; reinstall after the pilot, before
scale-up** (numerically identical, just quiet logs). Census bonus facts:
harness runs `skip_conv = TRUE` with a fixed 5-iteration chain (Shinichi's
2026-09 budget design — `converged: FALSE, drift: 1e+07` in diag is that
sentinel, not a failure), and `ovr_categorical = TRUE` (his convention,
pre-dating the multinomial default switch).

**2026-10-28, PILOT COMPLETE (73/73 cells, zero errors) → scale-up
launched.** Full-pilot verdict: per-value 95% coverage 0.946–0.961 at every
n × trait; zRMSE improves with n and λ; `frac_B0 = 0` everywhere; discrete
accuracy rises with n (bin 0.78→0.82, cat3 0.72→0.75). Slope attenuation at
ρ=0.5 shrinks with n (λ=1: −0.11 → −0.05 → −0.05 at n=100/300/1000) but
intervals shrink faster (width 1.32 → 0.14), so ALL 6 coverage misses (of
146 estimand cells) are hairline low-side misses at ρ=0.5, 4 of them at
n=1000 — the exact regime Shinichi's report flags ("except complete
phylogenetic signal" at n=1000); ρ=0 cells are spotless (no invented
signal) and the complete-data reference missed 0/73 (perfectly calibrated
DGP+model). Flat tables: `results_2026-10/pilot_{estimands,cells,disc_cells}.csv`.
**Scale-up:** fixed package reinstalled (65e6fa7 live), full campaign
relaunched 2026-10-28 at launcher defaults (50/50/25 seeds = 750 cells;
73 pilot cells resume-skipped; ~1,050 core-h ≈ 6.3 days at CORES=7;
detached `nohup … caffeinate -i`, PID in
`results_2026-10/campaign_launch_2026-10-28.log`; stop any time — fully
resumable). Next after n=100 tier: interim read; after completion:
aggregate + build the BACE-only accuracy/bias/coverage figures into
`ms/ms.qmd` (Dan: "we build these figures into the paper with just BACE").
Win-builder: 2026-10-28 local R CMD check on post-fix source **0e/0w/0n**;
fresh tarball uploaded to win-builder R-devel (supersedes the 2026-10-27
round — Dan's results email will reference this build; CRAN still ON HOLD
until the PR #15 reply is posted and Dan gives the go).

**2026-10-27, software-comparison table (Track C):** `ms/ms.qmd` "Gap in
current tooling" now carries a 16-row feature table (`#tbl-software`) —
phylogeny / proper MI / multilevel / five variable types — replacing the
floating `ms/MI_table.pdf` (Shinichi's draft; now superseded, file left in
place). Every cell verified against the package's own reference manual
(one reader agent per package + independent spot-checks of the decisive
pages); page-cited evidence: [software-comparison-evidence.md](software-comparison-evidence.md).
13 new Crossref-verified bib entries. Headlines: no other tool combines a
tree with proper MI; mice has NO count method; miceadds multilevel cannot
impute nominal; only smcfcs (single-level), JointAI (no phylogeny) and hmi
(removed from CRAN 2022) cover all five types.
- **2026-08: six root causes of the Study B performance failures found, fixed,
  and regression-tested.** Full write-up: [investigation-2026-08.md](investigation-2026-08.md).
  Headline: `bace_final_imp` never reinstated response NAs, so every final
  "posterior-predictive" model conditioned on the convergence chain's point
  imputations as observed data (improper MI → the gaussian 0.76 / poisson 0.69
  under-coverage); plus the poisson mean-of-exponentials point estimate &
  fixed 1e6 rate clip (the catastrophic counts), gaussian-copied poisson
  priors, R nu=2 (the 2026-04 "fix" — actually inflates residual variance and
  biases λ̂ down; reverted to nu=0.002), a degenerate multinomial sim DGP
  (Bayes ceiling ~0.50 — categorical was near-optimal all along), and a
  marginal-AIC poisson/gaussian type-detection race that silently imputed
  strong-signal count variables as gaussian.

## Done this cycle (2026-08)
- **Package fixes** (all with regression tests in `test-performance-fixes.R`):
  NA reinstatement in `.one_final_imp` runs; poisson `post_median` point
  estimate + data-adaptive rate guard (3·max_obs+5); family-specific priors
  (gaussian R V=1 nu=0.002; poisson R V=0.1 nu=1, PX alpha.V=1);
  `ovr_categorical` default → FALSE (true multinomial: better calibrated,
  ~2.6× faster); OVR models kept in `all_ovr_models`; deterministic
  non-negative-integer → poisson type rule.
- **Simulation machinery**: `sim_bace` non-degenerate multinomial DGP
  (`beta_resp_liab` per-liability slopes, independent per-liability random
  draws, returns `response_probs` for Bayes-ceiling reporting); dev/13 v3
  metrics separating model calibration from finite-draw artifacts
  (`cover95_t`, prob-based discrete sets, Spearman/log recovery,
  per-rep `bayes_bal_acc`).
- **Spot-checked on the exact worst pre-fix replicates**: poisson rep 587
  max imputation 1,039,177 → 17; poisson MAR rep 1 coverage 0.464 → 0.964;
  gaussian MAR rep 1 coverage → 0.895(q)/1.00(t).

## Re-validation (2026-08-09, full re-runs complete)
Study B = 9,997/10,000 replicates ok, 5.3 h wallclock. Same seeds/config as the
2026-07 run, so differences are pure fix effect. Report re-rendered
(`simulation-report.html`); pre-fix archive in
`dev/simulation_results/archive_2026-07_prebugfix/`.

| Type | Recovery pre→post | Cover(q) pre→post (ceiling 0.845) | Cover t/prob (target .95) |
|---|---|---|---|
| gaussian | 0.83/0.88 → same | 0.756/0.784 → **0.844/0.845 (= at ceiling)** | **t: 0.949/0.947** |
| poisson | 0.58/0.62 → **0.69/0.74** (log-scale 0.75/0.78) | 0.691/0.880 → 0.867/0.902 | t: 0.916/0.932 |
| binary | 0.86 → same | 0.944 → 0.98 (draw-set, inflated) | prob: 0.941/0.938 |
| ordinal | 0.88 → same | 0.968 → 0.99 | prob: 0.959/0.957 |
| categorical | acc 0.50/0.52 → **0.81/0.78**; bal-acc 0.63/0.65 vs Bayes ceiling 0.71/0.72 | 0.751/0.769 → 0.984 | prob: 0.984 |

Poisson catastrophes: P(|std bias|>1) 14.6%→6.5% (MAR), 8.8%→2.8% (MCAR); max
|std bias| **217,838 → 5.3**. Categorical runtime ~30-42s → ~22s (true
multinomial default). Study A: MAR slope coverage 0.925→**0.950**, MCAR 0.975;
complete-case MAR mean-bias still corrected (−0.259→+0.011); CI width 1.59→1.83
(restored between-imputation variance); cell correlation unchanged.

**Remaining honest limitations** (for the ms): poisson MAR t-interval coverage
0.916 (log-link extrapolation on MAR-hidden high-leverage cells is
intrinsically hard; mild positive mean bias ~0.15-0.18 SD, median ~0.04-0.08);
categorical prediction sets slightly conservative (0.984); 3/10,000 replicate
failures (2 length-zero-argument errors, 1 singular mixed-model equations) —
low-priority follow-up.

## Tutorial / vignette modernisation (2026-09-18)
- **New workflow figure**: `vignettes/figs/bace_function_map.{svg,png}` — a
  bayesTLS-Figure-S1-style function map (stage chips, solid main pipeline,
  dashed optional/advanced routes, `bace()` wrapper band with the auto-retry
  loop, and the `with_imputations()` → `pool_mi()` Rubin branch). Edit the SVG,
  regenerate the PNG with `rsvg::rsvg_png(svg, png, width = 4200)`.
- **`vignettes/bace.qmd` updated for the post-2026-08 API**: `sim_bace()` now
  returns the species column as **`Species`** (capital S) — every chunk's
  `ran_phylo_form`/`species_col` was updated (the old lowercase forms all
  errored); deterministic type-detection table added to the types callout
  (integer ≥ 0 → poisson; store integer-coded continuous as doubles);
  `ovr_categorical = FALSE` default documented; `n_final = 50` default
  documented (executed examples pin `n_final = 10` for build speed, with a
  note); new top-level section *Downstream analysis with Rubin's rules*
  (`with_imputations()` + `pool_mi()`, smoke-tested end-to-end with a poisson
  GLM); the `species = TRUE` no-replicates error demo now actually errors
  (`error: true` chunk); outcome-dependent convergence prose made robust.
- **`bib/BACE.bib`**: added 8 entries the vignette cited but the bib lacked
  (pagel1999, blomberg2003, fritz2010, zhou2010, gelman2013, rubin1987,
  barnard1999, vanbuuren2018). All keys verified to resolve via a dummy render.
- **Rendered and published (2026-09-18)**: full render completed against the
  current package (all executed chunks ok, figure + Rubin section + error demo
  verified in the output); `vignettes/index.html` copied to root `index.html`
  and pushed (commit 7502548).
- **CI incident + fix (2026-09-18)**: the tutorial commits initially broke
  `R-CMD-check` — new files under `vignettes/` (figs/, an accidentally
  committed `bace.rmarkdown` Quarto intermediate) leaked into the package
  tarball because `.Rbuildignore` only excluded `vignettes/*.qmd`, triggering
  the fatal-under-`error-on: warning` "no VignetteBuilder" WARNING. Fixed
  (commit 4aec764): `.Rbuildignore` now has `^vignettes$` (tutorial is a
  Pages site, not a package vignette), the intermediate is untracked, and
  `vignettes/.gitignore` ignores `*.rmarkdown` / `/index.html` /
  `index_files/`. Verified by tarball inspection; CI green again.
- **Local build gotcha**: `/usr/local/bin/pandoc` on this machine is a broken
  x86 binary ("Bad CPU type"); render the vignette with quarto's bundled ARM
  pandoc: `export RSTUDIO_PANDOC=/Applications/quarto/bin/tools/aarch64` (and/or
  prepend that dir to PATH) before `quarto render vignettes/bace.qmd`. The
  rendered self-contained `vignettes/index.html` is copied to the repo root
  `index.html` (the GitHub-Pages tutorial). A full render takes ~1.5 h of MCMC.

## Open / priority (details in roadmap.md)
- **Track B** — production reference benchmark + competitor arm
  (missForest / Rphylopars). **Track C** — manuscript (restore missing bib
  assets; draft Methods/Results). **Release** — NEWS.md, trim exported internals.
- Optional follow-ups from the investigation: multiple stochastic sweeps per
  final run for multi-variable-missing data; strengthen the convergence
  diagnostic (Geweke silently disabled); fit-once-draw-many optimization when
  only one variable is missing.
