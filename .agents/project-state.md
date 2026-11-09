# Project state (updated 2026-11-09)

Keep this current — it's the shared "where are we" so nobody re-derives it.

## 2026-11-09: simulation campaign rebuilt for the new design

Plan for co-authors: [plans/simulation-plan-2026-11.md](plans/simulation-plan-2026-11.md)
— **awaiting Shinichi's approval**; it has four open decisions in its §8.
Code: `dev/archive-2026/rubin_study_2026-09-26/code/` (start at its `README.md`, or run
`Rscript code/walkthrough.R` for one dataset with every step printed).
Lane suite **455 pass / 0 fail**.

Design now: n ∈ {100, 500} × λ ∈ {0.3, 0.7, 0.9} × three information regimes
(ρ_traits 0 / 0.5 / 0.5-plus-a-covariate-at-0.8) × {mcar, mar, clade},
1000 reps = 54,000 datasets ≈ 72,000 core-hours ≈ 16 days on 200 cores.

**Two silent failures found and fixed — worth remembering, both cost nothing
to make and would have cost weeks:**
1. **The installed BACE was not the repo's.** Repo at 1.0.0, R library still
   at 0.1.0.9000 (built 2026-11-07). Everything ran on the old package and
   looked entirely normal. `rubin_cell.R` and the launcher now refuse below
   `BACE_MIN_VERSION`, and every result file records `bace_version`.
2. **Two scenarios wrote the same filename.** The runner skips a dataset whose
   `.rds` exists (that is what makes long runs resumable), so the "typical"
   and "good covariate" regimes — identical but for the covariate — would have
   collided and a third of the design would have been skipped silently.
   Filenames now come from one function, `cell_tag()` in `rubin_lib.R`, with a
   test over all 54 scenarios.

**Three measured facts the manuscript should use:**
- Complete-case survival: **MCAR 8%** (= 0.7⁷, so it fails by *attrition*, not
  bias), **MAR ~14%**, **clade ~68%** (so clade is where the *bias* story is).
- The realised covariate correlation scatters hugely — requested 0.8 gives
  mean 0.799 but sd 0.096 over 300 datasets at n = 120; requested 0.3 gives
  sd 0.222. Report it as a distribution, not as "0.8".
- A single auxiliary at 0.8 with *all seven* traits is not positive definite:
  the bound is r² < (1 + (K−1)ρ)/K, i.e. r < 0.756 at K = 7, ρ = 0.5.

**Open, for Shinichi:** at n = 100, ~23% of the fixed effects in the *final*
phase (the one that feeds inference) have effective sample size < 100. BACE
1.0.0's new `final_nitt`/`final_thin` can address it. A ~20 core-hour pilot
settles it against a 72,000 core-hour campaign — recommended before launch.

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
Sweeps calibration (dev/archive-2026/15) DONE 2026-10-27: **default sweeps = 3 confirmed**
(matches the 8-pass reference on every metric at ~half the cost); the 1-pass
arm quantified the old bug — slope biased up ~15% of truth via incoherent
y–x pairs, between-imputation spread understated ~13–20% (full table in
investigation-2026-10.md). Study A re-run DONE 2026-10-27: numerically
identical to the 2026-08 run (max |Δ slope| ≈ 1e-8; 240/240 ok) — no
regression, auto-collapse confirmed on the y-only path. R CMD check --as-cran
PASSED 2026-10-27 (only the "New submission" NOTE + coveralls-badge 403 false
positive). dev/archive-2026/13 categorical cells DONE 2026-10-27 (numbers above).
**Whole 2026-10 validation campaign COMPLETE (2026-10-27):** calibration,
Study A, categorical cells, SIMULATION_REPORT, R CMD check, tutorial
re-render + publish (render now ~7 min — vectorised prediction removed the
old 1.5 h bottleneck; bib restored after reference-manager overwrites, see
95997ef), tarball rebuilt + uploaded to win-builder R-devel (results email
pending). CRAN submission remains ON HOLD until Dan posts the PR #15 reply
and gives the go.

**2026-10-27, Shinichi's Rubin study imported:** `dev/archive-2026/rubin_study_2026-09-26/`
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
"Private; not for circulation" — the import was held locally until
2026-10-28, when Dan explicitly authorised pushing the full history to the
public repo (confirmed via a direct question naming the bundle; pushed at
f462ce3). **Standing push policy since then: a 3-hourly checkpoint commits
and pushes SIMULATION FILES ONLY (results_2026-10 rds + campaign logs,
enforced by pathspec); everything else is committed/pushed only when a
human says so.**

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
win-builder upload predates a80586d. (A re-upload followed on 2026-10-28,
but it was not authorised — see the standing rule below.)

**2026-10-27, campaign scoped to BACE-only and LAUNCHED (Dan's call):** goal
is accuracy / low bias / nominal coverage for Current BACE across the full
grid, not a competitor race (freq arms remain one `ARMS=` env away; the
castor Table-1-row question is moot). Pilot on this Mac: 72 cells = all 18
scenarios × 5/5/2 seeds (n=100/300/1000), CORES=7, output
`dev/archive-2026/rubin_study_2026-09-26/results_2026-10/` (git-tracked since dd05d7b,
2026-10-28 — Dan committed the pilot rds/logs/flat tables himself; commit new
cells with the same pattern as the campaign progresses). **n=100 tier early read (30/30 cells, zero errors, zero
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

**2026-10-29, n=100 tier COMPLETE (300/300 cells) — interim read:** headline:
**nominal coverage everywhere that was broken pre-fix.** Per-value 95% coverage,
pooled over datasets per scenario: continuous (c1+c2) 0.946–0.970, counts
0.948–0.964, proportions 0.951–0.969 (logit scale). Compare 2026-09 (same DGP,
200 seeds): as-shipped BACE 0.901, bace_chain 0.954, freqA 0.941 → Current BACE
0.955 matches the chained variant and beats the frequentist arm, with the
improper-MI signature gone. Downstream slope/cor coverage 0.96–1.00 per scenario
(slightly conservative; complete-data reference itself 0.90–0.98); overall bias
−0.02, almost all carried by the known λ=1 × ρ=0.5 attenuation cell (slope bias
−0.104 there, the pilot watch-item, needs the honest ms sentence). Discrete:
accuracy bin/cat3/ord 0.78/0.67/0.61 ≈ 2026-09 values (Bayes-ceiling-bound);
credible-set coverage 0.93–0.96 ≈ bace_chain, better calibrated than as-shipped;
cat3 NA-draw rate unchanged (~39/dataset vs ~35 in 2026-09 — known, not new).
**Correction to the pilot's λ=1 optimism:** the full tier has **15/100 λ=1
failures** (4 at ρ=0, 11 at ρ=0.5), all `MCMCglmm "Mixed model equations
singular: use a (stronger) prior"` — same error class as 2026-09 (which had
≥47/400 inline + 71 as-shipped-failed at λ=1 n=100), so the hardened package
fails less often and *records* the failure, but λ=1 is not failure-free as the
11-cell pilot suggested. Consequence: surviving λ=1 ρ=0.5 cells (39/50) are not
a random subsample — the same selection caveat we raised about Shinichi's λ=1
numbers applies to ours and belongs next to any λ=1 figure. Possible remedy if
Dan wants it: harness-level retry with a stronger prior via the new `priors`
argument (the exact remedy MCMCglmm's error message names). Campaign log scan:
only the benign pre-65e6fa7 "50+ warnings" lines + the 30 singular-MME lines;
nothing else. Interim script:
scratchpad `interim_n100.R` (session a86edbd4); n=300/1000 tiers still running.

**2026-11-01 update — singular-MME failures quantified across the whole campaign
(661/750 done):** 20 failed datasets total (3.0%), ALL λ=1, ALL the bace_fit arm,
ALL "Mixed model equations singular". Rate falls with n: n=100 15/100 (ρ=0: 4,
ρ=0.5: 11), n=300 5/100 (2, 3), n=1000 0 so far (6 datasets done). Every cell has
its full 50 (or 25) seed datasets — failures are recorded inside the rds
(`$errors$bace_fit`), no seeds were replaced. **Survivorship is real, measured:**
on the 204 finished λ=1 datasets, the complete-data reference fit shows failed
draws have larger realised slopes (0.37 vs 0.24, Wilcoxon p = 0.036) and larger
SEs (p = 0.036) than survivors — failures hit the strong-signal draws, so λ=1 BACE
cell summaries condition on a slightly weaker-signal subset, and the −0.104
attenuation estimate at λ=1 × ρ=0.5 is if anything optimistic. Scripts:
scratchpad `tally_failures.R`, `survivorship_check.R` (session a86edbd4). For the
ms: report per-cell n_ok + one-sentence caveat next to any λ=1 figure. Note
`bace()` auto-retry does NOT help here (retries add runs, not stronger priors);
roadmap candidate: auto-strengthen the prior on singular-MME retry.

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
a fresh tarball went to win-builder R-devel the same day, but that upload
was **not authorised** (Dan, 2026-10-28: "You're sending without
authorisation to win-builder. Stop doing that."). **Standing rule: no
win-builder upload — or any other external submission — without Dan's
explicit go-ahead at the time. The next round happens only when BACE is
certainly fully fixed, just before CRAN submission.** The one results email
from the 2026-10-28 upload can be ignored. CRAN still ON HOLD until the
PR #15 reply is posted and Dan gives the go.

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
  draws, returns `response_probs` for Bayes-ceiling reporting); dev/archive-2026/13 v3
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
`dev/archive-2026/simulation_results/archive_2026-07_prebugfix/`.

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

## 2026-11-04: CAMPAIGN COMPLETE (750/750) → aggregated → ms figures built
- Full campaign finished 2026-11-04 (Dan committed the last cells, 659015a).
  **Final failure ledger: 22/750 datasets (2.9%), all λ=1, all bace_fit, all
  "Mixed model equations singular"** — per-cell n_ok: ρ=0 46/48/24, ρ=0.5
  39/47/24 at n=100/300/1000 (of 50/50/25). The survivorship caveat from
  2026-11-01 stands (failures hit strong-signal draws; λ=1 summaries slightly
  optimistic).
- **Aggregation**: new `code/rubin_full_aggregate.R` (flattens all 750 rds →
  `results_2026-10/{full,agg}_{cells,disc_cells/disc,estimands/down}.csv`,
  `agg_failures.csv`, `sanity_full.txt`). Headlines: per-value 95% coverage
  0.94–0.97 in every scenario for c1/c2/cnt/prp; discrete credible-set coverage
  0.89–1.00 (dips only ordinal λ=0.3 n=100); downstream (PGLS slope + phylo cor,
  Rubin-pooled M=20) coverage 0.94–1.00 everywhere, |bias| ≤ 0.04 except
  λ=1×ρ=0.5 slope attenuation −0.104 → −0.045 → −0.024 at n=100/300/1000
  (fmi 0.51–0.68 there, so coverage stays 0.96–0.98). Discrete traits scored at
  n ≤ 300 by launcher design. At λ=1 cat3 is effectively binary (DGP degeneracy).
- **ms/ms.qmd "Simulation results"**: placeholder replaced with a drafted
  3-paragraph results text (marked DRAFT for Shinichi & Dan) + two executable
  figure chunks reading the agg CSVs: `fig-sim-values` (per-value coverage /
  zRMSE / discrete accuracy, 7 traits) and `fig-sim-downstream` (bias +
  coverage for slope & correlation, the attenuation + honest-interval story).
  Captions carry the n_ok counts and the singular-failure caveat. Verified by a
  full `quarto render ms.qmd --to docx` (clean; only the pre-existing empty
  tbl-studyb warning). The old 2026-09 "beats frequentist arms" claim is left as
  a bracketed note — it comes from the paired 2026-09 study, not these runs.
- **System-library BACE reinstalled from source** post-campaign (2026-11-04) —
  the TODO below is done; interactive sessions now have the current printouts.
- 3-hourly checkpoint auto-commits retired (campaign over).

## Tutorial render gotcha: render uses the INSTALLED package (2026-11-01)
- The tutorial's setup chunk loads BACE with `pacman::p_load(BACE)` — i.e. the
  **installed** copy, not the source tree (the `install_github` line was removed
  deliberately so renders can't silently swap the package mid-campaign). After
  editing `R/`, the rendered HTML will NOT reflect the change until the package
  is reinstalled. This bit us: the 2026-10-31 printout updates (accessor-first
  messages, acronym fix, commit 8a8a6aa) were missing from two subsequent
  renders; Dan spotted the stale printouts in the published page.
- Correct procedure while a simulation campaign runs against the system library
  (do NOT overwrite it — live workers lazy-load from the installed rdb and the
  campaign must keep running the code it started with): install the source into
  a scratch library and point the render at it —
  `R CMD INSTALL --library=$LIB . && R_LIBS=$LIB quarto render vignettes/bace.qmd`.
  Then verify the output actually changed (e.g. grep the new printout text)
  before copying to root `index.html` and committing. Fixed render pushed as
  0733031. The render is fully seeded, so a re-render is byte-identical except
  for genuine changes — a small diff is expected, not suspicious.
- **TODO when the Rubin campaign ends**: `R CMD INSTALL .` into the system
  library so interactive sessions also get the current code.

## Open / priority (details in roadmap.md)
- **Track B** — production reference benchmark + competitor arm
  (missForest / Rphylopars). **Track C** — manuscript (restore missing bib
  assets; draft Methods/Results). **Release** — NEWS.md, trim exported internals.
- Optional follow-ups from the investigation: multiple stochastic sweeps per
  final run for multi-variable-missing data; strengthen the convergence
  diagnostic (Geweke silently disabled); fit-once-draw-many optimization when
  only one variable is missing.

## 2026-11-04: manuscript simulation sections tightened and checked

The source sections in `ms/ms.qmd` (original lines 216–325) were edited against
campaign code, aggregate CSVs and saved RDS files. This supersedes the earlier
claims that discrete coverage was uniformly nominal, the categorical trait was
always binary at λ=1, and failure selection necessarily made results optimistic.
Numeric prediction-interval coverage is 0.9402–0.9718; discrete-set coverage is
0.8861–0.9955 (lowest: ordinal, n=100, λ=0.3, ρ=0.5). Among the 250 λ=1 datasets,
`length(unique(as.character(x$truth$cat3)))` gives 198 with three realised classes,
43 with two, and 9 with one. All 22 failed fits were at λ=1; manuscript wording
now states that those summaries condition on successful fitting. Stronger priors
were suggested by the error message but were not tested as a remedy here.

Methods now distinguish five fixed initial passes from confirmed convergence:
`skip_conv=TRUE` disables retries, not assessment. The first two final-phase
passes use 11,250 iterations; the retained pass uses 50,000 (10,000 burn-in,
thinning 25). Scoring definitions now match `code/rubin_lib.R`,
`code/rubin_cell.R` and `code/rubin_discrete.R` under
`dev/archive-2026/rubin_study_2026-09-26/`. Both slope and correlation attenuate at λ=1, ρ=0.5.
No simulation results or package code were changed.

## 2026-11-05: performance-measure definitions clarified

`ms/ms.qmd` now gives explicit zRMSE and accuracy equations, defines coverage
separately, and documents fractional credit for tied modal categories. The figure
caption matches these definitions: zRMSE=1 means RMSE equals the observed SD;
fixed accuracy references represent uniform guessing, not majority-category
predictions. `code/oracle_floor.R` uses known parameters and observed-cell latent
values, averaging four trees and two masks per tree. Its references use population
scaling and are optimistic comparisons, not hard bounds on plotted dataset-level
summaries. No scoring code, simulation data or plotted values were changed.

## 2026-11-05: manuscript PDF and DOCX rerendered

Regenerated `ms/ms.pdf` and `ms/ms.docx` with the revised metrics text and native
equations. The two simulation plots now use `cairo_pdf` for LaTeX output, preserving
Greek plot labels. Portrait figures use a 6.9-inch DOCX output width to avoid
clipping inside Quarto figure tables. Full PDF and DOCX page previews inspected;
remaining Word-preview template quirks include duplicate footer page numbers and
a comparison-table row split across pages.

## 2026-11-05: bird case study and illustrated manuscript panels

Replaced the bird placeholder in `ms/ms.qmd` with a concise account grounded in
Ayumi's saved `vignettes/bird_social/results/` files. The example covers 10,563
species, with social records missing for 4,265; random and whole-family social
validation have category-averaged accuracies of 54% and 36%. Mixed-variable
validation of 2,113 measured hand-wing indices gives r=0.9646 and RMSE=3.9981.
Panel c shows wing-shape validation; panel d shows accuracy within each social
category, separating the two holdout sets. These use separate saved validation
fits; no imputation runs were repeated.

Panel d now includes four transparent scientific illustrations selected from
observed BIRDBASE records: Alcedo atthis, Tyto alba, Dacelo novaeguineae and
Cacatua galerita, in category order. Assets, source references and exact built-in
image-generation prompts are in `vignettes/figs/bird_social/`. They sit outside
the plotting area, are not to scale, and the caption identifies the coding
examples. Standalone exports are in `ms/figures/`; the PDF (18 pages) and Word
manuscript (17 pages) were regenerated, and the new figure/caption checked in
both formats. Concurrent user edits to the title and HWI intervals/wording were
preserved.

User instruction: edit project files directly; do not create scratch copies in
the repo. Temporary render copies from this work were removed. Keep future QA
previews and scratch work outside the repo.

## 2026-11-05: Stacked-posterior arm complete; discrete-coverage experiment launched

- **Stacked-posterior arm done** (580/580 datasets, 0 MCMCglmm failures;
  `results_2026-10/stack_down.csv`, `agg_stack.csv`, generator
  `code/rubin_stack_posthoc.R`). BACE-native workflow (with_imputations ->
  pool_posteriors -> get_pooled_model), analysis model MCMCglmm(c2 ~ c1, phylo).
  Paired vs Rubin on identical datasets: coverage 0.972 vs 0.974, point estimates
  differ 0.0025 on average, stacked intervals ~6% narrower (biggest gap at
  lambda = 1, where Rubin's small-M df correction is most conservative). The
  lambda=1 x rho=0.5 attenuation (-0.10 at n=100 -> -0.04 at n=300) is identical
  in both arms: it lives in the imputations, not the pooling.
- **Discrete set-coverage investigation** (Dan's "50k too short?" question):
  (1) set rule exonerated - perfectly calibrated M=20 imputer scores >= 0.95
  under the campaign's exact rule (oracle demo); (2) generic chain length
  disfavoured - well-covered lambda=1 cells mix WORSE (ess_med 224 vs 730) and
  continuous traits in the same runs are nominal; (3) remaining hypotheses =
  OvR construction (cat3) and cutpoint uncertainty (ord).
  `code/coverage_experiment.R` now running: multi (ovr_categorical=FALSE) and
  long3x (150k/30k/75) arms, n=100 lambda=0.3 rho={0,0.5}, seeds 1-25, paired
  to campaign datasets. ~14 h on 7 cores; cache results_2026-10/coverage_exp_cache/.
- **Oracle references in ms**: `code/oracle_floor.R` -> `oracle_floor.csv`
  (mask-aware conditional-variance floors for zRMSE, Bayes ceilings + chance for
  accuracy); fig-sim-values panels (b)/(c) now draw per-scenario dashed lines;
  0 violations across all 108 observed cells.
- Hourly stacking checkpoint cycle retired (purpose fulfilled).

## 2026-11-06: Coverage experiment decided it — long runs go multinomial

- **Results (99/100 runs; paired to campaign datasets, pooled over rho):**
  - `long3x` (OvR kept, 3x chain 150k/30k/75, same 1,600 retained): coverage
    unchanged for all three discrete traits (bin +0.004, cat3 -0.003, ord +0.005,
    all << 1 SE) despite better mixing (ess_med 932) and 3x cost (89 vs 32
    min/dataset). Dan's "50k too shallow" hypothesis is refuted by controlled test.
  - `multi` (ovr_categorical = FALSE, campaign MCMC): cat3 set coverage
    0.902 -> 0.923 (paired +0.021, SE 0.009), accuracy 0.51 -> 0.53; bin/ord
    unchanged (clean internal control — the flag doesn't touch them).
  - Residual: ord ~0.90 / cat3 ~0.92 at n=100, lambda=0.3 remains below 0.95 with
    set rule, chain length, and OvR now all ruled out/quantified. Best remaining
    candidate: honest model error of the chained linear-predictor conditionals at
    small n + weak signal. Binary/continuous nominal in same runs; lambda=1 covers 0.98.
- **Decision (Dan, 2026-11-06): all long-run simulations use the multinomial**
  (`ovr_categorical = FALSE`, the shipped package default since the 2026-08 fixes;
  note the flag names the ALTERNATIVE — TRUE = one-vs-rest, FALSE = multinomial).
  Harness updated for Shinichi: `campaign_gnn_off_lib.R` run_bace() now passes
  FALSE; `rubin_bace.R::fit_bace_mi()` default flipped to FALSE (rubin_cell.R
  inherits it); `mi_bace_chain()` follows the parent fit via `fb$ovr`. Package
  itself needed no change. **Never resume a pre-2026-11-06 results cache across
  this change** — the 2026-10 campaign results were produced with OvR TRUE.

## 2026-11-07: Feature 1 implemented — tree-per-imputation (multiPhylo)

Plan: `.agents/plans/next-features-2026-11.md` (built on a 4-reader code map,
adversarially reviewed twice; feature 1 + its 4a prerequisite now LANDED).

- **What shipped** (version bumped to 0.1.0.9000; full NEWS bullet written):
  `phylo` accepts phylo | multiPhylo | list-of-phylo in `bace()`, `bace_imp()`
  (uses trees[[1]] — the convergence diagnostic needs a fixed target) and
  `bace_final_imp()` (run m fit under tree_index[m], deterministic in the run
  index, zero RNG consumed). New `R/tree_functions.R`: per-tree inverseA cache
  (`.get_ainv`/`.compute_ainv`, cleared at each phase entry; key = content hash
  incl. node labels), `.validate_trees` (tip coverage; negative AND zero-length
  branch errors with the 1e-8 remedy — inverseA rejects zero edges
  unconditionally, so this breaks nothing that worked; runs inverseA per tree
  up front so ultrametricity failures surface named at entry; strips node
  labels; root-depth CV warning), `.select_tree_indices` (subsequence /
  balanced recycling), `.hash_obj` (serialize-v2-minus-header, R-version
  stable — feature 3 will reuse it). `bace_final` gains $tree_index,
  $tree_hashes, $tip_labels. `pool_posteriors` detects trees_vary from
  tree_index (NOT colnames — same-tip trees have identical name SETS in
  different orders), stacks fixed/VCV/CP/Liab as before, aligns tip BLUPs BY
  NAME, drops node BLUPs ("\.Node\d+$" with an exact-tip-set fallback),
  records BACE_pooling$trees_vary + n_node_blup_cols_dropped; single-tree
  path byte-identical. `with_imputations(tree=)` accepts a tree list (tree m
  to fit m, balanced recycling; consistency message against $tree_hashes).
- **Verified**: new test file tests/testthat/test-multiphylo.R (68
  expectations) incl. the load-bearing IDENTITY test — on one shared converged
  object with the same seed, multiPhylo(A,B) runs 1,3 are bitwise-identical to
  single-tree-A runs 1,3 and runs 2,4 to single-tree-B (proves the tree
  reaches each run AND no stray RNG), serial==parallel under multiPhylo, the
  skip_conv=TRUE forwarding spy, pooling alignment on permuted-name fakes with
  a dotted tip label. Full suite: **1492 pass / 0 fail** (was 1424).
  R CMD check --as-cran running at time of writing.
- **Function map updated** (vignettes/figs/bace_function_map.svg + PNG
  regenerated): input node "raw data + tree(s) … one phylo, or a multiPhylo
  posterior tree sample"; bace_final_imp chip "with a multiPhylo, run m is fit
  under tree m"; with_imputations note "tree = one phylogeny for every fit, or
  one per imputation". ms.qmd/vignette embed the PNG automatically.
- **Still owed for feature 1**: the dev/ statistical gate study (replicated
  gold-standard arm + one-sided B-deflation test + sweeps dose-response, ~100
  seeds — design pre-registered in the plan) before the ms claims tree
  uncertainty is propagated correctly at default sweeps.

## 2026-11-07 (later): multiPhylo gate study run — leakage REFUTED, plateau not yet shown

`dev/archive-2026/16_multiphylo_gate.R` -> `dev/archive-2026/simulation_results/multiphylo_gate/`
(committed 87089b0). 100/100 seeds, ~1 min/seed, pre-registered verdicts:

- **A1 (primary, the feared anti-conservative deflation): PASS by a wide
  margin** (p = 7e-14). The shipped shortcut's between-tree variance component
  is NOT deflated relative to a gold standard that re-runs the chained phase
  on every tree with replicates — it is ~2x LARGER (tau2 0.049 vs 0.023), and
  downstream slope B is 8.5% HIGHER (ratio 1.085): the shortcut errs
  CONSERVATIVE, never anti-conservative. Interpretation: the shared tree-1
  start is not pulling tree groups together; each tree's model responds to
  the shared conditioning values differently, adding (conservative)
  between-tree spread that gold's per-tree starts average away into its
  (subtracted) chained-noise term.
- **A2 (positive control): SENSITIVE** — sweeps = 1 is significantly deflated
  vs sweeps = 10 (p = 0.007), so the design can detect initialization leakage
  when it exists. The default sweeps = 3 is what kills it.
- **A3 (plateau, sweeps 3 ~ sweeps 10 within 10%): FAIL** — point estimate
  -6% of the sweeps-10 level, CI [-18%, +5%]: equivalence not demonstrated at
  100 seeds (direction: slightly LESS of the conservative inflation, i.e.
  toward gold, not away from it).
- **Overall registered verdict: FAIL (via A3 only).** Practical reading for
  the ms/docs: tree uncertainty IS propagated (if anything over-propagated vs
  the expensive gold standard); no evidence of the anti-conservative failure
  mode; the only open question is cosmetic (sweeps 3 vs 10 differ by ~6% of
  tau2 in the conservative-to-less-conservative direction). Options for Dan:
  (a) accept + document (recommended: the conservative direction is safe);
  (b) bump default sweeps for multiPhylo runs; (c) extend the study
  (resumable cache; seeds are cheap at 1 min each).

## 2026-11-07 (later still): Feature 2 implemented — bace_validate() + get_imputation_intervals()

Per plan feature 2 (Dan asked for plain-language docs, simple functions,
accuracy). New `R/bace_validate.R`:

- **bace_validate()**: hides `prop` of each variable's OBSERVED cells (never
  the cluster; never empties a factor class; <10 observed cells -> skip +
  warn), runs the full bace() pipeline on the masked data, scores against the
  hidden truth. Metrics match the ms EXACTLY: zRMSE = RMSE of across-M mean /
  SD of remaining observed (mean-filler baseline printed beside it); coverage
  = per-cell Rubin t-interval, with a simulated fair expectation for
  low-count poisson (<2 mean rate); accuracy/balanced accuracy vs
  commonest-class baseline; Brier from MI-averaged prob_preds; optional
  per-variable logit scoring (ms clipping 1e-4). reps>=2 gives
  between-rep SEs. Scoring lives in its own seam (.score_validation) so every
  formula is unit-tested on hand-computed answers. Honest-docs: MCAR caveat
  (cannot detect MNAR), Amelia overimpute() comparison (ours masks a share at
  once = harsher), cost note (retry ladder multiplies). seed fixes masks
  exactly; scores only approximately (chained phase unseeded by design —
  doc corrected after a test caught my over-claim).
- **get_imputation_intervals()**: exported per-cell Rubin t-intervals from
  $all_datasets (nominal for any n_final >= ~20); the n_final roxygen in
  bace()/bace_final_imp() now points here instead of teaching the ~92%
  empirical-quantile recipe.
- **Verified**: 89 new expectations (hand-checked zRMSE/coverage/Brier/
  balanced-accuracy/logit arithmetic incl. adversarial prob-row matching with
  replicated species and pre-existing NAs; calibrated-imputer coverage sim
  0.92-0.98 band; mask bookkeeping + RNG restore). Full suite **1581 pass /
  0 fail**; `R CMD check --as-cran` **0/0/0** (one roxygen-markdown bracket
  autolink fixed). dev/archive-2026/17_validate_metric_gates.R at realistic scale (60 sp,
  prop 0.3, 5 reps/type): **ALL PASS** (gaussian zRMSE well under baseline,
  coverage in band; categorical accuracy 0.57 > 0.53 baseline).
- Vignette benchmarking section rewritten around bace_validate(); function
  map gains a "check it works" box (SVG + PNG regenerated). NEWS updated.
- **Outstanding**: the 3-agent adversarial review of this code died on the
  account's 5-hour usage limit before producing findings; it re-runs after
  the 12:20pm reset (timer set). Treat feature 2 as implemented but not yet
  adversarially reviewed.

## 2026-11-07 (validate follow-up): bace_validate() is now fit-first; map updated

Dan asked why validation sat at the START of the workflow and whether it could
take a finished bace object. Answer (now in the docs): a finished run cannot
be scored directly (nothing held out), so the check always RE-RUNS on masked
data — but `bace_validate(fit)` now takes the `bace_complete` object and
reuses its data, formulas, tree(s), chained-iteration count and settings, so
nothing is re-typed. Implementation: bace_validate is an S3 generic
(.bace_complete / .character / .list methods over one workhorse);
bace_final_imp now stores `$settings` (ran_phylo_form as given, nitt/thin/
burnin, species, ovr_categorical, nitt_cat_mult, priors) and `$trees_used`
(additive fields — also future checkpoint-fingerprint inputs). Older objects
without settings get a clear error. Function map: bace_validate() box retitled
("check a finished run"; arrow now data -> box) and a third accessor box
get_imputation_intervals() added under get_imputed_data() (accessor column
compacted to fit; pool_mi untouched). Both figure captions (vignette + ms)
name the new functions. Vignette leads with bace_validate(fit). Tests: 1598
pass / 0 fail (settings-handover mock test + end-to-end fit-method test).
Adversarial review of feature 2 still pending the usage-limit reset.

Correction to the map note above (Dan, same day): bace_validate() is now its
OWN dashed box in the bottom helper row, downstream of the bace() band (dashed
arrow from the band's bottom edge), not attached to the raw-data node; the
left node is back to "simulate test data" (sim_bace/sim_tree). PNG regenerated
and visually checked.

## 2026-11-07 (validate review): adversarial review done, all findings fixed

3-agent review (statistics / break-it-by-running-it / docs accuracy) of
bace_validate + get_imputation_intervals. Verified sound: zRMSE, t-interval,
logit path, count floors all reproduced exactly against ms + rubin_lib; the
fit-first handover faithful (mock + real fits); Brier row-matching safe.
Fixed from findings (all reproduced by the reviewers before fixing):
- MUST: list-formula fits with a complete predictor-only covariate crashed
  bace_validate(fit) — masking now targets only response (imputed) variables.
- MUST: roxygen \% under markdown double-escaped and CORRUPTED the rendered
  help pages (bace_validate, get_imputation_intervals, bace_final_imp) — all
  now "percent"/plain text (house style).
- MUST: bace_final_imp's new $settings/$trees_used return fields were
  undocumented — @return now lists them (+ miss_dat/n_final/sweeps/seed).
- Tie votes now get FRACTIONAL credit (ms rule; was first-alphabetical —
  ~11% of maximally-uncertain binary cells were decided by label order at
  even M); baseline uses the same rule.
- Count coverage expectation: always simulated from the variable's own rates
  (the mean<2 gate misreported mixtures on both sides), fixed internal seed
  (printed expectation is reproducible), n_sim scaled to ~20k total draws.
- Date-style seeds (e.g. 20261107) no longer overflow (mask seeds derived
  via RNG, not seed*1000); character vectors of formulas refused with a
  clear message; a formula method added; scale typos warn; settings
  overrides on the fit method warn instead of being silently dropped; type
  flips caused by masking warn; numeric masking always leaves >= 2 observed
  cells; summaries drop NA reps per metric; intervals degrade clearly at
  M=1 / all-factor; print wording fixes (SE label, zero-width instances,
  Brier 0-2 anchor, commonest OBSERVED class, per-repetition).
- NEWS leads with the fit-first form and scopes the "manuscript metrics"
  claim honestly; vignette's "well calibrated" quantile-interval claim
  corrected. Suite: **1626 pass / 0 fail**; R CMD check **0/0/0**.

## 2026-11-08: New functions integrated into both case studies

Dan asked to integrate bace_validate()/get_imputation_intervals() into
Patrice's and Ayumi's case studies, with two rules: outputs must match what
the case studies already do, and the ms case-studies figure must reproduce.
Mapping (2-agent workflow, archives str()-verified) + integration done:

- **What CANNOT be replaced, and why (documented in each case study):**
  both validation DESIGNS are structured — Patrice re-runs Pottier et al.
  2025's five fixed, species-clustered, condition-stratified folds; Ayumi
  holds out 630 random species + 631 species in 38 whole families (and
  measured-only hwi cells). bace_validate() masks random observed cells
  only (no fixed sets/folds/blocks; every response variable masked) and
  lacks their metrics (raw-unit RMSE, bias, RPS, recall/confusion,
  calibration, family bootstrap, bespoke baselines). The hand-coded designs
  stay authoritative; each case study now says so in one plain sentence.
  Also: NO archived object is a bace_final/bace_complete (plain summary
  lists only; results/tmp deleted; runs at dc79be1 predate $settings), so
  the functions cannot be CALLED on archives — the t-interval FORMULA
  integrates instead, recomputed from archived draws/means/sds, no refits.
- **What integrated (verified from archives before editing):**
  (1) Interval construction everywhere: amphibian CV + full-run intervals,
  bird hwi coverage, ms fig panel (a) switched from empirical quantiles to
  the Rubin t-interval (= get_imputation_intervals, proven identical to the
  ms panel-c hand formula to machine precision). Coverage: amphibian CV
  0.925 -> 0.941 (353/375; r/RMSE/bias unchanged 0.864/1.437/-0.451); bird
  hwi 0.856 (empirical 90%) -> 0.950 (exactly nominal) — the tutorial's
  "~82% expected even if calibrated" apology deleted. Widths ~+12%.
  (2) Displayed (eval:false) code now shows get_imputation_intervals(fit).
  (3) Bird holdout gains a Brier line (0.519) tied to bace_validate
  vocabulary; amphibian CV gains the zRMSE bridge (0.400, computed inline).
  (4) Future-rerun patches: 03_run_full.R (amph) and 03_summarise.R (bird)
  now compute t-intervals (field names kept for compatibility).
- **ms figure verified reproduced** (page 14 eyeballed): all four panels,
  identical headline stats; caption now calls panels (a) and (c) the same
  95% multiple-imputation interval (was "95% credible interval" over a ~92%
  quantile interval — an internal inconsistency, now fixed).
- NOT done (optional, flagged to Dan): running bace_validate() as a NEW
  complementary analysis on the case data (~4 h/rep amphibian experimental
  subset; ~hours bird) — a fresh analysis, not an integration.
- README quick-start updated with both functions (same day, earlier).

## 2026-11-08: CRAN-ready at 0.2.0 (not submitted)

Full checklist: `.agents/plans/cran-submission.md`. Headline: `R CMD check
--as-cran --run-donttest` is **0 errors / 0 warnings / 0 notes**, tarball
`BACE_0.2.0.tar.gz` (223 KB) builds clean with no dev/ms/vignettes leakage.
Version bumped 0.1.0.9000 -> 0.2.0 (0.1.0 is already on GitHub and pinned by
the bird case study, so it cannot be reused for different code).

Two genuine defects found and fixed, both of which CRAN would have caught:
(1) graphics/grDevices/stats/utils were imported in NAMESPACE and parallel::
and tools:: called in code, with NONE declared in DESCRIPTION Imports;
(2) get_imputation_intervals() had no examples, and the bace_validate()
example ran 6.4 s (over CRAN's 5 s flag) -- now 0.5 s and 2.9 s.
Also updated: NEWS header, cran-comments.md (incl. why 0.2.0), inst/WORDLIST.

NOT sent: no win-builder, no submission -- both need Dan's explicit go.
Open question for Dan: the package ships no vignette (tutorial is a 12 MB
website); acceptable but a reviewer may ask.

## 2026-11-08 (later): feature 4b implemented — per-phase MCMC settings

`bace()` gains `final_nitt`, `final_thin`, `final_burnin` (default NULL =
reuse the chained-phase values, byte-identical to previous behaviour).
Resolved once in bace() with `%||%` immediately after the unused-dots
warning, so BOTH bace_final_imp() call sites (converged and
skip_conv/non-converged) get the same values; `nitt_cat_mult` still applies
downstream as usual. This removes one of the two reasons the case studies
bypass bace() (the other, memory/chunking on 18k rows, remains).

Tests: new `tests/testthat/test-final-mcmc-args.R` (13 expectations) —
defaults identical to chained; final_* reach only the final phase; forwarding
on the skip_conv branch; per-formula lists; and a real tiny fit proving the
phases stored (400-100)/2 = 150 vs (600-100)/2 = 250 draws. Suite 1639/0.

Docs: roxygen @param, NEWS bullet under 0.2.0, README feature list. Also
corrected an now-outdated sentence in the bird case study that said bace()
"uses one set of MCMC settings for both phases".

### 4b adversarial review (2026-11-08) — one silent breaking change caught

3-agent review (correctness-by-running / docs consistency / test adequacy).
The forwarding mechanism itself survived: with final_* NULL the arguments
reaching both phases are byte-identical to HEAD (verified by diffing spies
against a HEAD copy), the retry ladder and phylo_signal paths leak nothing,
and nitt_cat_mult applies once per phase. Mutation testing confirmed the new
tests catch all five seeded forwarding bugs. Fixed from findings:

- **MUST: the three arguments had been inserted mid-signature** (between
  `burnin` and `runs`), silently rebinding positional calls:
  `bace(..., 600, 2, 100, 3, 4)` read runs=3,n_final=4 as
  final_nitt=3,final_thin=4 with runs/n_final back at defaults — no error.
  Moved to the END of the formals; positional compatibility re-verified
  against the 0.1.0 order. Residual, documented in NEWS: `fixformula` can no
  longer be abbreviated to `f`/`fi` (`fix` and longer still work).
- **MUST: bace_validate(fit) re-ran the CHAINED phase at FINAL chain length.**
  bace_final_imp's $settings necessarily records the settings IT received.
  bace() now stores both phases on its result; bace_validate reads chained
  from there (falling back for older objects) and forwards final_* through
  .bace_validate_run and all three entry methods. Verified end to end.
- Bad final_* now validated in bace() BEFORE the chained phase runs, naming
  the argument the user wrote (was: hours wasted, error said "nitt").
- Warning when final_nitt < nitt (the convergence check says nothing about
  the final chains).
- Docs corrected: nitt_cat_mult multiplies nitt/burnin but never thin; the
  worked example now scales thinning with the chain (it was a multi-GB memory
  trap); NEWS no longer claims both case studies bypassed bace() for this
  reason (only the bird one did — the amphibian splits for memory); the bird
  case study's stale justification removed; tutorial gained a bullet and two
  "increase nitt" tips now point at the final-phase arguments; ms.qmd gained
  one sentence; plan header + roadmap updated.
- Tests: 1653 -> +7 (formals order pin, positional-call compat, early
  validation naming, short-chain warning, retry-ladder chained settings,
  settings round-trip incl. legacy objects, and the real fit now varies all
  three settings so honouring only final_nitt fails).

## 2026-11-09: the convergence check is broken — flag means nothing; fix in progress

Found while setting up the new simulation (dev/simulation_2026-11): BACE's
chained-equations convergence verdict came back FALSE on essentially every
run (2/20 on the 20-replicate test, 0/17 across every pilot setting), yet
interval coverage was nominal (0.95) and the verdict showed NO relationship
with bias, coverage, width or recovery across datasets. Diagnosed, through two demonstrations:

1. **Empirical**: 20 replicates, converged vs not: |bias| 0.178 vs 0.125,
   width 0.639 vs 0.637, coverage 1.00 vs 0.944. Nothing significant; signs
   inconsistent. The flag carries no information about result quality.
2. **Mechanical**: a series converged BY CONSTRUCTION (iid noise, no drift)
   passes the whole-dataset check 3.4% of the time at campaign settings
   (runs = 5). Matches the observed 2/20 and 0/17.

Root causes in `.assess_summary_convergence()` (R/convergence_diagnostics.R):
- **%-change criterion divides by the series' own mean**
  (`abs(diff(s)/(s[-n]+1e-10))`): for mean-centred traits (standardised
  gaussians — the common case) the denominator is ~0, the "percentage"
  explodes, and the criterion can NEVER pass (0% pass rate at trait mean 0;
  needs |mean| ≳ 2 to behave). This is why x1 and cat3 fail everywhere.
- **ACF criterion** correlates 4 points against a |r| < 0.3 cutoff — rejects
  ~60% of perfectly settled chains (sampling SD of r on 4 points ≈ 0.5).
- **Geweke is silently disabled** (`geweke_result$converged <- NA`,
  "Temporary turn-off"), so only 3 live criteria remain, one unpassable.
- Combination: per-var needs >50% of live criteria, overall needs >70% of
  variables (= 5 of 6). Hence ~3% overall pass rate on healthy chains.

User impact: bace() defaults (skip_conv = FALSE, max_attempts = 3) retry at
runs 5→8→11 = 24 run-units instead of 5 (~4.8× compute) and then falsely
warn. More chain cannot help (pilot: 5→15 runs and 50k→150k final nitt all
0/17 converged).

Also today, campaign economics settled:
- **n=500 cost measured: 109.5 min** (one full bace() at campaign settings,
  while sharing the machine with 5 pilot jobs) — UNDER the 2.2 h budget;
  campaign back to ~16.5 days on 200 cores. Result quality at n=500: bace
  bias +0.0006, se 0.062 vs complete-case 0.115, all arms covered.
- **Wall time varies ~5× between datasets** at identical settings (21.6 vs
  115.6 min at n=100). Not yet explained; matters for per-job wall limits.
- Chain-length pilot (17/40 done, parked): coverage 17/17, bias ≈ 0 at ALL
  settings — runs 5/10/15 × final_nitt 50k/150k change nothing but cost.
  Verdict: the campaign settings (runs 5, nitt 50k) are already enough for
  nominal coverage; no evidence longer chains buy anything.

**FIXED same day — BACE 1.0.1.** The check was replaced wholesale (see
`.agents/investigation-convergence-2026-11.md` for the full record):
scale-free stability check on the mean AND SD of imputed values, first half
of the chain treated as warm-up, three named pathologies (drift /
oscillation / stuck, each with its remedy); bace() now retries only on
drift. Calibrated: settled chains pass >=99% at runs=10 wherever their mean
sits; each pathology caught >=99%; all 18 real harvested chain-series pass
(their measured coverage was 0.95). Old per-criterion / energy / wasserstein
methods defunct with explanatory errors; pct_change_threshold/lag/alpha
removed (drift_threshold/stuck_ratio added); plot type "steps" replaces
"pct_change". Suite: 1775 pass / 0 fail. The old campaign's convergence
percentages (28/51/82%) measured the bug and are struck from the plan; the
chain-length pilot is retired (its 17 finished runs: coverage 17/17 at every
setting — campaign settings already suffice).
