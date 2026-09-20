# CRAN release plan (v0.1.0)

Created 2026-09-19. Status of getting BACE from "CI-green research package" to a
CRAN submission. Baseline: `R CMD check --as-cran` is already **Status: OK**
(0 errors / 0 warnings / 0 notes) on CI across macOS, Windows, Ubuntu
(release, devel, oldrel-1); the whole check runs ~2 min (MCMC tests use tiny
`nitt`; tarball excludes vignettes/dev/ms/bib).

## Step 1 — DESCRIPTION + NEWS + \value docs  (2026-09-19, in progress)

- [x] DESCRIPTION overhaul: Title fixed to canonical "by Chained Equations"
  wording + domain; Version 0.0.0.9000 → 0.1.0; full Description paragraph
  with method references and DOIs (Hadfield 2010; Zhou & Reiter 2010; Rubin
  1987; Barnard & Rubin 1999); URL + BugReports added.
- [x] NEWS.md created (0.1.0 feature summary).
- [x] `@return` tags added for the 18 exported functions/methods whose Rd
  lacked `\value` (8 plot fns, 9 print/summary methods, pipe re-export);
  each claim adversarially verified against the function body. Two
  verification catches fixed by adding explicit `invisible(NULL)`:
  `plot_energy_convergence` and `plot_convergence_summary` ended in
  `grid()`, which since R 4.4 invisibly returns grid-line positions, so
  "returns NULL invisibly" would have been false. `devtools::document()`
  re-run; tests 1337 pass / 0 fail.
- Deliberately skipped: `\value` for `.extract_imputed_datasets` /
  `.mcmcglmm_fixef_idx` (become `@noRd` in Step 2) and `BACE-package.Rd`
  (overview page; no `\value` needed).
- [x] Step 1b (commit 73be417): CRAN conventions done. print methods return
  `x` invisibly; `summary.bace_pooled_MCMCglmm` returns the summary invisibly
  (its step-1 doc wrongly said "printed, not returned");
  `plot.bace_convergence` got a real \value + `invisible(NULL)`; ALL roxygen
  non-ASCII converted to ASCII; `tools::checkRd()` clean on all 90 Rd pages
  (also fixed Lost-braces itemize + a period-terminated title);
  `bace_option_defaults` docs corrected (claimed verbose FALSE / gelman 1,
  code returns verbose TRUE / gelman 2). Tests 1337/0.
- Remaining known cosmetics (deliberate): ~735 non-ASCII chars in R code
  strings/comments (cat() box-drawing etc.) — R CMD check --as-cran does not
  flag them (CI Status: OK); revisit only if CRAN incoming complains.

## Step 2 — Trim exported surface (roadmap C4)  (DONE 2026-09-19 except decision item)

- [x] All 19 dot-prefixed internals un-exported; NAMESPACE 49 → 30 exports
  (verified: exactly 19 `export()` lines removed, 0 added). Extended the
  sweep: EVERY dot-function roxygen block is now `@noRd` (docs stay in
  source), so all 41 `dot-*.Rd` pages are gone and man/ is purely
  user-facing — this also disposed of the two `\value`-less internal pages.
  `phylo_signal_summary.R`'s 14 internals already had `@noRd`.
- [x] Tests: `tests/testthat/helper-internals.R` aliases the 19 via
  `BACE:::` once, so the 78 bare call sites (test-prep_functions,
  test-species_effects, test-track-a-hardening) run unchanged; consistent
  with the 132 existing `BACE:::` uses in the suite. 1337 pass / 0 fail.
- [x] Verified no `BACE::.` double-colon usage anywhere; dev/ scripts using
  internals rely on `devtools::load_all()` and are unaffected for CRAN.
- [x] **Decision (2026-09-20, Dan): tight API.** All 15 sim-machinery /
  converter exports internalised (`@noRd`; man pages removed; helper-aliased
  for tests). Final public API = 23 exports: the pipeline (bace, bace_imp,
  bace_final_imp, assess_convergence, pool_posteriors, pool_mi,
  with_imputations), accessors (get_pooled_model, get_imputed_data),
  phylo_signal_summary, options (bace_options, bace_option_defaults),
  simulation (sim_bace, sim_tree), the 8 plot_* functions, and the pipe.
- [x] Gate: CI R-CMD-check --as-cran green on steps 1, 1b, 2.

## Step 3 — Examples  (DONE 2026-09-20)

- [x] Zero `\dontrun` remains. 9 MCMC pages got minimal verified `\donttest`
  pipelines (adapted from a pre-validated ~5 s template; sim_bace n=40 +
  nitt=600); pool_mi / with_imputations / bace_options got fully runnable
  unwrapped examples (no MCMC); phylo_signal_summary gained its first
  example. All examples EXECUTED locally via
  `devtools::run_examples(run_donttest = TRUE)`: exit 0, ~9 s wall total —
  comfortably inside CRAN limits even though --as-cran runs donttest.
- [x] Fixed a roxygen link warning (`[k,k]` linkified in
  phylo_signal_summary docs).

## Step 3 — Examples

- [ ] Replace `\dontrun{}` (17 pages) with runnable examples where fast, or
  `\donttest{}` for MCMC-heavy ones (CRAN policy: `\dontrun` only for code
  that genuinely cannot run).

## Step 4 — Final gates  (2026-09-20)

- [x] **Maintainer decision (Dan)**: `cre` switched to Daniel Noble.
  Authors@R reordered to match the paper team: Noble (aut, cre), Mizuno
  (aut), Pottier (aut), Drobniak (aut), Nakagawa (aut) — Ayumi Mizuno and
  Patrice Pottier added as package authors per Dan.
- [x] `inst/CITATION` updated: package title corrected ("by Chained
  Equations ... Phylogenetic Comparative Data"), five authors in the order
  Noble, Mizuno, Pottier, Drobniak, Nakagawa; Hadfield (2010) kept as
  second entry. Parses and renders cleanly. Add the preprint DOI as a third
  entry when available.
- [x] `cran-comments.md` written (environments, 0E/0W/0N, first-submission
  notes).
- [x] Spelling: `Language: en-GB` declared; genuine typo fixed
  (`mypkg.*` -> `BACE.*` template leftover in bace_options docs);
  `inst/WORDLIST` (90 jargon/surname terms); spell check now clean.
- [x] `urlchecker::url_check()`: all URLs correct.
- [x] win-builder round 1 (R-devel, 2026-09-20): **2 NOTEs**.
  (a) top-level `cran-comments.md` + stray local `Rplots.pdf` leaked into
  the tarball — FIXED (commit f7755d3: both .Rbuildignore'd, plus repaired
  a fused `^notes/` pattern; rebuilt tarball verified clean);
  (b) "possibly misspelled" DESCRIPTION words are reference surnames
  (Hadfield, Zhou, Reiter) + "phylogenetically" — unavoidable, explained in
  cran-comments.md.
- [x] win-builder rounds 2-3 (2026-09-20): BOTH still reported the
  'Rplots.pdf' / 'cran-comments.md' top-level NOTE — **proven to be
  win-builder-side residue, not the package**. Evidence chain: (i) the
  uploaded tarball (sha1 0f9d0acac80742e0f6d04b770357e4190c3826bf, 101
  members) verifiably contains neither file; (ii) the 06:35 result
  demonstrably checked that tarball (Drobniak-second Authors@R in its built
  zip); (iii) the identical tarball through a pristine local
  `R CMD check --as-cran` produces NO non-standard-files complaint — only
  the unavoidable incoming NOTE (new submission + reference surnames).
  Mechanism: round-1's genuinely dirty tarball unpacked into win-builder's
  reused per-package work dir; later unpacks overwrite matches but never
  delete strays (a check run cannot invent cran-comments.md). All other
  Windows R-devel checks pass (examples, tests 16 s, Rd, manual).
- [ ] **Decision (Dan, 2026-09-20): wait ~72 h** for win-builder's auto-
  cleanup (files removed ~2026-09-23), then re-upload for a clean log and
  submit. Re-upload: rebuild from the repo (state is fully committed) with
  `R CMD build --no-manual .` and
  `curl -T BACE_0.1.0.tar.gz ftp://win-builder.r-project.org/R-devel/`;
  verify the tarball manifest first (expect 101 members, no Rplots.pdf /
  cran-comments.md). Note: win-builder FTP was intermittently unreachable
  on 2026-09-20; retry politely if it times out.
- R-hub skipped deliberately: CI already covers ubuntu (devel/release/
  oldrel-1) + macOS + windows with --as-cran, and R-hub v2 needs
  interactive GitHub auth. Revisit only if CRAN flags a platform we
  don't cover (e.g. special Solaris/clang builds).
- [ ] **Submit**: after win-builder comes back clean —
  https://cran.r-project.org/submit.html with the built tarball and
  cran-comments.md; Dan confirms the email. (Or `devtools::release()`
  interactively.)

## Open decisions

- **Maintainer (`cre`)**: currently Shinichi. CRAN corresponds only with the
  maintainer (including the submission confirmation email). If Dan drives the
  release, consider switching cre to Dan. NOT changed unilaterally.
- Version 0.1.0 chosen (conventional first release); bump to 1.0.0 with the
  manuscript if preferred.
- CRAN track is independent of Tracks B/C (benchmark/manuscript); a preprint
  DOI can be added to Description/CITATION in a later patch release.
