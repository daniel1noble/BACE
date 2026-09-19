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
- [ ] **Decision (open)**: borderline non-dot exports — `mnom_liab2cat`,
  `ordinal_liab2cat`, `generate_default_beta_matrix`,
  `print_sim_bace_summary`, `sim_bace_gaussian/poisson/binary`, `sim_tree`:
  keep exported (documented, benign) or internalise for a tighter v0.1 API?
- [x] Gate: CI R-CMD-check --as-cran on the commit.

## Step 3 — Examples

- [ ] Replace `\dontrun{}` (17 pages) with runnable examples where fast, or
  `\donttest{}` for MCMC-heavy ones (CRAN policy: `\dontrun` only for code
  that genuinely cannot run).

## Step 4 — Final gates

- [ ] `cran-comments.md` (test environments, R CMD check results, first
  submission note).
- [ ] `spelling::spell_check_package()`, `urlchecker::url_check()`.
- [ ] win-builder (`devtools::check_win_devel()`) + R-hub.
- [ ] `inst/CITATION` (cite package + Hadfield 2010; add preprint DOI when
  available).
- [ ] Submit; `cre` must confirm the submission email.

## Open decisions

- **Maintainer (`cre`)**: currently Shinichi. CRAN corresponds only with the
  maintainer (including the submission confirmation email). If Dan drives the
  release, consider switching cre to Dan. NOT changed unilaterally.
- Version 0.1.0 chosen (conventional first release); bump to 1.0.0 with the
  manuscript if preferred.
- CRAN track is independent of Tracks B/C (benchmark/manuscript); a preprint
  DOI can be added to Description/CITATION in a later patch release.
