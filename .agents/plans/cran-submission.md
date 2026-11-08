# CRAN submission checklist — BACE 0.2.0

Prepared 2026-11-08. Nothing has been sent anywhere: win-builder and the CRAN
submission both need Dan's explicit go at the time (standing rule).

## State: ready

`R CMD check --as-cran --run-donttest` on macOS, R 4.5.3: **0 errors, 0
warnings, 0 notes**. Tarball `BACE_0.2.0.tar.gz` (223 KB) builds clean and
ships only DESCRIPTION, NAMESPACE, NEWS.md, R/, README.md, build/, inst/,
man/, tests/ — no dev/, ms/, vignettes/, .agents/ or index.html.

## Why 0.2.0 and not 0.1.0

0.1.0 was distributed from GitHub and is pinned by a published analysis
(the bird case study tells readers to install commit dc79be1). This release
adds exported functions, so it takes a new number rather than reusing 0.1.0
for different code. `cran-comments.md` states this so the reviewer is not
confused by a "new submission" at 0.2.0.

## Fixed during preparation (2026-11-08)

- **Undeclared dependencies** (would fail CRAN's stricter configurations):
  NAMESPACE imported from graphics, grDevices, stats, utils and the code
  calls parallel:: and tools::, none of which were in DESCRIPTION. All six
  now declared in Imports.
- **`get_imputation_intervals()` had no examples** — added (0.5 s).
- **`bace_validate()` example took 6.4 s**, over CRAN's 5 s flag on machines
  slower than ours — trimmed to 2.9 s (lighter MCMC settings, both call
  forms kept).
- Version/NEWS/cran-comments updated to 0.2.0; new terminology added to
  inst/WORDLIST (zRMSE, Honaker, ultrametricity, BLUPs, CMRG, MCC, imputer,
  multiPhylo, preds).

## Verified clean

URLs (only the coveralls badge 403s to robots — resolves in a browser, and
cran-comments explains it); no \dontrun; no examples writing to the user's
filesystem; CITATION uses modern bibentry with meta$Version; no R-version-
specific syntax (no native pipe, no lambda shorthand), so no Depends: R
needed and GHA oldrel-1 passes; non-ASCII only in code comments; tests use
at most 2 cores and skip_on_cran() guards the heavy ones.

## Steps at submission time (each needs Dan's go)

1. Push the release commit; confirm the 5-platform GitHub Actions check is
   green.
2. win-builder R-devel (and optionally R-release):
   `devtools::check_win_devel()`. Expect the spell-check note on surnames
   and terminology; compare against cran-comments.
3. Optionally `rhub::rhub_check()` for extra platforms.
4. Submit `BACE_0.2.0.tar.gz` at https://cran.r-project.org/submit.html with
   cran-comments.md pasted into the comments field.
5. On acceptance: tag v0.2.0 on GitHub, and start a 0.2.0.9000 development
   version in DESCRIPTION/NEWS.

## Known reviewer question to be ready for

The package ships **no vignette** (the tutorial is a ~12 MB rendered website
with real MCMC runs, far too heavy for CRAN; the URL field points to it).
This is acceptable and common, but a reviewer may ask. If one insists, the
cheapest answer is a short static vignette that explains the workflow with
pre-computed output rather than live fits.
