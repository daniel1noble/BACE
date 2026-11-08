# CRAN submission checklist — BACE 1.0.0

Updated 2026-11-09. Standing rule: win-builder and the CRAN submission each
need Dan's explicit go at the time. Both were authorised on 2026-11-09.

## State

`R CMD check --as-cran --run-donttest` on macOS, R 4.5.3: **0 errors, 0
warnings, 0 notes** at 1.0.0. Tarball `BACE_1.0.0.tar.gz` (241 KB) ships only
DESCRIPTION, NAMESPACE, NEWS.md, R/, README.md, build/, inst/, man/, tests/ —
no dev/, ms/, vignettes/, .agents/ or index.html.

## Why 1.0.0

Dan's call (2026-11-09), replacing the earlier plan to submit at 0.2.0. This
is the first CRAN release and the user-facing API is now stable, so it takes
a 1.x number rather than continuing the 0.x series; 0.1.0 was distributed
only from GitHub and is pinned by the bird case study (commit dc79be1).
`NEWS.md` and `cran-comments.md` both say so.

Consequence worth remembering: under semantic versioning 1.0.0 is a promise
of API stability, so a future breaking change (for instance one demanded by
manuscript review) should go to 2.0.0 rather than 1.1.0.

## Fixed during preparation (2026-11-08/09)

- **Undeclared dependencies** (would fail CRAN's stricter configurations):
  NAMESPACE imported from graphics, grDevices, stats, utils and the code
  calls parallel:: and tools::, none of which were in DESCRIPTION. All six
  now declared in Imports.
- **`get_imputation_intervals()` had no examples** — added (0.5 s).
- **`bace_validate()` example took 6.4 s**, over CRAN's 5 s flag on machines
  slower than ours — trimmed to 2.9 s (lighter MCMC settings, both call
  forms kept).
- Version bumped to 1.0.0 in DESCRIPTION and the NEWS header, with a short
  note in NEWS explaining the jump and naming the one non-additive change
  (`fixformula` can no longer be abbreviated to `f`/`fi`).
- `inst/WORDLIST` gained `Zhou` and `phylogenetically`. Note this file is
  read only by the `spelling` package — it does **not** suppress CRAN's
  "Possibly misspelled words in DESCRIPTION" note, which has its own
  dictionary. That note is explained in `cran-comments.md` instead.

## Verified clean

No `\dontrun`; no example writes anything at all (`checkpoint_dir` appears
only in `\usage`/`\arguments`, and defaults to `NULL`); the only writes in
`tests/` are inside `tmpdir_fresh()`, which sits under `tempdir()`; CITATION
uses modern `bibentry` with `meta$Version`, so it renders "1.0.0"
automatically and needs no edit per release; no R-version-specific syntax
(no native pipe, no lambda shorthand), so no `Depends: R` is needed and GHA
oldrel-1 passes; non-ASCII only in code comments; tests use at most 2 cores
and `skip_on_cran()` guards the heavy ones.

## Network constraint at Dan's location (important)

`devtools::check_win_devel()` **does not work** here: it uploads by FTP, and
the passive-mode data channel is blocked (the control connection opens and
returns 229, then the transfer times out). HTTPS to win-builder is fine.

Workaround that does work: POST the tarball to the win-builder web form.
Script kept in the session scratchpad; the recipe is

1. `GET https://win-builder.r-project.org/upload.aspx` with a cookie jar,
2. scrape `__VIEWSTATE`, `__VIEWSTATEGENERATOR`, `__EVENTVALIDATION`,
3. POST multipart with `FileUpload{N}` + `Button{N}=Upload File`.

Control mapping: **1 = R-release, 2 = R-devel, 3 = R-oldrelease**.
(4 = R-devel UCRT exists in the HTML but is inside an HTML comment, i.e.
retired — R-devel is already UCRT.)

A successful upload is confirmed by the response echoing back `File name:`,
`File Size:` and `Content type:` — roughly 235 bytes larger than the blank
form. The page shows no "thank you", so size-vs-baseline is the signal.

The CRAN submission endpoint itself (`https://xmpalantir.wu.ac.at/cransubmit/
index2.php`, which `devtools:::upload_cran` POSTs to) is HTTPS and **is**
reachable, so the FTP problem does not affect submission.

## Steps at submission time (each needs Dan's go)

1. Push the release commit; confirm the 5-platform GitHub Actions check is
   green.
2. win-builder R-devel, R-release and R-oldrelease via the web form above.
   Expect only the spell-check note on surnames and terminology.
3. Submit `BACE_1.0.0.tar.gz` at https://cran.r-project.org/submit.html with
   `cran-comments.md` pasted into the comments field. **CRAN then emails
   Dan a confirmation link that must be clicked** — the submission is not
   lodged until he does, which is a built-in final gate.
4. On acceptance: tag v1.0.0 on GitHub, and start a 1.0.0.9000 development
   version in DESCRIPTION/NEWS.

## Known reviewer question to be ready for

The package ships **no vignette** (the tutorial is a ~12 MB rendered website
with real MCMC runs, far too heavy for CRAN; the URL field points to it).
This is acceptable and common, but a reviewer may ask. If one insists, the
cheapest answer is a short static vignette that explains the workflow with
pre-computed output rather than live fits.
