# CRAN comments for BACE 0.1.0

## Test environments

* Local: macOS (Apple Silicon), R 4.5.3 (aarch64-apple-darwin20)
* GitHub Actions (`R CMD check --as-cran`, failing on any warning):
  * ubuntu-latest: R devel, R release, R oldrel-1
  * macOS-latest: R release
  * windows-latest: R release
* win-builder: R devel

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the standard "New submission" note. It also flags the coverage
badge URL (https://coveralls.io/github/daniel1noble/BACE?branch=main) as
possibly invalid (403): coveralls.io blocks automated requests; the badge
resolves normally in a browser.

## Comments

* First submission of BACE to CRAN.
* All examples are executable. Examples that fit MCMC models are wrapped in
  \donttest{} and use deliberately small simulated datasets and short chains;
  the full set (including \donttest) completes in roughly 10 seconds locally.
* The test suite (1424 tests) completes in about two minutes: MCMC-dependent
  tests use very small chains, and the heaviest scenarios are additionally
  guarded with skip_on_cran().
* Words flagged by spell checking are statistical/phylogenetics terminology,
  author surnames, and R package names (see inst/WORDLIST). In particular,
  the DESCRIPTION words flagged by win-builder — Hadfield, Zhou, Reiter —
  are surnames from the method references, and "phylogenetically" is
  standard phylogenetics terminology.
