# CRAN comments for BACE 1.0.0

## Test environments

* Local: macOS (Apple Silicon), R 4.5.3 (aarch64-apple-darwin20) —
  `R CMD check --as-cran`, including `--run-donttest`
* GitHub Actions (`R CMD check --as-cran`, failing on any warning):
  * ubuntu-latest: R devel, R release, R oldrel-1
  * macOS-latest: R release
  * windows-latest: R release
* win-builder: R-devel, R-release and R-oldrelease (Windows Server 2022 x64)

## R CMD check results

0 errors | 0 warnings | 1 note

The only note is the standard new-submission note from the CRAN incoming
feasibility check:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Daniel Noble <daniel.noble@anu.edu.au>'

New submission

Possibly misspelled words in DESCRIPTION:
  Hadfield (13:65)
  Reiter (18:46)
  Zhou (18:37)
  phylogenetically (23:34)
```

Hadfield, Zhou and Reiter are surnames of the authors of the method references
cited in the Description (Hadfield 2010 <doi:10.18637/jss.v033.i02>; Zhou and
Reiter 2010 <doi:10.1198/tast.2010.09109>). "phylogenetically" is standard
phylogenetics terminology. All four are listed in `inst/WORDLIST`.

The local and GitHub Actions checks report 0 errors | 0 warnings | 0 notes.

## Comments

* This is the first submission of BACE to CRAN. It is numbered 1.0.0 rather
  than continuing the 0.x series because the user-facing API is now stable and
  this release accompanies the methods manuscript; version 0.1.0 was
  distributed only from GitHub.
* All examples are executable. Examples that fit MCMC models are wrapped in
  `\donttest{}` and use deliberately small simulated datasets and short chains.
  The full example set, including `\donttest`, runs in well under a minute.
* The test suite completes in about 90 seconds: MCMC-dependent tests use very
  small chains, and the heaviest scenarios are additionally guarded with
  `skip_on_cran()`.
* The package writes files only when the user explicitly supplies a
  `checkpoint_dir`; no example, test or default path writes outside
  `tempdir()`.
