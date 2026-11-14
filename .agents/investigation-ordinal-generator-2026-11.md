# The simulated ordinal trait collapsed when phylogenetic signal was low

**Found:** 2026-11-13, by Shinichi, in the first runs on BACE 1.0.3: a four-level
ordinal trait came out with one level holding a single species and another empty
(seed 1 at signal 0.3: 0 / 47 / 52 / 1), and one dataset ran for an hour where its
neighbours took nineteen minutes.
**Fixed:** 2026-11-14, in BACE 1.0.4 — after a first fix was rejected by review.

---

## What went wrong

The simulation asks for an ordinal trait with four levels. At low phylogenetic signal
it was getting a two- or three-level trait with a stray species or two. Measured on
the 1.0.3 generator, 40 seeds per setting, 100 species, traits related at slope 0.5,
whole clades hidden:

| signal asked for | fewer than 4 levels in truth | a level with nobody observed, after masking |
|---|---|---|
| 0.3 | 27 / 40 | 31 / 40 |
| 0.7 | 7 / 40 | 9 / 40 |
| 0.9 | 2 / 40 | 4 / 40 |

So in a third of the design, the "does BACE handle ordinal traits" part of Claim 3 was
being tested on a trait that was not ordinal in any useful sense. (An earlier draft of
this note gave 25 / 18 / 4 and 6 / 2 / 2 for columns like these; those counted only the
levels that were populated in truth, which is not what the labels said. The review
caught it.)

## Was it the missingness?

No. The first column is the **truth**, before any value is hidden. Hiding whole clades
made it worse — a thin level is easy to empty — but it did not cause it.

## Why it happened

To make an ordinal trait, the generator draws a hidden continuous value for each
species and cuts it at thresholds. The thresholds were **fixed constants**: −2.25, 0 and
+2.25 for four levels.

But the spread of the hidden value is not fixed. In this simulation the ordinal's hidden
value is `0.5 × x1` plus its own species, phylogenetic and residual effects, and `x1`'s
own spread grows with the signal asked for. Its variance is

```
                 2 × 0.5²  +  0.2             0.7
   var  =  ------------------------------  =  -----------
                  1 − phylo_signal            1 − signal
```

(2 is `x1`'s variance before signal, 0.2 the ordinal's own; both are `sim_bace()`
defaults.) So the standard deviation is 1.00 at signal 0.3, 1.53 at 0.7 and 2.65 at 0.9
(realised 0.96 / 1.40 / 2.36, a touch lower because one tree is one draw), and for an
ordinal with no parent 0.53 / 0.82 / 1.41. At signal 0.3 about two species in a hundred
fall beyond ±2.25; with no parent, none. The thresholds never knew what scale they were
cutting.

The binary trait does not have this problem (it is drawn with a probability, not cut at
a constant), and the continuous and count traits have no thresholds. The count trait has
a related property worth knowing, though: its hidden scale also grows with signal, and
after the log link that makes it wildly overdispersed at signal 0.9 (variance-to-mean
about 1600, maximum 5,778, 40% zeros) — not a bug, but a very different count trait from
the one at signal 0.3.

## The hour-long run

Seed 1 of `n100_ps0.7_good_covariate_clade`: truth 1 / 43 / 55 / 1, and after clade
masking the ordinal was **observed as 0 / 29 / 39 / 0**. Two levels with no species in
them. MCMCglmm estimates a cutpoint between each pair of adjacent levels; between two
empty levels there is nothing to anchor it, and the sampler crawls. Same root cause as
the collapse, not a separate problem.

## The first fix, and why it was wrong

The first fix cut at the **quartiles of the realised hidden values**: every dataset then
has exactly 25 / 25 / 25 / 25, whatever the scale. It passed every test and both mock
runs. The adversarial review rejected it, with numbers:

- Fixing the sample counts makes the **hidden cells' level counts the exact complement
  of the observed cells'**. Across seeds the correlation between the observed share of a
  level and its hidden share was −1.00 in every setting tried. Under constant thresholds
  it is +0.2 to +0.9, because a tree that shifts together shifts its observed and hidden
  cells together — which is information an imputer can legitimately use.
- So "guess the commonest observed level" scores **below chance** on the hidden cells:
  0.11–0.21 against 0.25. And an imputer that estimates its cutpoints from the observed
  cells — BACE's threshold model does — is pushed the wrong way: a phylogenetic threshold
  model scored 0.03–0.12 lower under the quartile cut than under constant thresholds on
  identical hidden values and masks.
- The paper's only ordinal metric is class accuracy. A data-generating rule that biases it
  downward and moves the baseline is not something to leave undisclosed.

## The fix

Cut at **population thresholds**: constants chosen from the known distribution of the
hidden value, not from the dataset. The campaign supplies the exact quartiles of its hidden
scale, `qnorm(c(.25, .5, .75)) × sqrt(0.7 / (1 − signal))` (0.2 instead of 0.7 when the
traits are unrelated), so each level has 25 species in expectation at every signal, and
hidden and observed cells are draws from the same population. `sim_bace()` gained an
`ordinal_thresholds` argument for this; when nothing is supplied its default is the
quartiles of the variable's **own** hidden scale, which is exact for an ordinal with no
parent and, with a parent, gives the outer levels more than their share rather than fewer.
Never the old constants.

Measured on the same datasets as the first fix (40 seeds per setting, n = 100):

| | sample quartiles (rejected) | population thresholds (1.0.4) |
|---|---|---|
| levels in truth | exactly 25 each | vary; smallest 6 at signal 0.3, 1–3 at 0.7, 0 at 0.9 |
| empty level in truth, signal 0.9 | never | 2 / 40 |
| level unobserved after clade masking, signal 0.9 (0.7) | never | 2–4 / 40 (1 / 40) |
| same at n = 500 | never | 0 / 10 |
| correlation, observed vs hidden share | −1.00 | +0.2 to +0.9 |
| "commonest observed level" floor (chance 0.25) | 0.11–0.21 | 0.27–0.43 |

The price is honest: at strong signal a whole tree can shift together, so about 1
dataset in 20 at signal 0.9 and n = 100 has an ordinal level with no species at all, and
1 in 10–20 has a level unobserved after masking. A level nobody observes cannot be imputed
by any method, and those are the slow datasets. The campaign now records the level counts
of every discrete trait per dataset, before and after masking (`realised$levels_truth`,
`levels_observed`), so the paper can report ordinal accuracy with and without them.

The cut uses no random numbers under any rule. Verified on 200 datasets (five scenarios ×
40 seeds), 1.0.3 against the final 1.0.4:

| | identical before / after (by position) |
|---|---|
| tree | 200 / 200 |
| `x1`, `x2`, `cnt`, `driver`, `aux` | 200 / 200 |
| `bin` | 200 / 200 |
| which cells are hidden | 200 / 200 |
| `ord` | changes — the point |
| `cat3` | changes where `ord` does — it is simulated *from* `ord` |
| species names | same set, now in C order on every machine (see below) |

## What this changed

- **The version, to 1.0.4.** A fixed build that still called itself 1.0.3 would pass the
  campaign's guard. `run_one.R`, `run_all.sh` and the README now require 1.0.4.
- **Everything already run on 1.0.3 must be rerun.** Not only the ordinal results: the
  campaign imputes every trait from every other, so a different `ord` on the right-hand
  side changes the imputation of `x1`, `x2`, `cnt` and `bin` too. Only timings carry over.
- **The probe** prints the ordinal, binary and categorical codes themselves, and its
  Section A for 1.0.4 is committed as `reproducibility_probe_reference.txt`.

## How we stop it coming back

`tests/testthat/test-ordinal-generator.R`:

1. The cutter's contract under all three rules: population thresholds cut exactly where
   told and are validated; the quantile rule populates every level at any scale;
   monotone; no random numbers consumed; rounding-level noise cannot move a code; the
   old constants unchanged for those who ask; the empty-level warning names its cause.
2. `sim_bace()`'s default thresholds equal the own-scale formula to 1e-12, and explicit
   thresholds reach the variable (the codes equal the cut of the hidden value rebuilt from
   the returned random effects).
3. **A golden test for every trait type**: the campaign's exact recipe for one seed —
   tree, two continuous traits, a count, a binary, the ordinal at its population
   thresholds, and the categorical response — pinned to recorded values in any locale.
   The 1.0.3 golden test covered only the continuous draw; this one found the locale
   problem below.

## A second thing the golden test caught: species names depended on the locale

Running the new golden test under `R CMD check`, which uses the C locale for
tests, every value, count and code reproduced exactly, but the species names did
not: the first three species were `Z42dv l9qwa U5Lyc` in an interactive en_US
session and `yjMQ9 SxiO7 nSPob` under check. The generator sorted its random
species names with the session's collation, and the two locales order mixed-case
strings differently (en_US puts `1wIAO` before `1WsN2`, C puts it after).

The values, the tree's structure and the hidden cells were identical in both
locales. Only the strings attached to each position differed. But a dataset
whose species are called different things on different machines is not the same
dataset, and it makes two result files impossible to diff.

The obvious fix, sorting the names in C order on every machine, **silently
changed the values** in an en_US session: same sums, different cases
(`x1[1:3]` went from 2.95 / −3.09 / 2.36 to −1.63 / −0.69 / 1.40). The reason is
worth recording. The random effects are drawn in tree-tip order, and the matrix
that hands each case its species' effect was built from `factor(species)`, whose
level order follows the session's collation. The pairing was right only because
two independent sorts happened to agree. Change one and every case quietly
receives another species' effect, with nothing to warn.

So the fix is two lines, not one: names are sorted in C order on every machine
(`method = "radix"`), and the design matrix names its levels from the tree's tip
labels, so each case is matched to its tip **by name**. The golden values are
unchanged in every locale, the labels are now the same everywhere, and a test
runs the generator under the C and en_US locales and requires identical output.
This is the same defect class as the rest of the 2026-11 audits: a pairing that
was re-derived (by sorting twice) instead of recorded (by name).

## Checks before commit

Full suite on the final code: **2119 pass / 0 fail / 0 error / 0 warn / 0 skip**.
`R CMD check --as-cran`: **0 errors / 0 warnings / 2 NOTEs**, both about the environment
(the coveralls badge URL returns 403 under the "new submission" check; the local HTML
Tidy is old). Under check, which runs the tests in the C locale, 1551 pass / 0 fail.
Two rounds of adversarial review (four lenses, then three, each finding re-run by a
skeptic): the first rejected the sample-quantile fix and found the number and locale
defects described above; the second reviewed the population-threshold design, the
corrected numbers and the new tests. The second round's break-it reviewer found nothing
(nine checks passed, including the variance formula and the empty-level counts); its two
other reviewers were cut off by the account's usage limit, so the test-strength check was
done by hand instead: nine deliberate breakages of the new code (supplied thresholds
ignored; the default forgetting the phylogenetic variance, or the intercept; the cut using
≥ instead of >; an unmatched `beta_ix` falling back to a random draw; interaction names or
species names sorted by locale; the one-species guard removed; the design matrix paired by
locale order) each made at least one test fail, and the untouched copy passed before and
after.

## Mock runs on the final build

`run_one.R` on the installed 1.0.4 (R 4.5.3, this Mac), campaign settings (M = 20, nitt
50,000, burnin 10,000, thin 25, runs 7), for the two datasets that failed Shinichi.

**Seed 1 of `n100_ps0.3_typical_clade` — the collapse.** On 1.0.3 its ordinal was
0 / 47 / 52 / 1 in truth and 0 / 30 / 37 / 1 observed. Now 15 / 32 / 27 / 26 in truth
(population cuts at ±0.674 and 0), observed 10 / 20 / 17 / 21, 32 cells hidden. Completed in
**22.5 minutes** (1,348 s wall, 1,343 s in `bace()`), no failure. Ordinal accuracy on the
hidden cells 0.48 (chance 0.25; it was 0.39 on the rejected sample-quartile build, as the
review predicted), categorical 0.45, binary 0.50; `x1` zRMSE 0.54; per-cell interval coverage
`x1` 0.96, `x2` 0.97, `cnt` 1.00. Downstream slope of `x2` on `x1` (truth 0.5): oracle 0.47,
complete-case 0.41, BACE 0.45, all three intervals covering. The chained-equations stability
check (recorded, not acted on — the harness runs `skip_conv = TRUE`) returned `ok` for four
variables and `oscillation` for `cat3` and `ord`; the campaign records that verdict per
dataset, so its frequency will be measured rather than guessed.

**Seed 1 of `n100_ps0.7_good_covariate_clade` — the hour-long run.** On 1.0.3 its
ordinal was 1 / 43 / 55 / 1 in truth and **0 / 29 / 39 / 0 observed**, two levels with
nobody in them. Now 15 / 29 / 38 / 18 in truth (cuts at ±1.03 and 0), observed
9 / 20 / 30 / 9, 32 cells hidden. Completed in **23.9 minutes** (1,434 s wall, 1,432 s in
`bace()`) against roughly an hour before, no failure, stability check `ok` for all six
variables. Ordinal accuracy on the hidden cells 0.28, categorical 0.60, binary 0.59; `x2`
zRMSE 0.69, `x1` 1.15, `cnt` 1.90 — a whole clade hidden at n = 100 is a hard dataset, which
is what this arm is for; one dataset's ordinal score (32 cells) says little on its own, and
the campaign measures the distribution. Per-cell coverage `x1` 0.86, `x2` 1.00, `cnt` 0.92.
Downstream slope (truth 0.5): oracle 0.58, complete-case 0.51, BACE 0.57, all covering.

Both runs used R 4.5.3 through its own binary, because R on this Mac was upgraded to
4.6.2 with an empty package library at midday (two earlier attempts died on "no package
called ape" for that reason, after a transient typo in the runner had killed the pair
before them at the scoring step, 24 minutes in). The runner itself is unchanged. Once the
4.6 library was populated, the golden recipe on R 4.6.2 reproduced every value of the
4.5.3 run exactly — one more machine-independence check, this time across R versions.

## Related

- [investigation-simulator-reproducibility-2026-11.md](investigation-simulator-reproducibility-2026-11.md)
  — the day before: same seed, different data on different machines (1.0.3)
- [investigation-audit-2026-11.md](investigation-audit-2026-11.md) — the numerical audits
- `dev/simulation_2026-11/README.md` — "Things that will bite"
