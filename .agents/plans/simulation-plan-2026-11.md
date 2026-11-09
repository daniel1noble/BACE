# Simulation plan — extending the manuscript campaign

**Written for: Shinichi, Szymek and Dan.** Drafted 2026-11-09 by Dan + Claude.
Please comment before we launch, because the run takes about three weeks.

## What we are asking you to approve

The design in §3, and four decisions in §8. The code is written and tested
(§6.1), so approval is the only thing between here and launching.

In one paragraph: re-run the manuscript's simulation at **n = 100 and 500
species**, 1000 datasets each, across **three levels of phylogenetic signal**,
**three information regimes** (nothing to borrow / correlated traits / a good
auxiliary variable) and **three ways of being missing** (at random, driven by
an observed variable, whole clades). Score BACE against both the best-possible
benchmark and against what people actually do — dropping every species with a
gap. 54,000 datasets, about 72,000 core-hours, roughly 16 days on 200 cores.

**Four things need a decision before we start** (§8): what the trait–trait
correlation should be in the covariate regimes; whether three weeks is
acceptable; whether to spend ~20 core-hours first checking that the MCMC
chains are long enough at n = 100, where they currently look short; and
whether to stop giving the fully observed driver variable an imputation model
it does not need (saves about 8.5%).

---

## 1. What we are trying to show

Three things, in plain words.

**Claim 1 — BACE beats throwing data away.**
The usual thing people do with missing data is drop every species with a gap
("complete-case analysis"). That is fine if the gaps are purely random. It is
not fine when the gaps depend on something — for example, if poorly studied
species are both more likely to be missing a trait *and* different in that
trait. We want to show that when missingness is not random, complete-case
gives the wrong answer and BACE does not.

**Claim 2 — a good covariate rescues you when the phylogeny cannot.**
BACE borrows information from two places: related species, and other traits on
the same species. If phylogenetic signal is weak the first source dries up. We
want to show that is not fatal, so long as you have another trait that
correlates reasonably well with the missing one. This is the practical advice
the paper can give: what to put in your imputation model.

**Claim 3 — BACE is honest across all trait types.**
Little bias, and intervals that contain the truth 95% of the time, for
continuous, count, proportion, binary, ordinal and categorical traits alike.

---

## 2. What we already have

The manuscript's simulation is `dev/rubin_study_2026-09-26/`, run on the
corrected BACE package in 2026-10. It is complete: **750 datasets, 1,148
core-hours**, aggregated into `results_2026-10/`.

(1,148, not the 1,020 of early drafts and not the 1,177 of a later one. Every
result file records its own wall-clock in `walls[["bace_fit"]]`; nobody had read
it, so 1,020 was a guess apportioned from a total, and 1,177 came from
multiplying per-n *medians* by counts instead of adding the actual times up.
1,148.6 is the direct sum over all 750 files, and it is corroborated
independently: the campaign ran 7-wide for 166.3 hours, so 1,164 core-hours were
available and 1,148 were used — a 99% duty cycle. See §5.)

Design as it stands:

| | Current |
|---|---|
| Species | n = 100, 300, 1000 |
| Phylogenetic signal | Pagel's λ = 0.3, 0.7, 1.0 |
| Trait–trait correlation | ρ = 0, 0.5 |
| Traits | 7 (2 continuous, count, proportion, binary, 4-level ordinal, 3-level categorical) |
| Auxiliary covariate | **none** — see §3.2; the campaign has never run one |
| Missingness | **MCAR only**, 30% from each trait |
| Replicates | 50 at n = 100 and 300, **25** at n = 1000 |
| Arms | BACE (and `bace_chain`; Rphylopars/castor available but off) |
| Scoring | zRMSE, accuracy with fractional tie credit, coverage — the manuscript definitions |
| Downstream | PGLS slope + correlation, pooled by Rubin's rules |

**The pleasant surprise: most of what we want is already written.** Going
through the code rather than assuming:

| What we want | Status as of 2026-11-09 |
|---|---|
| MAR missingness | already coded (`--miss mar`), never run — **now in the grid** |
| Clade missingness | already coded (`--miss clade`), never run — **now in the grid** |
| Auxiliary covariate | `build_covs()` existed but was unreachable — **now wired and tested** (§6.1) |
| More replicates | **free** — the runner skips datasets that exist, so new seeds just append |
| Complete-case arm | did not exist — **now written and tested** (§6.1) |

All of it is written. What is left is agreeing the design and pressing go.

**What the move to n = 100 and 500 costs us in reuse.** Dan's decision (§3.0)
drops n = 300 and n = 1000, and n = 500 has never been run. Nothing at n = 300
or 1000 carries over, and only the n = 100 cells at λ ∈ {0.3, 0.7}, ρ ∈ {0, 0.5},
MCAR — 200 datasets, about 92 core-hours — can be reused at all. In practice the
new campaign starts from scratch. That is worth stating plainly: the 1,148
core-hours already spent are a sunk cost, not a head start.

---

## 3. What changes

### 3.0 Species counts: n = 100 and 500, nothing larger

**Decided (Dan, 2026-11-09).** The design is n = 100 and n = 500 only, at
**1000 replicates each**, with **all trait types scored at both**.

The reasoning is that these are the dataset sizes comparative biologists
actually work with. A 1000-species tree is not representative of the field, so
spending most of a compute budget on it buys reach into a regime few readers
occupy. Precision at realistic n is worth more than an asymptotic arm.

Two consequences to handle rather than argue with:

- **Two levels of n show a direction, not a curve.** With {100, 500} the paper
  can say accuracy and coverage do not degrade at larger trees; it cannot plot
  how they approach their limit. Only one sentence in the current ms cites
  n = 1000 by value (`ms.qmd` L336, the slope bias declining with n), so the
  rewrite is small — but the figures' x-axis becomes two points.
- **All traits at both n is nearly free.** Discrete scoring was previously
  capped at n ≤ 300. That was never a cost limit: the campaign record says
  "Cost: nothing extra in BACE; about 30 s per n = 300 dataset"
  (`records/2026-09-26-rubin-discrete-prerun.md` L44), and the reason it stopped
  was an approval hold — "n = 1000 on hold (Shinichi)" (same file, L3). At
  n = 500 it adds roughly a minute to a ~2 hour dataset, about 1%.

### 3.1 Drop λ = 1, use λ = 0.90 instead

**Every failure in the entire campaign came from λ = 1:**

| λ | datasets | failed | rate |
|---|---|---|---|
| 0.3 | 250 | 0 | 0% |
| 0.7 | 250 | 0 | 0% |
| **1.0** | 250 | **22** | **8.8%** |

and it is worst where there is least information to pin the variance down —
15% at n = 100, 5% at n = 300, 4% at n = 1000.

At λ = 1 the phylogenetic correlation matrix is pure Brownian motion and close
to singular; at λ = 0.90 it is 0.9V + 0.1I, which is comfortably invertible.
Moving to 0.90 should remove the failures outright and costs us nothing
scientifically — nobody's real data sits at exactly 1.

It also removes a weak spot in the paper. The ms currently has to say "failed
fits were excluded", and a 15% exclusion in one corner of the design is exactly
what a referee picks at.

**Cost:** none, as it happens. The change of species counts (§3.0) means the
whole campaign is re-run anyway, so dropping λ = 1 costs nothing extra — it is
simply one value in a loop. The ms text and figures that quote λ = 1 still need
updating.

One number worth keeping from the old campaign even though its datasets are
retired: the 22 failures were **all** at λ = 1, so if the re-run at 0.90 comes
back with zero failures we can say so, and drop the "failed fits were excluded"
caveat entirely.

### 3.2 The auxiliary covariate — correcting a wiring error in this plan

**Earlier drafts of this plan put Dan's three correlation levels on the wrong
axis.** This section replaces them.

`--rho` does *not* set one covariate's correlation. It fills a compound-symmetry
matrix — `Sigma_rho <- matrix(rho, K, K); diag(Sigma_rho) <- 1` in
`sim_latents()` — so with K = 8 it sets **all 21 trait-trait pairs to the same
value at once**. At ρ = 0.8 every missing trait would have six strong proxies.
BACE would score well, for a reason that tells a reader nothing about their own
data, and Claim 2 could never attribute the rescue to the covariate.

**What we actually want already exists.** `build_covs()` constructs exactly the
variable Dan described — think body size, measured for every species:

```
   cov  =  ρ_cov · (the target trait's latent)  +  √(1 − ρ_cov²) · e
```

one covariate, correlated exactly ρ_cov with **one** named trait, with the noise
term `e` itself evolving on the tree so the covariate is phylogenetically
structured rather than white noise.

- `ρ_cov` — the correlation we choose (0.8 for the good-covariate regime)
- `e` — phylogenetically structured noise, so the covariate looks like a real
  comparative trait rather than a random column

It is a predictor and never a target: the code keeps covariates outside `truth`,
so they are never masked and never scored. Positive-definiteness is automatic
for any ρ_cov between −1 and 1.

The authors recorded why the alternative fails — imposing covariate correlations
on `Sigma_rho` directly "is NOT positive definite at rho = 0 … min eigenvalue
−0.32 at r = 0.25".

**It has never been switched on.** `rubin_cell.R` calls `make_cell()` without
`n_cov`, which defaults to 0, so `build_covs()` returns `NULL`. Every result on
disk is covariate-free. The "+1 fully observed covariate" in older drafts of §2
was actually `d1`, the MAR *missingness driver*, hard-wired at 0.35 to all
traits — it causes the gaps, it is not the auxiliary.

**Three information regimes, not a factorial.** Reading Dan's original request
closely — *"correlation for an auxiliary (covariate) to be 0.8, others can be
0.5, and then a 0 correlation baseline"* — it describes three regimes:

| regime | ρ between traits | ρ of the covariate | what it represents |
|---|---|---|---|
| **floor** | 0 | — (no covariate) | nothing to borrow but the phylogeny |
| **typical** | 0.5 | — (no covariate) | related traits, no special predictor |
| **good covariate** | 0.5 | **0.8** | one well-measured variable, e.g. body size |

That is one axis with three levels, so **the grid stays at 27 scenarios per n
and the cost in §5.3 is unchanged**. It also makes Claim 2 a direct comparison:
*typical* versus *good covariate* at λ = 0.3 isolates what the covariate buys
when the phylogeny cannot help.

**Why it has to be one variable rather than all of them — arithmetic, not
preference.** For one variable correlated r with K traits that are themselves
mutually correlated ρ, the joint correlation matrix is a valid one only while

```
                 1 + (K − 1) · ρ
       r²   <   ─────────────────
                        K
```

At K = 7 traits and ρ = 0.5 that is r < 0.756. So "0.8 with everything" is not
a design we could run even if we wanted it — it is not a positive-definite
distribution. 0.8 with *one* trait is both runnable and the realistic case.

**A caveat to keep saying in the paper, now with numbers.** Two of them.

First, these correlations are on the latent Gaussian scale, before traits become
counts, categories and so on. The correlation you would measure on an observed
binary or ordinal trait is lower.

Second, and more surprising: the correlation **realised in any one dataset
scatters very widely**, because the covariate and the trait evolved on the same
tree, so the effective sample size is far below n. Measured over 300 datasets at
n = 120 (`tests-rubin/test-covariate.R` pins these):

| requested ρ_cov | mean realised | sd | what one dataset can look like |
|---|---|---|---|
| 0.8 | 0.799 | 0.096 | anywhere from about 0.60 to 0.95 |
| 0.3 | 0.314 | 0.222 | can come out negative |

The construction is unbiased — both means sit within about one standard error of
the target — but the paper should report the realised correlation as a
**distribution**, not as "0.8", or a reader will assume every dataset had a
covariate that good. This is itself a finding worth a sentence: *in a
comparative dataset, how useful your auxiliary variable turns out to be is
itself uncertain, and much more so than the sample size suggests.*

### 3.3 Add the two missingness mechanisms that already exist

- **MCAR** — gaps at random (what we have).
- **MAR** — gaps driven by a fully observed driver variable, which the code
  keeps correlated at 0.35 with the scored traits so that "MAR" really is MAR
  and not MCAR-on-an-unrelated-variable. Implemented:
  `P(miss) = plogis(a + log(3) × z(driver))`, with `a` solved so the expected
  fraction is 30%.
- **Clade** — whole clades go missing, because whole groups are understudied.
  Implemented by drawing clades of 5–15% of tips until the target fraction is
  reached, then trimming the overshoot.

**Clade is the one that matters most for this paper.** The traits have
phylogenetic signal, so relatives resemble each other and clades genuinely
differ. Removing a clade does not remove a random sample, it removes a
systematically different chunk — so complete-case should fail badly. It is also
the hardest case for BACE, because BACE's whole premise is borrowing from close
relatives and a missing clade is exactly where the close relatives are gone. If
BACE still recovers the truth there, that is a far stronger result than MCAR can
give, and it is specific to *phylogenetic* imputation rather than imputation in
general.

### 3.4 Add a complete-case arm

The only real new code. For each dataset, drop every species with any missing
trait, run the same downstream PGLS, and report the same estimands. No MCMC, so
it is cheap. Claim 1 lives or dies on this arm.

Worth noting what already exists and is *not* this: the campaign scores against
a `complete` arm, but that is the **complete-data oracle** (the truth before we
made holes), not complete-case. We need both — the oracle is the benchmark,
complete-case is the competitor.

### 3.5 What the complete-case arm will actually show — measured, before we run it

Worth stating now, because it changes how Claim 1 should be written. Complete-case
analysis can fail in two quite different ways, and the three mechanisms separate
them cleanly. Measured on the implemented code (`tests-rubin/test-completecase.R`
pins the MCAR figure):

| mechanism | species surviving | which failure this is |
|---|---|---|
| **MCAR** | **8%** (0.7⁷) | **attrition.** With 7 traits each 30% missing independently, a species survives only by dodging all seven draws. At n = 100 that is about 8 species, and often too few to fit anything at all. |
| **MAR** | ~14% | attrition, plus the survivors are a biased sample |
| **Clade** | ~68% | **bias.** Little attrition — the species outside the missing clades are untouched — but the survivors are systematically different. |

**The consequence for the paper.** Under MCAR, complete-case does not give a
*biased* answer, it gives *no usable answer*. That is a true and useful finding,
but it is a different sentence from "complete-case is biased", and a referee
will notice if we blur them. **Clade is where the bias story lives**, and that is
the mechanism the manuscript should lead on.

So the complete-case arm records, for every dataset, how many species survived,
and refuses to score at all below a floor of 10. The **refusal rate is itself a
headline result** under MCAR, not a diagnostic to be tidied away.

---

## 4. What we measure

Unchanged — these are the manuscript's definitions and we should not move them.

**zRMSE**, for continuous, count and proportion traits: the error in the
imputed values, in units of the trait's own spread.

```r
y_hat <- rowMeans(imputations)        # mean imputed value per cell
zrmse <- sqrt(mean((y_hat - truth)^2)) / sd(observed_values)
```

0 means perfect recovery; 1 means no better than filling in the average.
Proportions are clipped to [0.0001, 0.9999] and scored on the logit scale.

**Accuracy**, for binary, ordinal and categorical traits: how often the most
frequently imputed category was the true one, with **fractional credit for
ties** — if two categories tie and one is correct, that cell scores 0.5.

**Coverage**, for everything: how often the 95% interval or category set caught
the truth. Numeric intervals are the Rubin t-interval
`mean ± t(M−1, 0.975) × sqrt((1 + 1/M) × s²)`; discrete sets are the smallest
set of categories holding 95% of the imputed draws.

**Downstream**: bias and 95% coverage of the PGLS slope and the phylogenetic
correlation, pooled by Rubin's rules with Barnard–Rubin degrees of freedom —
now reported for **three arms**: oracle, complete-case, BACE.

We keep the idealised reference lines the ms already computes. Without them a
count zRMSE of 0.6 reads as poor when it may be the best achievable.

---

## 5. How big this gets, honestly

### 5.1 What a dataset actually costs — measured, not estimated

Earlier drafts of this plan quoted 22 min / 59 min / 4.1 h per dataset at
n = 100 / 300 / 1000. Those were an **apportionment**: a campaign total divided
up by a guessed scaling. They were wrong by a quarter.

Every one of the 750 result files records its own wall-clock, in
`walls[["bace_fit"]]`, written by `run_arm()` in `rubin_cell.R`. All 750 share
identical settings (nitt 50000, burnin 10000, thin 25, runs 5, M 20, arms bace)
on one host, so **n is the only thing that varies**. Reading them gives the real
distribution (728 clean runs; 22 crashed fits excluded):

| n | datasets | median | mean | IQR | max |
|---|---|---|---|---|---|
| 100 | 285 | **27.6 min** | 27.2 min | 26.3–29.0 | 31.8 min |
| 300 | 295 | **77.9 min** | 77.2 min | 76.0–79.9 | 107 min |
| 1000 | 148 | **4.48 h** | 4.31 h | 4.16–4.60 | 5.28 h |

So the old figures understated cost by 26% at n = 100, 32% at n = 300 and 9% at
n = 1000, and the campaign's true total is **1,148 core-hours** (the direct sum
over all 750 files; an intermediate draft said 1,177, which was medians times
counts rather than the actual sum). Cost is essentially **linear in n** — fitted
exponent 0.974 on datasets, 0.979 on means, 0.988 on medians — not sub-linear.

### 5.2 n = 500, which has never been run

Three independent routes, all on the current package, all agreeing:

1. **Matched short runs.** Identical cells at n = 100/300/500/1000, three seeds
   each, varying only n (coefficient of variation under 2%).
2. **A fixed-versus-iteration decomposition.** Each n run at two chain lengths,
   solving `t(n) = a(n) + b(n) × nitt` to separate one-off setup cost `a` from
   per-iteration cost `b`. Production uses long chains, so production scaling
   follows `b`, and this is what makes a short-chain extrapolation defensible
   instead of merely plausible.
3. **Anchoring** the measured `b(500)/b(n)` ratios on the production medians
   above.

| anchored on | implied n = 500 |
|---|---|
| production n = 100 | 2.14 h |
| production n = 300 | 2.15 h |
| production n = 1000 | 2.40 h |

**Plan on 2.2 h per dataset at n = 500**, with 2.1–2.4 h as the honest band.

### 5.3 What the chosen design costs

27 scenarios per n (3 λ × 3 ρ × 3 mechanisms) × 1000 replicates × 2 levels of n
= **54,000 datasets**.

| | datasets | core-hours | 150 cores | 200 cores |
|---|---|---|---|---|
| n = 100 @ 1000 reps | 27,000 | 12,400 | | |
| n = 500 @ 1000 reps | 27,000 | 59,400 | | |
| **total** | **54,000** | **71,800** | **~22 days** | **~16.5 days** |

Wall-clock includes ~10% for scheduling tails. Sensitivity to the n = 500
estimate is mild — at the pessimistic 2.4 h it is 17.7 days on 200 cores.

**n = 500 is 83% of the bill.** If the three weeks is too long, that is the only
place worth cutting, and replicates there are the cheapest thing to trim:
1000 at n = 100 with 250 at n = 500 is **6.3 days on 200 cores** and still gives
the n = 100 arm its full 0.7% precision, with 2.8% at n = 500 — enough for that
arm's job, which is to show the n = 100 result does not fall apart on bigger
trees. Recorded as an option, not a recommendation; the decision was 1000.

### 5.3b The largest unquantified term: this was all measured on one Mac

Every one of the 750 datasets ran on the same machine — `Daniel Noble Mac Pro`,
`aarch64-apple-darwin20`, checked across all 750 result files, no exceptions. The
150–200 cores we are proposing to use are **x86 cluster nodes**. Nothing in §5.1
or §5.2 says anything about how fast those are.

This is not a footnote on a three-week estimate. A plausible cluster factor of
0.7–1.3× turns "16 days on 200 cores" into **12 to 21 days**. The 2026-09
campaign did run on clusters (0.61–0.87× the Mac Pro per dataset) but with the
old package doing different work, so it is not a valid calibration.

Two further things the single-machine measurement cannot see:

- **Contention.** Every recorded time already includes seven MCMCglmm processes
  competing on one box. Contention at 150–200 cores on a cluster node is a
  different thing, in either direction.
- **Per-task startup.** R startup, library loading and tree simulation happen
  before the first log line, so they are not in `walls[["bace_fit"]]` at all.
  Seconds each — but paid 54,000 times on a SLURM array, plus scheduler latency.
  Budget a few percent.

**The fix is cheap and worth doing with the §6.4 pilot:** run 5–10 n = 500
datasets on the actual target cluster. That calibrates the whole estimate *and*
would be the only direct measurement of n = 500 in existence — every figure in
§5.2 is an interpolation.

### 5.4 Memory and disk are not constraints

Worth recording because it was raised and turned out to be a false alarm. The
cluster scripts request 4 GB per task at n = 100 and 10 GB at n = 300, which
suggested roughly 15 GB at n = 500 and therefore ~3 TB across 200 concurrent
jobs — infeasible. Those are **padded scheduler requests, not usage**. A
measured n = 500 cell peaks at about **0.6 GB**, so 200 concurrent jobs need
roughly **120 GB**. Not a constraint.

Disk is likewise small: existing result files average well under a megabyte, so
54,000 of them is tens of gigabytes. The one thing to watch is 54,000 files in a
single flat directory, which is awkward for the aggregation scripts — worth
splitting the output by n.

### 5.5 Precision, for the record

Coverage is a proportion, so its Monte Carlo error is `sqrt(0.95 × 0.05 / R)`:

| replicates | Monte Carlo SE | an observed 0.95 is consistent with |
|---|---|---|
| 50 | 0.031 | 0.89 – 1.00 |
| 250 | 0.014 | 0.92 – 0.98 |
| 1000 | 0.007 | 0.94 – 0.96 |

At 1000 replicates the paper can state coverage to better than a percentage
point, which is what Claim 3 needs. Morris, White & Crowther (2019,
*Statistics in Medicine* 38(11):2074–2102, doi:10.1002/sim.8086) set out this
calculation for choosing the number of repetitions; it is not yet in `bib/` and
should be if we make the argument in the ms.

### 5.6 Which machine

Both launchers exist already and encode their own limits:

- **`rubin_totoro.sh`** — Totoro, no queue, one cell per core, hard-capped at
  **150 cores** (D-143); it counts the user's busy cores and refuses to exceed
  it. Its memory guard (`free ≥ PAR × GB + 100`) is satisfied easily at 0.6 GB
  a cell.
- **`rubin_campaign_nibi.sh`** — Nibi, SLURM array on `def-snakagaw_cpu`, one
  core per task, concurrency set by `THROTTLE`. Nibi caps a user at 1,000
  submitted array tasks, so 54,000 datasets must go up in batches with many
  seeds per task; the script already computes the block size that fits and
  refuses an oversized submission. It also still hard-codes `rho in 0 0.5`, has
  no `--miss` flag, and by its own header requires Shinichi's sign-off.

---

## 6. The work, in order

### 6.1 Campaign code — written and tested (2026-11-09)

All done, `devtools`-style test suite green: **455 pass, 0 fail** across the
simulation lane, including 57 assertions in three new test files. The code lives
in `dev/rubin_study_2026-09-26/code/`; start at its `README.md`, which gives a
reading order, or run `Rscript code/walkthrough.R`, which takes one dataset
through every step and prints what it did at each one.

| # | what | where |
|---|---|---|
| 1 | **The auxiliary covariate reaches the model.** `--n_cov` / `--rho_cov` / `--cov_target` added; `fit_bace_mi()` now puts the covariate on the right-hand side of every trait's formula and never on the left. Verified: all 8 formulas carry `cov1`, none has it as a response. | `rubin_covariate.R` (new), `rubin_bace.R`, `rubin_cell.R` |
| 2 | **The complete-case arm.** Drops every species with a gap, runs the same two estimands with degrees of freedom from the *survivors*, records how many survived, and refuses below 10 rather than reporting a meaningless number. | `rubin_completecase.R` (new) |
| 3 | **The new grid.** n ∈ {100, 500}; λ ∈ {0.3, 0.7, 0.9}; the three information regimes; mechanisms `mcar mar clade`; output split into `n100/` and `n500/`. Dry-run confirms 54 scenarios. | `rubin_campaign_local.sh` |
| 4 | **Discrete scoring at both n**, no longer capped at n ≤ 300. | launcher, `rubin_discrete_truth.R` |
| 5 | **The completeness checker** no longer keeps its own copy of the design — it reads the launcher's task list, so it cannot drift from it again. | `rubin_missing.R` |
| 6 | **Reference lines** computed on the new grid. | `oracle_floor.R` |
| 7 | **Fewer files.** `code/` went from **32 tracked files to 19**, four of which are new (the README, the walkthrough, and the two new arms); 17 files of finished 2026-09 tooling (HTML report builders, the pre-run stage, one-off experiments) **moved** to `code/archive-2026-09/`. Nothing deleted. Tests went 13 → 16. | `code/README.md` (new) |

**Two silent failures found while doing this, both now impossible:**

**(a) The wrong BACE.** The repository was at **1.0.0** while the R library
still had **0.1.0.9000** installed, built two days earlier. Every run — including
the first rehearsals of this very plan — used the old package. Nothing warned
anybody; the results looked entirely normal. Both `rubin_cell.R` and the launcher
now refuse to start below a stated `BACE_MIN_VERSION`, and every result file
records the version that actually produced it.

**(b) Two scenarios writing the same filename.** The runner skips any dataset
whose `.rds` exists — that is what makes a three-week run restartable. But the
"typical" regime (ρ = 0.5, no covariate) and the "good covariate" regime
(ρ = 0.5, with one) share λ, ρ, n, M and seed. Without a covariate token in the
name they collide, and the second regime would be **skipped as already done and
never run** — no error, no warning, and a third of the design silently missing.
Filenames are now built in one place (`cell_tag()`), and a test asserts that all
54 scenarios get distinct names.

### 6.2 Manuscript code — this is the part that bites

**Changing n from {100, 300, 1000} to {100, 500} breaks `ms.qmd` silently.**
Not with an error — the figures would render, and would quietly contain no
data. These must all be fixed *before* the run finishes, or we will not notice:

| Where | What happens |
|---|---|
| `ms.qmd` L355 | `n_f = factor(n, c(100, 300, 1000))` — **every n = 500 row becomes `NA` and is dropped from all five panels of both figures.** ggplot only warns. This is the dangerous one. |
| `ms.qmd` L243 | `reps_big <- max(agg_fail$n_total[agg_fail$n == 1000])` → `max(numeric(0))` = `-Inf`; the figure caption then reads "… and -Inf at n = 1,000" |
| `ms.qmd` L254 | `bias_n1000 <- ...[sl_l1r5$n == 1000]` → `character(0)`; the inline value at L336 renders empty or errors |
| `ms.qmd` L382 | `subset(orc, ... & n <= 300)` — drops the n = 500 oracle ceiling lines |
| `ms.qmd` L242 | `reps_small <- max(...[n <= 300])` — would only ever see n = 100 |
| `oracle_floor.R` L35 | `NS <- c(100L, 300L, 1000L)` — the floor/ceiling reference lines would not be computed at n = 500 at all |

**Prose to update:** L287 ("n = 100, 300 or 1,000"), L271, L322 ("discrete traits
only at n ≤ 300" — no longer true), L334 (accuracy improving with n, now on two
points), L336 (the only sentence citing n = 1000 by value), L675 (anchors on
"~300 species"). L89, the abstract's "even a small number of species (~100)",
is unaffected.

### 6.3 Then run it

Re-run everything: as noted in §2, almost nothing carries over. Resumability is
already handled and proven — `rubin_cell.R` skips any dataset whose `.rds`
exists, so the run can be stopped and relaunched freely, and replicates can be
raised later without redoing earlier seeds.

---

### 6.4 One thing to settle before we commit three weeks: are the chains long enough?

**New information, dated after this plan was first drafted.** BACE 1.0.0
(2026-11-08) added `final_nitt`, `final_thin` and `final_burnin`, so the final
imputation runs can have different MCMC settings from the convergence phase. Its
own release notes explain why: *"The two phases do different jobs — the chained
phase only needs point imputations to settle, while the final runs supply the
draws that are pooled for inference — so they often want very different
chains."* The campaign currently uses one setting for both
(nitt 50,000, burn-in 10,000, thin 25 → 1,600 retained draws).

**Measured on the existing results, which used those settings.** The share of
fixed effects whose effective sample size is below 100:

| n | % of parameters under-mixed | BACE's own convergence verdict |
|---|---|---|
| **100** | **47%** | converged in **29%** of fits |
| 300 | 20% | 54% |
| 1000 | 3% | 81% |

That pools both phases, and the convergence phase genuinely does not need good
mixing. So we separated them on a fresh full-length n = 100, λ = 0.3 fit:

| phase | parameters | median ESS | % under 100 | worst |
|---|---|---|---|---|
| convergence | 95 | 1430 | 13% | 2.7 |
| **final — this is what feeds inference** | 504 | 1073 | **23%** | 2.9 |

So about a quarter of the parameters in the models that produce the published
numbers are poorly mixed at n = 100 — and n = 100 is half of this design.

**How much this matters is genuinely unclear**, and we should not overstate it.
Each final run contributes a *single* posterior-predictive draw, and the
imputation uncertainty comes from the spread *across* the 20 runs, not from
within any one chain. A badly mixed chain still gives a draw; it gives one from
a less well-converged posterior. The 2026-10 campaign ran at these settings and,
after the defects were fixed, performed well.

**But it is cheap to settle and expensive to be wrong about.** 20 datasets at
n = 100, run at the current settings and again with longer final chains, is
about **20 core-hours against a 72,000 core-hour campaign** — under 0.03%. If
coverage and zRMSE are unchanged, we launch with a defensible answer to the
referee who asks. If they are not, we found it before spending three weeks
rather than after. **Recommendation: run the pilot.** (Decision 3 in §8.)

---

## 7. What we expect to see

Written down before we run it, so we cannot talk ourselves into whatever
appears.

| | MCAR | MAR | Clade |
|---|---|---|---|
| Complete-case, slope | unbiased | mildly biased | biased |
| Complete-case, correlation | unbiased | **biased** | **badly biased** |
| BACE, slope | unbiased | unbiased | unbiased |
| BACE, correlation | unbiased | **unbiased** | **unbiased** |
| BACE coverage | ≈ 0.95 | ≈ 0.95 | ≈ 0.95, maybe a little low |
| BACE interval width | narrow | narrow | **wider** — correctly so |

The clade column is the one to watch. BACE's intervals *should* widen there,
because a missing clade genuinely is harder than missing singletons.
Wider-but-still-covering is the right answer. Narrow-and-still-covering would be
suspicious; narrow-and-not-covering would be a real problem, and far better
found now than after submission.

For Claim 2, across ρ × λ: BACE's zRMSE should worsen as λ falls **only when ρ
is low**. At ρ = 0.8 the λ = 0.3 result should be close to the λ = 0.9 result.
The ρ = 0, λ = 0.3 corner is the floor — genuinely no information — and BACE
should be honest about it with wide intervals rather than confident wrong
answers.

---

## 8. Things to decide

**Already settled** (Dan, 2026-11-09): 1000 replicates; λ = 0.90 replacing 1.0;
two levels of n (100 and 500); discrete scoring at both n; the three correlation
levels sit on the auxiliary covariate, not on all the traits at once.

### 1. What should the trait–trait correlation be in the two covariate regimes?

§3.2 proposes 0 for the floor regime and 0.5 for the other two, which is the
plain reading of Dan's *"others can be 0.5"*. **For Shinichi**, because it sets
how much work the covariate has to do: the higher it is, the less a good
covariate adds, and Claim 2 is measured against exactly that difference.

### 2. Is three weeks of cluster time acceptable?

~16.5 days on 200 cores, ~22 on 150 (§5.3). Dan has decided on 1000 replicates;
this is here so Shinichi can object before it starts rather than after. If it is
too long, §5.3 records the cheapest thing to trim (replicates at n = 500, which
is 83% of the bill).

### 3. Run the ~20 core-hour pilot on chain length first? — *we recommend yes*

§6.4. At n = 100, about a quarter of the fixed effects in the models that feed
inference have an effective sample size below 100. BACE 1.0.0 has just added the
means to fix that (`final_nitt`, `final_thin`). Checking costs under 0.03% of
the campaign and gives us a defensible answer either way. Launching without it
means finding out in three weeks' time, if at all.

### 4. Stop giving the fully observed driver its own imputation model? — *saves ~8.5%*

The MAR driver `d1` is complete by construction, so there is nothing in it to
impute, yet the campaign has always fitted a full model for it. Moving it to the
right-hand side of the other traits' models leaves those models **unchanged** —
`d1` is a predictor in them either way — for measured **8.5% less compute**
(14.7 s → 13.5 s on a matched n = 80 fit). On this campaign that is roughly
**6,000 core-hours, about 1.3 days**.

Implemented behind `--driver_as_predictor`, **off by default** so runs made
before 2026-11 still reproduce. Turning it on is a one-character change to the
launcher. It is a methodological tidy-up rather than a pure saving, so it should
be a decision rather than something that happens quietly.

### 5. A courtesy flag, not a question

Discrete scoring beyond n = 300 was previously recorded as "on hold (Shinichi)"
in the 2026-09 campaign notes. It now runs at both n, at about 1% extra cost.
Worth an explicit nod rather than a silent switch.
