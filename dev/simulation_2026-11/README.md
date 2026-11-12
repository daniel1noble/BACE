# The simulation

Make comparative datasets where we know the true answer. Hide some values. Ask BACE to
fill them back in. Check whether it got them right, and whether it was honest about how
sure it was.

Full reasoning and the measured costs behind every number here:
[`.agents/plans/simulation-plan-2026-11.md`](../../.agents/plans/simulation-plan-2026-11.md).

---

## What we want to test

**1. BACE beats throwing data away.**
What people normally do with gaps is delete every species that has one. That is fine if
the gaps are purely random. It is not fine when the gaps depend on something — say,
poorly studied species are both more likely to be missing a trait *and* genuinely
different in it. We want to show that deleting gives the wrong answer in that case and
BACE does not.

**2. A good covariate rescues you when the family tree cannot.**
BACE borrows information from two places: related species, and other traits on the same
species. When species don't resemble their relatives much, the first source dries up. We
want to show that isn't fatal, as long as you have another trait that tracks the missing
one. This is the practical advice the paper can give: what to put in your imputation
model.

**3. BACE is honest about every kind of trait.**
Little bias, and its "95% sure" ranges really do contain the truth about 95% of the time
— for continuous, count, binary, ordinal and categorical traits alike.

---

## What will actually run

Four things are varied:

| | values | why |
|---|---|---|
| species | 100, 500 | does the answer hold on bigger trees? |
| how much relatives resemble each other | 0.3, 0.7, 0.9 | weak, moderate, strong |
| what else there is to borrow from | three settings, below | this is Claim 2 |
| how the gaps are made | three ways, below | this is Claim 1 |

**What else there is to borrow from:**

| name | meaning |
|---|---|
| `floor` | traits unrelated to each other, no helper variable — nothing to go on but the tree |
| `typical` | traits related to each other — the ordinary situation |
| `good_covariate` | traits related, *plus* one fully observed helper (think body size) |

**How the gaps are made:**

| name | meaning |
|---|---|
| `mcar` | purely at random |
| `mar` | driven by another variable we can see |
| `clade` | whole related groups go missing — the realistic bad case |

30% of values are hidden in each trait.

That is **54 scenarios**, each run **1000 times** with different random data =
**54,000 datasets**.

Per dataset, at the campaign's real settings (M = 20, nitt = 50,000): **20-30 minutes**
at 100 species, **about 2 hours** at 500. That totals **52,000-72,000 core-hours**, or
roughly **12-17 days on 200 cores**. The range is machine speed, not uncertainty about
the design: the plan was costed on one Mac at 27.6 min per n = 100 dataset, and a faster
box measured 18-23 min, which is where the lower figure comes from. Time one dataset on
your own machine before committing to a schedule.

**The 500-species runs are 83% of that bill.** If three weeks is too long, that is the
only place worth cutting — and cutting replicates there is the cheapest trim. 1000 runs
at 100 species plus 250 at 500 species comes to about **6 days on 200 cores** and still
answers the question that arm is there to answer.

---

## What it is compared against

| | what it is | what it's for |
|---|---|---|
| `truth` | the complete data, before anything was hidden | the best anyone could do — **and the sanity check**: if its answer isn't roughly unbiased with ~95% coverage, something is wrong with the setup |
| `complete_case` | delete every species with a gap | what people do now. Claim 1 is this versus BACE |
| `bace` | fill the gaps with `bace()` | the method the paper is about |

---

## What gets measured

**Did it get the values right?**
For continuous and count traits, the error in the filled-in values, measured in units of
that trait's own spread. 0 is perfect; 1 is no better than using the average.
For categories, how often the most-imputed category was the right one.

**Was it honest about how sure it was?**
How often the "95% sure" range actually contained the truth. Should be 95%.

**Does the downstream answer come out right?**
The relationship between two traits, fitted on the filled-in data, compared against the
truth — for all three arms.

---

## What we expect to see

Written down before running, so we can't talk ourselves into whatever appears.

| | random gaps | driven gaps | whole groups missing |
|---|---|---|---|
| deleting — relationship | fine | mildly wrong | wrong |
| deleting — correlation | fine | **wrong** | **badly wrong** |
| BACE — relationship | fine | fine | fine |
| BACE — correlation | fine | **fine** | **fine** |
| BACE — honesty | ~95% | ~95% | ~95%, maybe slightly low |
| BACE — range width | narrow | narrow | **wider** — correctly so |

**The last column is the one to watch.** BACE's ranges *should* get wider when whole
groups go missing, because that genuinely is a harder problem. Wider but still catching
the truth is the right answer. Narrow and still catching it would be suspicious. Narrow
and missing it would be a real problem — and far better found now than after submission.

For Claim 2: BACE should only get worse as relatives stop resembling each other **when
there's nothing else to borrow from**. With a good helper variable, weak phylogenetic
signal shouldn't hurt much.

---

## Running it

```bash
# 1. Install the package being measured. The scripts refuse to start on an older one.
cd <BACE repo root> && R CMD INSTALL .

# 2. One dataset, to watch it work end to end (20-30 min at n = 100,
#    about 2 h at n = 500, at the campaign's real settings)
Rscript run_one.R --id n100_ps0.7_good_covariate_clade --seed 1 --out results

# 3. See what the whole campaign would do, without doing it
DRY_RUN=1 bash run_all.sh        # writes results/tasks.txt — every command, readable

# 4. Run it. Stop and restart whenever you like.
CORES=8 bash run_all.sh

# 5. Turn the results into tables
Rscript collect.R results
```

### The code, one file per step

| file | step |
|---|---|
| `R/01_scenarios.R` | the list of scenarios — print `make_scenarios()` to see exactly what will run |
| `R/02_simulate.R` | one dataset: traits, tree, and the three ways of making gaps |
| `R/03_fit.R` | run BACE, and the delete-the-gaps comparison |
| `R/04_score.R` | the scoring |
| `run_one.R` | one dataset end to end — the unit of work |
| `run_all.sh` | runs `run_one.R` once per dataset |
| `collect.R` | gathers everything into five CSVs |

---

## Things that will bite

**The wrong version of BACE.** An old copy sitting in the R library produces results that
look completely normal and mean nothing. This has already happened once, on 2026-11-09.
Both `run_one.R` and `run_all.sh` refuse to start below **1.0.3**, and every result file
records the version that made it. If you change the package, reinstall it *and* bump the
version — the guard can only catch what the version number tells it.

**Interruptions.** A wall-time limit or a dead node can hit any dataset, and at 500
species that's 2.2 hours lost. BACE checkpoints by default, so a resumed run costs about
5% of the time and gives an identical answer. They're deleted as soon as a dataset
succeeds — only interrupted ones leave anything behind, which is exactly what should be
resumed. `--no_checkpoint` turns it off.

They are **big**, and the size grows in direct proportion to the number of species:

| | per running dataset | × 200 concurrent |
|---|---|---|
| n = 100 | 1.5 GB | 300 GB |
| n = 500 | 7.5 GB | **~1.5 TB** |

It does not accumulate, but that much scratch space has to exist while the run is going.
The scenarios are sorted by tree size, so the two halves run one after the other: the
300 GB and the 1.5 TB are separate peaks, not a sum, and **the 1.5 TB belongs entirely to
the n = 500 half**.

If 1.5 TB is awkward, the choice is between disk and time:

| | disk at peak | n = 500 phase | cost |
|---|---|---|---|
| 200 at once, checkpointed | 1.46 TB | 11.2 days | baseline |
| 100 at once, checkpointed | 0.73 TB | 22.5 days | **+11 days** |
| 200 at once, no checkpoints | 0 | 11.2 days | re-run whatever is interrupted |

Dropping checkpoints looks like the worse idea and is usually the better one.
Checkpointing is insurance: you pay disk on *every* run to avoid a loss that only happens
to *some*. At 500 species a lost dataset costs one 2-hour re-run, so even a 20 per cent
interruption rate across the whole half is about 10,800 core-hours — roughly 2 days at 200
cores, against the 11 days that halving concurrency costs.

**The exception, and it is the reason checkpointing exists: a wall-time limit shorter than
a single dataset.** If jobs are killed at, say, 1 hour, then *every* n = 500 run is
interrupted, the rate is 100 per cent rather than 5, and checkpointing is the only thing
that lets the half finish at all. Check the queue's wall-time limit before deciding.

**The R version used to change the data.** `sample()` changed in R 3.6.0, so before
2026-11-12 the same seed produced *different* trees and hid *different* cells on different
R versions — not slightly different, unrelated. Measured on seed 1: tree length 230.4 vs
182.0, ordinal counts 2/36/51/11 vs 4/31/49/16. `run_one.R` now pins the generator
explicitly, and every result file records the R version, platform and package versions. If
you are comparing results between machines, run `reproducibility_probe.R` on both first —
it takes seconds, runs no MCMC, and tells you whether you are even analysing the same data.
Note the pin fixes the *data* for R >= 3.6; it does not make results identical across
different MCMCglmm or ape versions, which is why those are recorded.

**Stopping and restarting is always safe.** Any dataset whose result file already exists
is skipped, so you can interrupt the campaign freely, and add more replicates later
without redoing earlier ones.

**Failures don't lose data.** If any stage fails, the failure is recorded and the dataset
still writes its file. `collect.R` counts them into `failures.csv`.

---

## Two numbers that are easy to misread

**Phylogenetic signal** is set on the hidden scale the models actually work on, not on the
values you see. The two are the same thing for a continuous trait and very different for a
categorical one, so "asked 0.7, measured 0.17" is not a failure. Measured over 25 seeds,
all asked for 0.7:

| trait | mean realised lambda | how often exactly 0.00 |
|---|---|---|
| `x1` continuous | 0.68 | 0% |
| `cnt` count | 0.45 | 0% |
| `ord` ordinal | 0.54 | 4% |
| `bin` binary | **0.17** | **16%** |

A binary trait reading 0.00 is ordinary; the same from `x1` is not, and would be worth
chasing. Every dataset records what it actually realised.

**`trait_beta` and `aux_beta` are regression slopes, not correlations**, and a slope buys
less correlation than its own value: slope 0.5 gave a correlation of 0.45, slope 0.8 gave
0.62. Every dataset records its own realised correlation, so the paper reports what the
data had rather than what was asked for.

---

## Open — needs a decision, not a code fix

1. **BACE's ranges are drawn slightly too narrow.** The plan says to use the
   Barnard–Rubin small-sample adjustment; the code doesn't apply it (`pool_mi()` is
   called without `df_fun`). Costs about half a percentage point of coverage at 100
   species, less at 500, and it makes BACE look slightly *worse* than it is. Small enough
   to note in the paper rather than fix.
2. **Counts come out about 8% low when the big values are the missing ones.** BACE caps
   imputed counts to stop runaway values, and that cap is set from the counts it can see
   — which are the small ones. It warns when this happens from 1.0.2, but the cap itself
   is unchanged. Affects the count trait under driven missingness.
3. **The delete-the-gaps arm deletes on all six traits**, though the analysis only uses
   two. That's harsher than any real analyst would be, and the penalty depends on the
   mechanism, so some of the Claim 1 pattern comes from the comparison rule rather than
   from the missingness. A deliberate choice — worth revisiting before the paper, not
   before the run.
4. **Unconfirmed:** if a category trait loses a whole level (possible when whole groups go
   missing), imputed values may land on the wrong category label. One reviewer reproduced
   it, another couldn't. Worth testing before committing the cluster time.
