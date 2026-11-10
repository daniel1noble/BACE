# Handover note for Shinichi — the simulation campaign

Draft for Dan to send. Written 2026-11-10.

---

Hi Shinichi,

The simulation is ready to run. Everything is on `main` at
https://github.com/daniel1noble/BACE, in the folder `dev/simulation_2026-11/`.

## What to read

**Start here:** `dev/simulation_2026-11/README.md` — about five minutes. It says what
the simulation does, what it will run, what gets measured, and what we expect to see
(written down before running, so we can't talk ourselves into whatever appears).

**If you want the reasoning behind any number:**
`.agents/plans/simulation-plan-2026-11.md` — the full design document, including how
every cost estimate was measured rather than guessed.

**If you want to know what we fixed first:**
`.agents/investigation-audit-2026-11.md`. We audited the package before running
anything and found ten bugs that gave wrong numbers silently — no error, no warning.
Two would have changed published results. All are fixed, each with a test that checks
the number rather than just that the code runs.

## How to run it

The package must be installed first, and it must be **1.0.2 or newer**. The scripts
refuse to start on anything older, because an old copy produces results that look
completely normal and mean nothing — that happened to us once already.

```bash
cd <BACE repo root>
R CMD INSTALL .

cd dev/simulation_2026-11

# 1. One dataset, end to end. About 2 minutes. Watch it work.
Rscript run_one.R --id n100_ps0.7_good_covariate_clade --seed 1 --out results

# 2. See exactly what the full campaign would do, without doing it.
DRY_RUN=1 bash run_all.sh          # writes results/tasks.txt — every command, readable

# 3. A small rehearsal of the identical thing (2 replicates instead of 1000).
REPS=2 CORES=8 bash run_all.sh

# 4. The real thing.
CORES=200 bash run_all.sh

# 5. Turn results into tables.
Rscript collect.R results
```

Stop and restart it whenever you like — any dataset already done is skipped. You can
also raise the number of replicates later without redoing earlier ones.

**Size:** 54 scenarios × 1000 replicates = 54,000 datasets, about 72,000 core-hours,
roughly **16–17 days on 200 cores**. The 500-species runs are 83% of that. If that's too
long, cutting replicates there is the cheapest trim: 1000 at n=100 plus 250 at n=500 is
about **6 days on 200 cores** and still answers the question that arm exists to answer.

## What comes out

Five CSVs in `results/`: per-dataset slopes and recovery, three summary tables, and
`failures.csv` counting anything that went wrong. Nothing is ever silently dropped — a
failure is recorded and the dataset still writes its file.

## Four things we'd like your view on

None of these are bugs. They're choices, and they affect what the paper can claim.

**1. The ordinal trait collapses when phylogenetic signal is low.** It is specified as
four categories. At signal 0.3 it comes out with three, one of them holding a single
species. So the "does BACE handle ordinal traits" part of Claim 3 is, in those
scenarios, being tested on something closer to a binary trait. This is the data
generator, not BACE. We think it's worth fixing before the run — it weakens a third of
Claim 3 in a third of the scenarios.

**2. BACE's intervals are drawn very slightly too narrow.** The plan specifies the
Barnard–Rubin small-sample adjustment; the code doesn't apply it. It costs about half a
percentage point of coverage at n=100 and less at n=500, and it makes BACE look slightly
*worse* than it is, never better. Our inclination is to note it in the paper rather than
change anything.

**3. Count imputations come out about 8% low when the large counts are the missing
ones.** BACE caps imputed counts to stop runaway values, and the cap is set from the
counts it can see — which, under that kind of missingness, are the small ones. From
1.0.2 it warns when the cap binds, but the cap itself is unchanged. This affects the
count trait under driven missingness.

**4. The complete-case comparison deletes species using all six traits**, although the
analysis model uses only two. That is harsher than any real analyst would be, and the
penalty depends on the missingness mechanism — so some of the Claim 1 pattern comes from
the comparison rule rather than from the missingness itself. A deliberate choice, but
worth revisiting before we write it up.

Happy to talk any of these through.

Dan
