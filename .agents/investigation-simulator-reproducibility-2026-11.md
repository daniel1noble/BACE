# The simulator gave different data on different computers

**Found:** 2026-11-12, by Shinichi, running the same scenario on three machines.
**Fixed:** 2026-11-13, in BACE 1.0.3.

---

## What went wrong

The same seed produced **different trait data on different computers**.

Not slightly different. Different:

| machine | `x1[1:3]` | realised signal |
|---|---|---|
| Mac (Accelerate) | 7.169, −1.485, 0.273 | 0.827 |
| Totoro (OpenBLAS) | 2.047, −1.537, 2.636 | 0.875 |
| Rorqual (AOCL) | 4.074, −1.417, 4.599 | 0.850 |

Each set was perfectly reasonable on its own. They were just not the same data,
so nobody could reproduce anybody else's results, and we could not tell a
harmless difference from a real bug.

## How we knew where to look

This is the useful part. On all three machines these were **identical**:

- the random numbers themselves (`runif`, `sample`)
- the family tree, down to its total branch length
- which cells were hidden

Only the trait values differed. So the random-number generator was fine, the
tree-building was fine, and the missingness was fine. Whatever went wrong had
to be in the one step between "random numbers" and "traits".

That narrowed a vague "different machines give different answers" into a single
function.

## Why it happened

To make related species resemble each other, the simulator takes a set of
independent random numbers and reshapes them using the family tree. Doing that
needs something like a **square root** of the tree's relationship matrix.

The catch: **a matrix has many valid square roots**, in the same way that √4 is
both +2 and −2. Both are correct. Pick either and the maths still works.

The old code used `MASS::mvrnorm`, which builds its square root with `eigen()`.
Eigenvectors have exactly this freedom — their signs are arbitrary, and nothing
in the standard says which one a maths library should return. Accelerate,
OpenBLAS and AOCL each made a different, equally valid choice.

So three machines took the same random numbers and reshaped them three
different ways. All correct. All different.

We confirmed it by reproducing the symptom on one machine: flipping the signs
of half the eigenvectors turned a value of −0.92 into +0.97, with the spread of
the data unchanged.

It also explains why an earlier small test looked fine. On a tiny 3x3 matrix
every library happens to agree; the freedom only shows up at realistic size.

## The fix

Use a **Cholesky factor** instead. For this kind of matrix there is only one,
so there is no choice for a maths library to make differently. Machines now
agree to rounding error rather than picking between alternatives.

Three things were checked before accepting it:

1. **The draws are still correct.** Over 4,000 draws they reproduce the
   requested variance and the tree's correlation structure.
2. **Imputation is no worse.** Accuracy against known truth is statistically
   identical before and after — mean correlation 0.510 against 0.528 over 30
   seeds, p = 0.85. An early 8-seed check suggested a drop; that was noise, and
   acting on it would have meant reverting a correct fix.
3. **Nothing downstream moved.** The new draw uses exactly the same number of
   random values, so trees, tip labels and missingness patterns are unchanged.
   Only the traits differ, which is the point.

## What this changed

**Every simulated dataset.** Unavoidable: the old values were machine-dependent,
so there was nothing worth preserving.

**The version, to 1.0.3.** This matters. A fixed build would otherwise still
have called itself 1.0.2, and the campaign's "needs at least 1.0.2" guard would
have accepted the broken one. Same trap as 1.0.1 → 1.0.2. The guards in
`run_one.R` and `run_all.sh` now require 1.0.3.

## How we stop it coming back

`tests/testthat/test-simulator-reproducibility.R` includes a **golden test**: a
fixed tree and seed with the expected numbers written down, at a tolerance of
1e-8 — loose enough for rounding, tight enough to catch a sign flip.

That means **running the test suite on any machine is a cross-machine check**.
If it passes on Totoro and on Rorqual, the simulator is machine-independent
there. The test says in a comment not to "fix" a failure by widening the
tolerance or regenerating the numbers.

The suite also pins the *method*, not just the answer, so reverting to
`mvrnorm` fails even if the distribution is still right.

## Two things still open

**This fixes the simulator, not the imputation.** Whether MCMCglmm also
diverges across machines is a separate question. It is now testable in
isolation, because the input data will finally match.

**Package versions are still unpinned.** There is no `renv`, so two machines can
run different MCMCglmm or ape versions. Each result file now records them, but
recording is not preventing.

## Related

- [investigation-audit-2026-11.md](investigation-audit-2026-11.md) — the
  numerical audits of the package itself
- `dev/simulation_2026-11/reproducibility_probe.R` — run on any machine,
  no MCMC, seconds
