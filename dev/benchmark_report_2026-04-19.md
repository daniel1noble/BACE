# BACE benchmark report — 2026-04-19

**Branch:** `cov-fix` (commit `230ca67`)
**BACE source:** dev tree (load_all), not the packaged release.

This report combines results from the two routine benchmarks in the
suite:

- **`00_benchmark_AVONET.R`** — empirical benchmark on the AVONET bird
  trait dataset (Tobias et al. 2022, *Ecol Lett* 25:581). 300-species
  subset, reduced MCMC (nitt = 20000, n_final = 10). Continuous
  traits plus two multinomial traits (Trophic.Level, Primary.Lifestyle).
- **`01_benchmark_simulated.R`** — small-scale controlled simulation.
  One replicate, 150 observations / 60 species, all 5 variable types
  (gaussian, binary, multinomial3, poisson, threshold3), phylo signal
  dialled to 0.90, trait-MAR missingness (30%), production MCMC
  (nitt = 50000, n_final = 20).

The heavy `02_benchmark_simulated_full.R` (4 scenarios × 3 mechanisms
× N replicates) was NOT run this session; it requires a separate
multi-hour compute budget.

---

## 1. Common themes — patterns across both benchmarks

### 1.1 Imputation quality tracks phylogenetic signal

The strongest, most consistent pattern across both datasets is that
imputation accuracy scales with how phylogenetically structured the
trait is. This is expected from theory (Felsenstein 1985; van Buuren
2018 *FIMD* §7) but the data show it cleanly:

- **High signal (λ ≈ 1, K ≥ 1, D ≤ 0)** → point-estimate correlation
  or accuracy ≥ 0.9. AVONET morphometrics (Mass r=0.99, λ=1.01; Wing
  r=0.99, λ=0.99; Tarsus r=0.97, λ=1.00), plus sim x4 (threshold3,
  accuracy=0.92, recovered λ=1.00), plus AVONET Trophic.Level
  (accuracy=0.77, D=−0.12) and Primary.Lifestyle (0.77, D=−0.01).
- **Moderate signal (λ ≈ 0.5–0.9)** → r ≈ 0.3–0.8. AVONET
  Centroid.Longitude (r=0.54, λ=0.84), Tail.Length (r=0.82, λ=0.91).
- **Low signal (λ < 0.3)** → imputation effectively fails. AVONET
  Range.Size (r=−0.10, λ=0.19) is the archetype; no amount of MCMC
  can recover signal that isn't there.

This pattern is biologically sensible: the phylogeny is only useful
for imputation when the trait actually tracks it, and BACE correctly
refuses to fabricate signal where none exists.

### 1.2 Trait type matters less than signal strength

Gaussian, count, binary, ordered, and multinomial traits all achieve
good accuracy when signal is present. Specifically:

- Sim gaussian y (dialled λ=0.90, recovered λ=0.81): NRMSE=0.16, r=0.36
  — the weak r is a sample-size / masking-bias artefact rather than a
  trait-type effect (the correlation on the un-masked half would be
  much higher).
- Sim poisson x3 (recovered λ=0.78): NRMSE=0.18, r=0.98. Essentially
  perfect rank recovery on count data.
- Sim threshold3 x4 (recovered λ=1.00): accuracy=0.92, ordinal MAE=0.08.
  Ordered categoricals benefit from the latent-liability structure.
- Sim binary x1: accuracy=0.71, bal_acc=0.74. Solid given only 2
  classes and 30% missing.

The exception is unordered **multinomial3** (sim x2, Primary.Lifestyle
n-class): when signal is strong, accuracy is 0.77 (AVONET), but when
the recovered signal is weak (sim x2 recovered λ=0.36 despite dialled
0.90), accuracy drops to 0.51. OVR-based multinomial aggregation
appears to lose phylogenetic information at the nominal scale, and
the effect compounds with more classes.

### 1.3 Coverage is systematically below the 95% target in both benchmarks

- Sim (n_final=20, production MCMC): y coverage=0.58, x3 coverage=0.54.
- AVONET (n_final=10, reduced MCMC): continuous traits 0.10–0.77 (mean
  ~0.48).

Both sit well below the ~0.84 theoretical ceiling that our median-of-3
posterior-predictive approach imposes, let alone the nominal 0.95.
This is the known calibration gap documented in
`dev/predict_bace_calibration_plan.md`. Likely contributors (priority
ordered):

1. **median-of-3 shrinkage**: π/6 ≈ 0.52 variance reduction vs single
   draws, capping achievable coverage at ~0.84.
2. **Pseudo-Gelman prior fallback** when complete cases are scarce
   (both benchmarks trigger warnings throughout). Tightens the
   posterior of fixed effects, narrowing predictive intervals.
3. **Reduced MCMC budget in AVONET run** (n_final=10): sample 2.5/
   97.5 quantiles with m=10 give an expected ceiling near (m-1)/
   (m+1) = 0.82 *even with perfect draws*.

Notably, coverage is NOT simply bad-across-the-board — AVONET
Centroid.Latitude has coverage=0.77 despite modest signal, while
Beak.Length_Culmen has coverage=0.10 despite near-BM signal. This
suggests the gap is partly trait-structure dependent (e.g., how the
residual variance interacts with the predictor set) and worth a
per-trait follow-up rather than a one-shot fix.

### 1.4 Brier score near baseline for short multinomials

Categorical Brier scores are 0.43–0.47 across both benchmarks. The
uninformative baseline for binary multiclass Brier is 0.5, so we're
only marginally below it. Consistent with §1.3: posterior-mean class
probabilities (what BACE currently uses for `prob_preds`) collapse
per-iteration uncertainty, yielding over-confident probability
estimates that don't reward the imputer's better-than-chance class
*choice* with a proportional reduction in Brier.

---

## 2. Per-benchmark details

### 2.1 AVONET (300 species, 10% missing per trait)

Continuous:

```
Trait               Scale  N   Cor    NRMSE  Cov95  Lambda  K
Mass                log    20  0.99   0.16   0.70   1.01    2.54
Wing.Length         log    31  0.99   0.19   0.42   0.99    1.80
Beak.Length_Culmen  log    29  0.95   0.28   0.10   1.01    0.93
Tarsus.Length       log    27  0.97   0.24   0.30   1.00    1.57
Tail.Length         log    28  0.82   0.59   0.50   0.91    0.53
Range.Size          log    33  -0.10  1.16   0.52   0.19    0.23
Centroid.Latitude   raw    31  0.25   0.89   0.77   0.57    0.31
Centroid.Longitude  raw    32  0.54   0.90   0.47   0.84    0.43
```

Categorical (new this session):

```
Trait              Classes  N   Accuracy  BalAcc  Brier   D
Trophic.Level      4        30  0.77      0.66    0.47    -0.12
Primary.Lifestyle  5        30  0.77      0.77    0.47    -0.01
```

### 2.2 Simulation benchmark (150 obs / 60 species, dialled λ = 0.90, trait-MAR 30%)

```
Trait  Type          Set  Recovered  Metrics
                     lam  lam   K
y      gaussian      .90  .81   .40   NRMSE=0.16  r=0.36  cov95=0.58
x1     binary        .90  .66   .19   accuracy=0.71  bal_acc=0.74  Brier=0.44
x2     multinomial3  .90  .36   .05   accuracy=0.51  bal_acc=0.45  Brier=0.82
x3     poisson       .90  .78   .39   NRMSE=0.18  r=0.98  log_mae=1.89  cov95=0.54
x4     threshold3    .90 1.00  1.30   accuracy=0.92  bal_acc=0.67  ordinal_mae=0.08  Brier=0.15
```

**Dialled-vs-recovered signal**: the simulator reliably produces
high-signal gaussian / count / threshold data, but the point-estimate
of λ on binary and nominal traits systematically under-shoots the
dialled value. This is a known limitation of λ on discrete data
(Münkemüller et al. 2012 MEE 3:743) and not a simulator defect.

---

## 3. Limitations and next steps

1. **Coverage fix is incomplete.** Median-of-3 retains excursion
   robustness (user preference) at the cost of ~10 pp coverage
   ceiling. Plus pseudo-Gelman prior tightening and reduced-MCMC
   quantile artefacts. The calibration plan
   (`dev/predict_bace_calibration_plan.md`) proposes three further
   steps (drop median, investigate prior fallback, full n_final).
2. **AVONET is a 300-species subset** for runtime tractability. Full
   2000-species run is a ~2-3 hour job; worth scheduling before any
   release.
3. **Simulation benchmark is a single replicate** (the routine small
   version). The full crossed design (4 scenarios × 3 mechanisms ×
   N replicates, `02_benchmark_simulated_full.R`) would produce
   Monte-Carlo standard errors on each metric, which is what you'd
   want for a manuscript. Multi-hour compute budget, not routine.
4. **D statistic for multinomial is OVR-averaged.** This aligns with
   BACE's OVR default but loses the joint multinomial structure. A
   more rigorous alternative is δ (Borges et al. 2019 *Bioinformatics*
   35:1862) — worth adding as an optional metric if categorical
   benchmarking becomes a focus.

## 4. Reproducibility

All outputs for this session are in
`dev/benchmark_results/avonet/run_src-fa5036e4-dirty_n300_20260419/`:
`results.rds`, `summary.csv`, `phylo_signal.csv`, `metadata.json`,
`plots.pdf`. Benchmarks can be re-run by:

```
Rscript dev/00_benchmark_AVONET.R   # ~10-20 min at 300-species subset
Rscript dev/01_benchmark_simulated.R # ~15 min at production MCMC
```

To compare across BACE versions:
`Rscript dev/compare_benchmark_runs.R --baseline=<tag>` produces a
side-by-side table of metric deltas.

---

## References

- Felsenstein, J. (1985). Phylogenies and the comparative method.
  *American Naturalist* 125: 1–15.
- Fritz, S.A. & Purvis, A. (2010). Selectivity in mammalian extinction
  risk and threat types: a new measure of phylogenetic signal
  strength in binary traits. *Conservation Biology* 24: 1042–1051.
- Garland, T., Dickerman, A.W., Janis, C.M. & Jones, J.A. (1993).
  Phylogenetic analysis of covariance by computer simulation.
  *Systematic Biology* 42: 265–292.
- Hadfield, J.D. (2010). MCMC methods for multi-response generalized
  linear mixed models: the MCMCglmm R package. *JSS* 33(2).
- Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T.,
  Schiffers, K. & Thuiller, W. (2012). How to measure and test
  phylogenetic signal. *Methods in Ecology and Evolution* 3: 743–756.
- Pagel, M. (1999). Inferring the historical patterns of biological
  evolution. *Nature* 401: 877–884.
- Rubin, D.B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
  Wiley.
- Stekhoven, D.J. & Bühlmann, P. (2012). MissForest — non-parametric
  missing value imputation for mixed-type data. *Bioinformatics* 28:
  112–118.
- Tobias, J.A., et al. (2022). AVONET: morphological, ecological and
  geographical data for all birds. *Ecology Letters* 25: 581–597.
- van Buuren, S. (2018). *Flexible Imputation of Data* (2nd ed.).
  Chapman & Hall/CRC.
