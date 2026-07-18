# pigauto — the companion package, and how BACE differs

**pigauto** (github.com/itchyshin/pigauto; "itchyshin" = Shinichi Nakagawa, a BACE
co-author) is BACE's **companion package**, co-designed to divide the phylogenetic-
imputation problem. They cite and interoperate. The manuscript should position them
as complementary, not rivals.

## The division of labour (both packages already encode it)

- pigauto's `pool_mi()` **rejects MCMCglmm fits** and points users to
  `BACE::bace() + BACE::pool_posteriors()` for the integrated Bayesian workflow.
- pigauto's `fit_baseline_bace()` calls `BACE::bace_imp()` as one of its baselines.
- pigauto keeps a full copy of BACE in-tree as a comparison baseline, and adopted
  BACE's **OVR one-vs-rest categorical trick** (credited: lifted their AVONET
  Trophic.Level accuracy ~42%→72%).

## How pigauto works (for contrast)

Torch-based, **no MCMC**. Core formula per cell:
`pred = (1 − r_cal)·baseline + r_cal·GNN_delta`, plus a linear covariate term.

- **Engine**: a graph-transformer denoising autoencoder (`ResidualPhyloDAE`). `delta`
  predicts the truth directly (not a residual). Phylogeny enters as attention bias
  (rate-aware Gaussian on cophenetic distances) + Laplacian spectral node features.
- **Baseline** (fully in-house, no Rphylopars): Brownian-motion conditional MVN for
  continuous/count/ordinal (GLS mean, REML variance, `SE = √(σ²(1−h_i))`); phylo label
  propagation / truncated-Gaussian liability threshold models for binary/categorical;
  **Henderson's sparse precision** (Hadfield & Nakagawa 2010) for the joint MVN.
- **r_cal**: per-trait gate calibrated on a held-out split; closes to 0 when the
  baseline is already optimal (so the GNN is a no-op on high-signal traits).
- **Uncertainty**: split-conformal intervals (distribution-free ≥95% coverage) — the
  primary CI; MI draws via conformal jitter / MC-dropout / mice-style PMM.
- Scales to 10k species; 8 trait types incl. proportion/zi_count/compositional;
  active-learning `suggest_next_observation()`.

## BACE vs pigauto — the substantive differences

| Axis | BACE | pigauto |
|---|---|---|
| Engine | fully Bayesian MCMC (MCMCglmm), single engine for imputation + inference | GNN + classical phylo baseline, deterministic ML/REML (no MCMC) |
| Imputations | proper posterior-predictive draws (Rubin-"proper", gold standard) | point estimate + conformal jitter / PMM (proper via calibration, not Bayesian coherence) |
| Multivariate coupling | chained equations condition every variable on all others | joint cross-trait EM off by default → per-column BM + one closed-form Σ |
| Uncertainty guarantee | posterior coverage (relies on model being right) | conformal (distribution-free guarantee) |
| Scale | MCMC per variable per iteration; smaller | tested to 10,000 species |
| Pooling | stacked posterior (mixture) `pool_posteriors`; Rubin via `pool_mi` (added) | Rubin's rules `pool_mi` (rejects MCMCglmm → points to BACE) |

pigauto even concedes the Bayesian point: its `lambda="bayes"` option is documented as
"the deterministic analogue of what BACE's MCMCglmm does … a fast numerical Riemann sum
rather than a chain."

## Manuscript framing

BACE = principled Bayesian, single-engine, best uncertainty calibration, smaller scale.
pigauto = flexible/scalable (GNN, 10k species), conformal intervals, multi-tree Rubin
pooling (Nakagawa & de Villemereuil 2019). **Complementary companion packages** with a
clean division of labour — cite each other.
