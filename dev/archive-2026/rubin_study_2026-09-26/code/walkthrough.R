# =================================================================================================
# walkthrough.R -- one dataset, every step, so you can check the simulation yourself.
#
# Best read by stepping through it line by line in your editor. To run the whole thing and see
# each line next to what it produced:
#
#   Rscript -e 'source("code/walkthrough.R", echo = TRUE, max.deparse.length = 1e4)'
#
# Small tree and short chains, so it finishes in about a minute. The NUMBERS are not the
# campaign's. Every STEP is.
# =================================================================================================

suppressPackageStartupMessages(library(ape))
# Where the other scripts are. Works from the study directory or from inside code/.
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (is.na(here) || !file.exists(file.path(here, "rubin_lib.R")))
  here <- if (file.exists("code/rubin_lib.R")) "code" else "."

for (f in c("campaign_gnn_off_lib.R", "rubin_lib.R", "rubin_covariate.R",
            "rubin_completecase.R", "rubin_bace.R", "rubin_discrete.R")) source(file.path(here, f))

# ---- settings. Change any of these and run it again. --------------------------------------------
N       <- 80       # species
LAMBDA  <- 0.7      # phylogenetic signal: 0 = none, 1 = pure Brownian motion
RHO     <- 0.5      # correlation between the traits
RHO_COV <- 0.8      # correlation between the auxiliary covariate and trait c1
MISS    <- "clade"  # how gaps are made: "mcar" at random, "mar" driven by d1, "clade" whole groups
SEED    <- 1


# =================================================================================================
# STEP 1.  Simulate the data.
#          make_cell() grows a tree, evolves seven traits along it, builds one auxiliary
#          covariate, then hides 30% of the trait values -- and keeps the truth so we can score.
# =================================================================================================
cell <- make_cell("types_mixed", N, SEED, miss_frac = 0.30, miss = MISS, lambda = LAMBDA,
                  rho = RHO, thresholds = "fixed", driver = TRUE, n_cov = 1L, rho_cov = RHO_COV)

str(cell, max.level = 1, give.attr = FALSE)   # what you get back: 10 slots

head(cell$truth, 4)                                        # the truth. BACE never sees this.
head(cell$df_miss[!complete.cases(cell$df_miss), ], 4)     # what BACE sees: NA where we hid a value

# We asked for traits correlated RHO. One dataset will not show exactly that: both traits evolved
# on the same tree, so the effective sample size is far below N and the realised value scatters.
cor(cell$truth$c1, cell$truth$c2)


# =================================================================================================
# STEP 2.  The auxiliary covariate -- the "body size" you measured for every species.
#          It is a predictor, never a target: it lives outside `truth`, so it is never hidden
#          and never scored. covariate_check() reports whether it came out as strong as asked.
# =================================================================================================
covariate_check(cell$covs, cell$truth, target_names = "c1", rho_cov = RHO_COV)

anyNA(cell$covs)                   # FALSE -- a covariate with gaps would be pointless
"cov1" %in% names(cell$truth)      # FALSE -- not a trait, so never imputed and never scored


# =================================================================================================
# STEP 3.  Where the gaps went.
# =================================================================================================
cell$realised_frac        # we asked for 0.30
colSums(cell$mask)        # hidden cells per trait. d1 is the MAR driver and is never hidden.


# =================================================================================================
# STEP 4.  The competitor: drop every species with a gap, analyse what is left.
#          Under "mcar" expect about 8 species in 100 to survive (0.7^7 -- a species must dodge
#          all seven draws). Under "clade" most survive, but they are a biased sample.
# =================================================================================================
kept <- complete_species(cell$df_miss)
length(kept) / N                   # the survival rate IS the result here, not a diagnostic

cc <- score_complete_case(cell$df_miss, cell$tree, truth_rho = RHO)
cc$status                          # "scored", or "too_few" when it refuses rather than guessing
cc$estimands                       # truth is RHO; `covered` says whether the interval caught it


# =================================================================================================
# STEP 5.  BACE: impute the gaps instead of discarding the species.
# =================================================================================================
fb <- fit_bace_mi(cell, M = 5L, nitt = 3000, burnin = 500, thin = 10, runs = 3)

# One model per trait: that trait on the left, every other trait plus the covariate on the right.
# cov1 appears on every right-hand side and on no left-hand side. That is the check that the
# covariate reached the model as a predictor and was not quietly given a model of its own.
cat(unlist(fb$fixformula), sep = "\n")

sets <- mi_bace_shipped(fb$outb, cell$df_miss)   # the M completed datasets

# The M datasets MUST differ. That difference is the missing-data uncertainty, and it is the
# whole reason for imputing M times rather than once. A spread of zero means it is being thrown
# away -- which is exactly the bug this study found in 2026-09.
i <- which(cell$mask[, "c1"])[1]                 # one hidden c1 value
imputed <- sapply(sets, function(d) d$c1[i])
imputed                                          # the five imputations of it
cell$truth$c1[i]                                 # what it really was
sd(imputed)                                      # must not be 0


# =================================================================================================
# STEP 6.  Score the imputations against the truth we kept back.
#            zRMSE    how far off, in units of the trait's own spread. 0 = perfect,
#                     1 = no better than filling in the average.
#            coverage how often the 95% interval held the true value. Should be about 0.95.
# =================================================================================================
score <- data.frame()
for (v in c("c1", "c2", "cnt", "prp")) {
  idx   <- which(cell$mask[, v])
  draws <- t(sapply(sets, function(d) d[[v]][idx]))
  ci    <- rubin_cell_intervals(score_scale(v, draws)$x)
  tr    <- score_scale(v, cell$truth[[v]][idx])$x
  sdt   <- sd(score_scale(v, cell$truth[[v]][!cell$mask[, v]])$x)
  score <- rbind(score, data.frame(
    trait    = v,
    hidden   = length(idx),
    zRMSE    = sqrt(mean(((tr - ci$mean) / sdt)^2)),
    coverage = mean(tr >= ci$lower & tr <= ci$upper),
    width    = mean((ci$upper - ci$lower) / sdt)))
}
score

# The discrete traits are scored differently -- accuracy with credit for ties, and a set of
# categories instead of an interval. See code/rubin_discrete.R.
score_discrete("bace", sets, cell$truth, cell$mask)[, c("trait", "accuracy", "set_coverage")]
