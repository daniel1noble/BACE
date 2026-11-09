# code/rubin_completecase.R
#
# THE COMPLETE-CASE ARM -- the thing most people actually do with missing data.
#
# WHAT IT IS, in one sentence: throw away every species that has a gap in any trait, then run the
# analysis on whatever species are left.
#
# WHY WE RUN IT: it is the competitor. The manuscript's first claim is that BACE beats throwing
# data away when the gaps are not random. That claim has to be measured against something, and
# this is the something. Without this arm the claim is an assertion.
#
# WHAT IT IS NOT: the campaign already has an arm called "complete". That one is the ORACLE -- the
# true, un-holed data, which nobody has in real life. It is the best-possible benchmark. This file
# is the WORST-common-practice comparison. We need all three:
#
#     oracle          the truth before we punched holes in it    (how well could anyone do?)
#     complete-case   drop species with any gap                  (what people do now)
#     BACE            impute the gaps                            (what we are proposing)
#
# THE TWO WAYS COMPLETE-CASE GOES WRONG, which this arm must be able to tell apart:
#
#   1. ATTRITION -- you simply run out of species. With 7 traits each missing 30% at random and
#      independently, a species survives only if it dodged every one of the 7 draws:
#      0.7^7 = 0.082. So about 8 species out of 100 survive. The analysis is not biased, it is
#      just hopeless, and sometimes impossible.
#   2. BIAS -- the survivors are not a random sample of the species. This is the interesting
#      failure, and it is what the MAR and clade mechanisms are for.
#
# Those are different failures and the paper should not confuse them, so EVERY result from this
# arm carries `n_complete` (how many species survived) alongside the estimate. A reader can then
# see whether a wide interval came from having 8 species or from the survivors being odd.
#
# WHEN WE REFUSE TO SCORE: fewer than `min_species` survivors. A PGLS slope on 4 species is a
# number, but it is not an estimate, and quietly including it would move the arm's average around
# for no good reason. We return status = "too_few" instead, and the campaign reports what fraction
# of datasets ended that way. That refusal rate is itself a result -- arguably the headline one
# under MCAR.
#
# HOW TO READ THIS FILE: two functions, run in order. complete_species() decides who survives;
# score_complete_case() analyses the survivors. Each is annotated step by step.
#
# Assumes rubin_lib.R has been sourced (est_pgls_slope_fast, est_phylo_cor).


#' FUNCTION 1 of 2 -- find the species with no gaps at all.
#'
#' @param df_miss the dataset with holes in it (make_cell()$df_miss)
#' @param ignore  columns that are not traits and so should not count against a species. The MAR
#'   driver `d1` is fully observed by construction, so this is belt-and-braces, but it makes the
#'   rule explicit rather than implicit.
#' @return character vector of species names (rownames) that are complete across all trait columns
complete_species <- function(df_miss, ignore = "d1") {

  ## STEP 1. Decide which columns count. A species is not disqualified by a column that was never
  ##         a trait in the first place.
  trait_cols <- setdiff(names(df_miss), ignore)

  ## STEP 2. Mark the rows with no NA anywhere in those columns. This single line IS the
  ##         complete-case rule; everything else in this file is bookkeeping around it.
  has_no_gap <- stats::complete.cases(df_miss[, trait_cols, drop = FALSE])

  ## STEP 3. Return their names, not their row numbers, because the tree is indexed by name.
  rownames(df_miss)[has_no_gap]
}


#' FUNCTION 2 of 2 -- run the campaign's two estimands on the surviving species only.
#'
#' The estimands and their intervals are deliberately IDENTICAL to the other arms -- PGLS slope of
#' c2 on c1 under Pagel's lambda, and the GLS-whitened phylogenetic correlation of c1 and c2 -- so
#' the three arms differ in the data they saw and in nothing else.
#'
#' A note on the tree. We do NOT prune it. `est_pgls_slope_fast()` and `est_phylo_cor()` both build
#' their correlation matrix as `cov2cor(vcv(tree))[sp, sp]`, i.e. they take the sub-matrix for
#' whichever species you hand them. Dropping tips from an ultrametric tree does not change the
#' covariance among the tips that remain, so sub-setting the full tree's matrix gives exactly the
#' same answer as pruning the tree first, with one less thing to get wrong.
#'
#' A note on the interval. Complete-case produces ONE analysis, not M imputations, so there is no
#' between-imputation variance to pool: the interval is the ordinary t interval, estimate +/-
#' t(df) * se. That is the same formula the imputation arms collapse to when the between-imputation
#' variance is zero, so the comparison is fair rather than rigged either way.
#'
#' @param df_miss the dataset with holes (make_cell()$df_miss)
#' @param tree the full phylogeny
#' @param truth_rho the population value both estimands equal in this DGP, so coverage is scored
#'   against it (the same convention rubin_cell.R uses for every other arm)
#' @param min_species refuse to score below this many survivors (default 10)
#' @param ignore non-trait columns, passed to complete_species()
#' @return list(estimands = a two-row data.frame shaped like the other arms' estimand tables,
#'              n_complete = integer, frac_kept = numeric, status = "scored" | "too_few" | "failed")
score_complete_case <- function(df_miss, tree, truth_rho, min_species = 10L, ignore = "d1") {

  ## STEP 1. Who survives? Count them, and record the survival fraction -- this number is a result,
  ##         not a diagnostic, because under MCAR it is most of the story.
  kept      <- complete_species(df_miss, ignore = ignore)
  n_cc      <- length(kept)
  frac_kept <- n_cc / nrow(df_miss)

  ## STEP 2. Refuse rather than pretend. Too few survivors and we return a status, not a number.
  ##         The campaign tallies how often this happens; a high rate is the finding.
  if (n_cc < min_species) {
    return(list(estimands = NULL, n_complete = n_cc, frac_kept = frac_kept, status = "too_few"))
  }

  ## STEP 3. Build the complete-case dataset: the original rows, minus everyone with a gap.
  ##         No imputation, no filling, nothing invented -- that is the whole point of this arm.
  cc <- df_miss[kept, , drop = FALSE]

  ## STEP 4. Fit the same two estimands the other arms report.
  ##         4a. PGLS slope of c2 on c1, with Pagel's lambda re-estimated on the survivors.
  ##         4b. Phylogenetic correlation of c1 and c2, reusing 4a's lambda -- exactly the
  ##             convention score_estimands() uses in rubin_cell.R, so the arms stay comparable.
  slope <- est_pgls_slope_fast(cc, tree)
  corr  <- est_phylo_cor(cc, tree, lambda = slope$lambda_hat)

  ## STEP 5. A fit can fail outright (too few survivors to pin lambda down, a degenerate column).
  ##         Record that as its own status; never let a non-finite number into the results table.
  if (!is.finite(slope$estimate) || !is.finite(slope$variance) || !is.finite(corr$z)) {
    return(list(estimands = NULL, n_complete = n_cc, frac_kept = frac_kept, status = "failed"))
  }

  ## STEP 6. Interval for the slope: ordinary t interval on (survivors - 2) degrees of freedom.
  se_slope <- sqrt(slope$variance)
  t_crit   <- stats::qt(0.975, df = n_cc - 2)
  lo_s <- slope$estimate - t_crit * se_slope
  hi_s <- slope$estimate + t_crit * se_slope

  ## STEP 7. Interval for the correlation: a normal interval on Fisher's z scale, where a
  ##         correlation is roughly normal, then transformed back to the r scale with tanh().
  se_z <- sqrt(1 / (n_cc - 3))
  lo_r <- tanh(corr$z - stats::qnorm(0.975) * se_z)
  hi_r <- tanh(corr$z + stats::qnorm(0.975) * se_z)

  ## STEP 8. Assemble a table shaped exactly like every other arm's, so the aggregation scripts
  ##         can rbind() it without knowing this arm exists. `covered` is the coverage indicator:
  ##         did the interval contain the population value?
  est <- rbind(
    data.frame(arm = "complete_case", estimand = "slope", estimate = slope$estimate, se = se_slope,
               lower = lo_s, upper = hi_s, df = n_cc - 2, fmi = NA_real_, m_ok = 1L,
               truth = truth_rho, complete_data = NA_real_,
               covered = lo_s <= truth_rho & truth_rho <= hi_s),
    data.frame(arm = "complete_case", estimand = "cor", estimate = corr$r, se = se_z,
               lower = lo_r, upper = hi_r, df = n_cc - 3, fmi = NA_real_, m_ok = 1L,
               truth = truth_rho, complete_data = NA_real_,
               covered = lo_r <= truth_rho & truth_rho <= hi_r)
  )

  list(estimands = est, n_complete = n_cc, frac_kept = frac_kept, status = "scored")
}
