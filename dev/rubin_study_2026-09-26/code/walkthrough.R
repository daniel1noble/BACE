# =================================================================================================
# walkthrough.R -- one dataset, every step shown, so you can check the simulation yourself.
#
# Run it:   Rscript code/walkthrough.R
#           Rscript code/walkthrough.R --n 150 --lambda 0.3 --miss clade --rho_cov 0.8
#
# It takes about a minute. It uses a small tree and short MCMC chains, so the NUMBERS are not the
# campaign's numbers -- but every STEP is the same step the campaign takes, in the same order,
# calling the same functions. If something here looks wrong, it is wrong in the campaign too.
#
# What it shows, in order:
#   1. the tree and the traits          -- what was simulated, and the correlations it really has
#   2. the auxiliary covariate          -- is it there, and is it as strong as we asked
#   3. the gaps                         -- how many, and where, under the chosen mechanism
#   4. complete-case                    -- who survives being dropped, and what that does
#   5. BACE                             -- what models it was given, and what it imputed
#   6. scoring                          -- the imputed values against the truth we kept
# =================================================================================================

suppressPackageStartupMessages({ library(ape) })

a <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default) { i <- match(paste0("--", flag), a); if (is.na(i)) default else a[i + 1L] }
N       <- as.integer(opt("n", 80))
LAMBDA  <- as.numeric(opt("lambda", 0.7))
RHO     <- as.numeric(opt("rho", 0.5))
RHO_COV <- as.numeric(opt("rho_cov", 0.8))
MISS    <- opt("miss", "clade")
SEED    <- as.integer(opt("seed", 1))

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
for (f in c("campaign_gnn_off_lib.R", "rubin_lib.R", "rubin_covariate.R",
            "rubin_completecase.R", "rubin_bace.R", "rubin_discrete.R")) source(file.path(here, f))

rule <- function(title) cat("\n", strrep("=", 86), "\n", title, "\n", strrep("=", 86), "\n", sep = "")
say  <- function(...) cat(sprintf(...), "\n", sep = "")

cat("BACE version in use:", as.character(utils::packageVersion("BACE")),
    " | settings: n =", N, " lambda =", LAMBDA, " rho =", RHO, " rho_cov =", RHO_COV,
    " missingness =", MISS, "\n")


# -------------------------------------------------------------------------------------------------
rule("STEP 1 of 6 -- simulate a tree and seven traits, and keep the truth")
# -------------------------------------------------------------------------------------------------
# make_cell() does three things at once: grows a tree, evolves correlated traits along it, and
# then hides some values. We ask it for one auxiliary covariate as well.
cell <- make_cell("types_mixed", N, SEED, miss_frac = 0.30, miss = MISS, lambda = LAMBDA,
                  rho = RHO, thresholds = "fixed", driver = TRUE,
                  n_cov = 1L, rho_cov = RHO_COV)

say("tree       : %d species, ultrametric = %s", ape::Ntip(cell$tree), ape::is.ultrametric(cell$tree))
say("traits     : %s", paste(names(cell$truth), collapse = ", "))
cat("\nwhat each trait is, and the first few true values:\n")
for (v in names(cell$truth)) {
  x <- cell$truth[[v]]
  say("  %-5s %-22s %s", v,
      if (is.ordered(x)) sprintf("ordered factor (%d)", nlevels(x))
      else if (is.factor(x)) sprintf("factor (%d levels)", nlevels(x))
      else if (is.integer(x)) "count (integer)" else "continuous",
      paste(utils::head(if (is.numeric(x)) format(round(x, 2), trim = TRUE) else as.character(x), 5),
            collapse = ", "))
}

# CHECK IT YOURSELF: we asked for traits correlated RHO. Does c1 really correlate RHO with c2?
cat("\ncheck: we asked for trait-trait correlation", RHO, "\n")
say("  realised corr(c1, c2) in this one dataset = %.3f", cor(cell$truth$c1, cell$truth$c2))
say("  (one dataset scatters widely because both traits evolved on the same tree --")
say("   the effective sample size is far below %d. The AVERAGE over datasets is %.2f.)", N, RHO)


# -------------------------------------------------------------------------------------------------
rule("STEP 2 of 6 -- the auxiliary covariate: the 'body size' you measured for everybody")
# -------------------------------------------------------------------------------------------------
chk <- covariate_check(cell$covs, cell$truth, target_names = "c1", rho_cov = RHO_COV)
if (is.null(chk)) {
  say("no covariate in this run (n_cov = 0).")
} else {
  print(chk, row.names = FALSE)
  say("")
  say("it is COMPLETE          : %s  (no gaps -- that is the point of an auxiliary)",
      !any(is.na(cell$covs$cov1)))
  say("it is NOT a trait       : %s  (so it is never hidden and never scored)",
      !("cov1" %in% names(cell$truth)))
  say("it IS phylogenetic      : built from noise that evolved on the same tree, not white noise")
}


# -------------------------------------------------------------------------------------------------
rule(sprintf("STEP 3 of 6 -- punch the holes, mechanism = %s", MISS))
# -------------------------------------------------------------------------------------------------
say("target missing fraction : 30%%")
say("realised                : %.1f%%", 100 * cell$realised_frac)
cat("\ngaps per trait:\n")
for (v in colnames(cell$mask)) say("  %-5s %3d of %d hidden (%.0f%%)", v, sum(cell$mask[, v]), N,
                                   100 * mean(cell$mask[, v]))
say("")
say("d1 is the variable that CAUSES the gaps under 'mar'. It is never hidden itself: %d gaps.",
    sum(is.na(cell$df_miss$d1)))
if (MISS == "clade") {
  say("under 'clade', whole groups of related species lose everything at once -- which is why")
  say("the survivors are a biased sample rather than a smaller random one.")
}


# -------------------------------------------------------------------------------------------------
rule("STEP 4 of 6 -- the competitor: drop every species that has a gap")
# -------------------------------------------------------------------------------------------------
kept <- complete_species(cell$df_miss)
say("species with no gap at all : %d of %d  (%.0f%%)", length(kept), N, 100 * length(kept) / N)
if (MISS == "mcar") {
  say("expected under MCAR        : %.0f%%  (0.7^7 = each species must dodge all seven draws)",
      100 * 0.7^7)
}
cc <- score_complete_case(cell$df_miss, cell$tree, truth_rho = RHO)
say("status                     : %s", cc$status)
if (!is.null(cc$estimands)) {
  cat("\nwhat complete-case concludes (the truth is", RHO, "):\n")
  print(cc$estimands[, c("estimand", "estimate", "lower", "upper", "df", "covered")],
        row.names = FALSE, digits = 3)
}


# -------------------------------------------------------------------------------------------------
rule("STEP 5 of 6 -- BACE: impute the gaps instead of discarding them")
# -------------------------------------------------------------------------------------------------
say("(short chains here so this finishes in a minute -- the campaign uses nitt 50000)")
fb <- fit_bace_mi(cell, M = 5L, nitt = 3000, burnin = 500, thin = 10, runs = 3)

cat("\nthe imputation models BACE was given -- one per trait, each on all the others:\n")
for (f in fb$fixformula) say("  %s", f)
say("")
say("the covariate appears on the RIGHT of every model : %s",
    all(grepl("cov1", unlist(fb$fixformula))))
say("the covariate is never a RESPONSE                 : %s",
    !any(grepl("^cov", vapply(fb$fixformula, function(f) sub(" ~.*", "", f), character(1)))))

sets <- mi_bace_shipped(fb$outb, cell$df_miss)
say("")
say("BACE produced %d completed datasets. They must DIFFER from each other -- that difference is", length(sets))
say("the missing-data uncertainty, and it is what makes the intervals honest.")
i <- which(cell$mask[, "c1"])[1]
vals <- vapply(sets, function(d) d[["c1"]][i], numeric(1))
say("  one hidden c1 value, imputed %d times : %s", length(vals),
    paste(sprintf("%.2f", vals), collapse = ", "))
say("  the truth we hid                      : %.2f", cell$truth$c1[i])
say("  spread across imputations (sd)        : %.3f   <- zero here would mean the uncertainty", sd(vals))
say("                                                     is being thrown away")


# -------------------------------------------------------------------------------------------------
rule("STEP 6 of 6 -- score the imputations against the truth we kept")
# -------------------------------------------------------------------------------------------------
say("zRMSE : how far the imputed values are from the truth, in units of the trait's own spread.")
say("        0 = perfect. 1 = no better than filling in the average.")
say("coverage: how often the 95%% interval contained the true value. Should be about 0.95.")
cat("\n")
for (v in c("c1", "c2", "cnt", "prp")) {
  idx <- which(cell$mask[, v]); if (!length(idx)) next
  draws <- t(vapply(sets, function(s) as.numeric(s[[v]][idx]), numeric(length(idx))))
  sc <- score_scale(v, draws); ci <- rubin_cell_intervals(sc$x)
  tr <- score_scale(v, cell$truth[[v]][idx])$x
  sdt <- stats::sd(score_scale(v, cell$truth[[v]][!cell$mask[, v]])$x)
  say("  %-4s %3d hidden | zRMSE %.2f | coverage %.2f | mean interval width %.2f sd",
      v, length(idx), sqrt(mean(((tr - ci$mean) / sdt)^2)),
      mean(tr >= ci$lower & tr <= ci$upper), mean((ci$upper - ci$lower) / sdt))
}

rule("DONE")
say("Every step above is the step the campaign takes. To run the real thing:")
say("  DRY_RUN=1 bash code/rubin_campaign_local.sh    # see the task list first")
say("  CORES=7  bash code/rubin_campaign_local.sh")
