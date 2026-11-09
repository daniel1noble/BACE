# code/oracle_floor.R
#
# Theoretical best-achievable (oracle) zRMSE per simulation scenario, for the
# zRMSE panel's per-scenario reference lines (Dan, 2026-11-05).
#
# The DGP (campaign_gnn_off_lib.R::sim_latents) draws the n x 8 latent matrix L
# with Vec(L) ~ N(0, Sigma_rho (x) C_lambda), so the minimum possible error for
# any missing cell is a multivariate-normal conditional variance:
#
#   precision  P = Sigma^{-1} (x) C^{-1}
#   v_missing  = diag( solve(P[miss, miss]) )   given the OBSERVED cells only
#                (30% MCAR per trait, covariate complete - same design as the
#                campaign), averaged over random trees and masks.
#
# This is a LOWER bound: the oracle sees the latent values of the discrete and
# noisy traits (cnt, prp, bin, ord, cat3), where a real method sees only their
# observed transformations. Count and proportion floors propagate the latent
# floor through the trait's own observation noise (Poisson; logit + N(0, 0.3^2)),
# scored on the same scale as the campaign (raw counts / logits), divided by the
# population SD of observed values.
#
#   Rscript code/oracle_floor.R [OUT_CSV]
#
# Output: one row per n x lambda x rho x trait (c1, c2, cnt, prp):
#   floor      - oracle zRMSE under the campaign's 30% MCAR design
#   floor_full - oracle zRMSE if every other value were observed (closed form
#                1 / (Sigma^{-1}_jj C^{-1}_ii); tighter, reference only)

suppressMessages(library(ape))
source(file.path("code", "campaign_gnn_off_lib.R"))   # phylo_corr()

a <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(a) >= 1) a[1] else "results_2026-10/oracle_floor.csv"

# The grid these reference lines are computed on MUST match the campaign's grid. If it does not,
# the floor/ceiling lines are simply absent at the missing n or lambda and the figure panel loses
# them silently -- ggplot drops unmatched rows with at most a warning.
# 2026-11 design (.agents/plans/simulation-plan-2026-11.md): n = 100, 500; lambda = 0.3, 0.7, 0.9;
# rho = 0 and 0.5 (the three information regimes use only these two trait-trait values).
NS <- c(100L, 500L); LAMBDAS <- c(0.3, 0.7, 0.9); RHOS <- c(0, 0.5)
TREES <- 4L; MASKS <- 2L; MISS <- 0.3; K <- 8L; DRIVER <- 8L
A_CNT <- 1.5; B_CNT <- 0.8; SD_PRP_ERR <- 0.3
set.seed(20261105)

sigma_rho <- function(rho) {
  S <- matrix(rho, K, K); diag(S) <- 1
  S[DRIVER, -DRIVER] <- 0.35; S[-DRIVER, DRIVER] <- 0.35
  S
}

# E[Var(count | rest)] with latent conditional variance v, marginal latent var 1
cnt_floor_var <- function(v, a = A_CNT, b = B_CNT)
  exp(a + b^2 / 2) + exp(2 * a + b^2 * (2 - v)) * (expm1(b^2 * v))
cnt_marg_var <- exp(A_CNT + B_CNT^2 / 2) + exp(2 * A_CNT + B_CNT^2) * expm1(B_CNT^2)

# thresholds = "fixed" in the campaign: population quantiles of the N(0,1) liability
CUTS <- list(bin = 0, ord = qnorm(c(.25, .5, .75)), cat3 = qnorm(c(1/3, 2/3)))
CHANCE <- c(bin = 0.5, ord = 0.25, cat3 = 1/3)
GH_U <- (seq_len(200) - 0.5) / 200   # equal-probability grid for E over the conditional mean

# Bayes ceiling for modal-category accuracy: liability | rest ~ N(m, v), m ~ N(0, 1 - v);
# the oracle picks the most probable category, so P(correct | m) = max_k p_k(m, v)
acc_ceiling <- function(v, cuts) {
  b <- c(-Inf, cuts, Inf)
  vapply(seq_along(v), function(i) {
    s <- sqrt(v[i]); m <- qnorm(GH_U) * sqrt(max(1 - v[i], 0))
    pk <- vapply(seq_len(length(b) - 1L), function(k)
      pnorm((b[k + 1L] - m) / s) - pnorm((b[k] - m) / s), numeric(length(m)))
    mean(apply(pk, 1L, max))
  }, numeric(1))
}

rows <- list()
for (n in NS) {
  n_miss <- round(MISS * n)
  for (tr in seq_len(TREES)) {
    tree <- ape::rcoal(n)
    tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
    for (lambda in LAMBDAS) {
      Cl <- phylo_corr(tree, lambda, "BM")
      Cinv <- chol2inv(chol(Cl))
      for (rho in RHOS) {
        Sinv <- solve(sigma_rho(rho))
        for (mk in seq_len(MASKS)) {
          # 30% MCAR per trait (traits 1..7), covariate complete - as in the campaign
          miss_i <- unlist(lapply(1:7, function(j) sample.int(n, n_miss)))
          miss_j <- rep(1:7, each = n_miss)
          # conditional covariance of missing given observed = inverse of the
          # missing-block of the precision P = Sinv (x) Cinv (never formed fully)
          Pmm <- Sinv[miss_j, miss_j] * Cinv[miss_i, miss_i]
          v <- diag(chol2inv(chol(Pmm)))
          v_full <- 1 / (diag(Sinv)[miss_j] * diag(Cinv)[miss_i])
          tn <- c("c1", "c2", "cnt", "prp", "bin", "ord", "cat3")
          for (j in 1:7) {
            vj <- v[miss_j == j]; vfj <- v_full[miss_j == j]
            rows[[length(rows) + 1]] <- if (j <= 4) data.frame(
              n = n, lambda = lambda, rho = rho, tree = tr, mask = mk, trait = tn[j],
              mse = switch(j, mean(vj), mean(vj),
                           mean(cnt_floor_var(vj)) / cnt_marg_var,
                           (mean(vj) + SD_PRP_ERR^2) / (1 + SD_PRP_ERR^2)),
              mse_full = switch(j, mean(vfj), mean(vfj),
                                mean(cnt_floor_var(vfj)) / cnt_marg_var,
                                (mean(vfj) + SD_PRP_ERR^2) / (1 + SD_PRP_ERR^2)),
              acc = NA_real_, acc_full = NA_real_)
            else data.frame(
              n = n, lambda = lambda, rho = rho, tree = tr, mask = mk, trait = tn[j],
              mse = NA_real_, mse_full = NA_real_,
              acc = mean(acc_ceiling(vj, CUTS[[tn[j]]])),
              acc_full = mean(acc_ceiling(vfj, CUTS[[tn[j]]])))
          }
        }
      }
    }
  }
  cat("n =", n, "done\n")
}

d <- do.call(rbind, rows)
ac <- aggregate(cbind(mse, mse_full) ~ n + lambda + rho + trait, d, mean, na.action = na.omit)
ac$floor <- sqrt(ac$mse); ac$floor_full <- sqrt(ac$mse_full)
ac$ceiling <- NA_real_; ac$ceiling_full <- NA_real_; ac$chance <- NA_real_
ad <- aggregate(cbind(acc, acc_full) ~ n + lambda + rho + trait, d, mean, na.action = na.omit)
ad$floor <- NA_real_; ad$floor_full <- NA_real_
ad$ceiling <- ad$acc; ad$ceiling_full <- ad$acc_full; ad$chance <- CHANCE[ad$trait]
cols <- c("trait", "n", "lambda", "rho", "floor", "floor_full", "ceiling", "ceiling_full", "chance")
agg <- rbind(ac[cols], ad[cols])
agg <- agg[order(match(agg$trait, c("c1","c2","cnt","prp","bin","ord","cat3")), agg$n, agg$lambda, agg$rho), ]
write.csv(agg, out_csv, row.names = FALSE)
cat("written:", out_csv, "\n")
print(agg, digits = 2, row.names = FALSE)
