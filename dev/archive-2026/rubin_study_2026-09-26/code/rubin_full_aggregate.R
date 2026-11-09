# code/rubin_full_aggregate.R
#
# Aggregate the COMPLETED 2026-10 BACE-only campaign (750 rds in results_2026-10/)
# into flat tables for the manuscript figures.
#   Rscript code/rubin_full_aggregate.R [RESULTS_DIR] [OUT_DIR]
#
# Unit of replication: one simulated dataset (n, lambda, rho, seed). Failures are
# recorded INSIDE the rds ($errors, arm absent from $cells/$estimands); no seeds were
# replaced, so per-cell n_ok is the honest denominator and n_fail the survivorship
# ledger (all failures are lambda=1 "Mixed model equations singular" in the bace arm).
#
# Outputs (CSV, written to OUT_DIR):
#   full_cells.csv      per-value metrics, dataset level (c1/c2/cnt on raw, prp logit)
#   full_disc_cells.csv discrete per-value metrics, dataset level (bin/ord/cat3)
#   full_estimands.csv  downstream slope/cor, dataset level (bace + complete arms)
#   agg_cells.csv       per-value summaries by trait x n x lambda x rho (+ MC SE)
#   agg_disc.csv        discrete summaries by trait x n x lambda x rho (+ MC SE)
#   agg_down.csv        downstream summaries by estimand x arm x n x lambda x rho
#   agg_failures.csv    per cell: n_total, n_ok, n_fail, error class
#   sanity_full.txt     file counts, failure ledger, headline numbers

a <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(a) >= 1) a[1] else "results_2026-10"
out_dir <- if (length(a) >= 2) a[2] else res_dir

fs <- list.files(res_dir, pattern = "^rubin_types_mixed_.*[.]rds$", full.names = TRUE)
cat(sprintf("reading %d rds files from %s\n", length(fs), res_dir))

cells <- list(); disc <- list(); est <- list(); fail <- list()
for (f in fs) {
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(x)) { cat("UNREADABLE:", f, "\n"); next }
  id <- data.frame(n = x$n, lambda = x$lambda, rho = x$rho, seed = x$seed)
  err_names <- names(x$errors)
  err_msg <- if (length(x$errors)) paste(vapply(x$errors, function(e) {
    m <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
    substr(m[1], 1, 60) }, ""), collapse = "; ") else ""
  fail[[f]] <- cbind(id, failed = length(x$errors) > 0,
                     err_arm = paste(err_names, collapse = ";"), err_msg = err_msg)
  if (!is.null(x$cells))      cells[[f]] <- cbind(id, x$cells)
  if (!is.null(x$disc_cells)) disc[[f]]  <- cbind(id, x$disc_cells)
  if (!is.null(x$estimands))  est[[f]]   <- cbind(id, x$estimands)
}
cells <- do.call(rbind, cells); disc <- do.call(rbind, disc)
est <- do.call(rbind, est); fail <- do.call(rbind, fail)
rownames(cells) <- rownames(disc) <- rownames(est) <- rownames(fail) <- NULL

se <- function(v) sd(v) / sqrt(length(v))

# per-value: pool c1+c2 is done at figure time; keep traits separate here
agg_by <- function(d, by, vars) do.call(rbind, lapply(split(d, d[by], drop = TRUE), function(g) {
  r <- g[1, by, drop = FALSE]
  for (v in vars) {
    r[[v]] <- mean(g[[v]], na.rm = TRUE)
    r[[paste0(v, "_se")]] <- se(g[[v]][!is.na(g[[v]])])
  }
  r$datasets <- nrow(g); r
}))

agg_cells <- agg_by(cells[cells$arm == "bace", ], c("trait", "n", "lambda", "rho"),
                    c("coverage", "zRMSE", "width", "interval_score"))
agg_disc <- if (!is.null(disc) && nrow(disc)) {
  agg_by(disc[disc$arm == "bace", ], c("trait", "n", "lambda", "rho"),
         c("accuracy", "brier", "ece", "set_coverage", "set_size"))
} else NULL

est <- transform(est, bias = estimate - truth, ciw = upper - lower,
                 covered = as.numeric(covered))
agg_down <- do.call(rbind, lapply(
  split(est, est[c("arm", "estimand", "n", "lambda", "rho")], drop = TRUE),
  function(g) data.frame(g[1, c("arm", "estimand", "n", "lambda", "rho")],
    coverage = mean(g$covered), coverage_se = se(g$covered),
    bias = mean(g$bias), bias_se = se(g$bias), rmse = sqrt(mean(g$bias^2)),
    ci_width = mean(g$ciw),
    fmi_median = if (all(is.na(g$fmi))) NA else median(g$fmi, na.rm = TRUE),
    datasets = nrow(g), row.names = NULL)))

agg_failures <- do.call(rbind, lapply(split(fail, fail[c("n", "lambda", "rho")], drop = TRUE),
  function(g) data.frame(n = g$n[1], lambda = g$lambda[1], rho = g$rho[1],
    n_total = nrow(g), n_ok = sum(!g$failed), n_fail = sum(g$failed),
    err_arms = paste(unique(g$err_arm[g$failed]), collapse = ";"), row.names = NULL)))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
wr <- function(d, nm) if (!is.null(d)) write.csv(d, file.path(out_dir, nm), row.names = FALSE)
wr(cells, "full_cells.csv"); wr(disc, "full_disc_cells.csv"); wr(est, "full_estimands.csv")
wr(agg_cells, "agg_cells.csv"); wr(agg_disc, "agg_disc.csv"); wr(agg_down, "agg_down.csv")
wr(agg_failures, "agg_failures.csv")

sanity <- c(
  sprintf("files read: %d; datasets with recorded failures: %d (expect 20, all lambda=1 bace_fit singular-MME)",
          length(fs), sum(fail$failed)),
  sprintf("cells rows: %d; disc rows: %d; estimand rows: %d",
          nrow(cells), if (is.null(disc)) 0L else nrow(disc), nrow(est)),
  capture.output(print(agg_failures[agg_failures$n_fail > 0, ], row.names = FALSE)),
  "--- bace per-value coverage by n (c1+c2 cell-weighted) ---",
  capture.output({
    cc <- cells[cells$arm == "bace" & cells$trait %in% c("c1", "c2"), ]
    pd <- aggregate(cbind(cv = coverage * n_cells, w = n_cells) ~ n + seed + lambda + rho,
                    data = cc, FUN = sum)
    pd$coverage <- pd$cv / pd$w
    print(aggregate(coverage ~ n, data = pd, FUN = function(v)
      round(c(mean = mean(v), se = se(v)), 4)), row.names = FALSE)
  }),
  "--- downstream (bace arm) coverage/bias by estimand x n ---",
  capture.output(print(do.call(rbind, lapply(
    split(est[est$arm == "bace", ], est[est$arm == "bace", c("estimand", "n")], drop = TRUE),
    function(g) data.frame(estimand = g$estimand[1], n = g$n[1],
      coverage = round(mean(g$covered), 3), bias = round(mean(g$bias), 4),
      datasets = nrow(g), row.names = NULL))), row.names = FALSE)))
writeLines(sanity, file.path(out_dir, "sanity_full.txt"))
cat(sanity, sep = "\n")
