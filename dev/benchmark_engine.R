# =============================================================================
# benchmark_engine.R - shared cross-dataset BACE benchmarking infrastructure
# =============================================================================
#
# Each dataset (pantheria, amphibio, bien, globtherm, leptraits) has a thin
# wrapper script (03..07_benchmark_*.R) that loads the bundled .rda pair and
# calls `benchmark_dataset()`. The engine handles type-aware masking,
# fits via bace(), and computes type-aware metrics:
#
#   continuous / count : NRMSE, MAE_fit, MAE_raw, Pearson r, coverage95
#   categorical        : accuracy, balanced_accuracy, multiclass Brier
#   binary             : same as categorical with K=2
#   ordinal            : categorical metrics + mean abs level distance
#
# Outputs land in dev/benchmark_results/<dataset>/run_<scope>_<date>/:
#   summary_metrics.csv  - long-format per-trait metrics
#   phylo_signal.csv     - pre-imputation Pagel lambda / Blomberg K / D
#   run_info.csv         - dataset size, MCMC config, runtime, convergence
#
# References for metric definitions:
#   NRMSE / MAE / coverage : van Buuren 2018 FIMD ch.5; Stekhoven & Buhlmann
#                            2012 Bioinformatics 28:112
#   Brier score            : Brier 1950 MWR 78:1; Gneiting & Raftery 2007
#                            JASA 102:359
#   Balanced accuracy      : Brodersen et al. 2010 ICPR; robust to imbalance
#   Pagel lambda / K       : Pagel 1999 Nature 401:877; Blomberg et al. 2003
#                            Evolution 57:717
#   Fritz-Purvis D         : Fritz & Purvis 2010 Conserv. Biol. 24:1042
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
})

# ----------------------------------------------------------------------------
# Tree prep — make any phylo MCMCglmm-acceptable:
#   (1) Make node labels unique — MCMCglmm::inverseA requires all
#       tip + node labels to be distinct. Source trees and taxonomy-
#       built trees (globtherm, leptraits) carry collapsed-singleton
#       internals labelled "" which duplicate to many copies of "".
#       make.unique-ing preserves informative names where present
#       while killing duplicates.
#   (2) Root the tree if unrooted — pantheria_tree historically
#       shipped unrooted, and MCMCglmm::inverseA refuses unrooted
#       phylogenies.
#   (3) Replace zero-length edges with a small positive epsilon —
#       these can appear after ape::keep.tip subsampling.
#   (4) Force ultrametric (unconditionally, AFTER step 3 — epsilon on
#       tip edges breaks ultrametricity if applied later).
# ----------------------------------------------------------------------------
.prep_tree <- function(tree) {
  if (!is.null(tree$node.label)) {
    nl <- ifelse(nzchar(tree$node.label) & !is.na(tree$node.label),
                 tree$node.label, "N")
    tree$node.label <- make.unique(nl)
  }
  if (!ape::is.rooted(tree)) {
    tree <- ape::root.phylo(tree, outgroup = 1L, resolve.root = TRUE)
  }
  if (any(tree$edge.length <= 0)) {
    pos <- tree$edge.length[tree$edge.length > 0]
    eps <- if (length(pos)) min(pos) / 1e3 else 1e-6
    tree$edge.length[tree$edge.length <= 0] <- eps
  }
  if (!ape::is.ultrametric(tree)) {
    if (requireNamespace("phytools", quietly = TRUE)) {
      tree <- phytools::force.ultrametric(tree, method = "extend")
    } else {
      depths    <- ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))]
      extend    <- max(depths) - depths
      tip_edges <- match(seq_len(ape::Ntip(tree)), tree$edge[, 2])
      tree$edge.length[tip_edges] <- tree$edge.length[tip_edges] + extend
    }
  }
  tree
}

# ----------------------------------------------------------------------------
# Trait type detection: mirrors R/prep_functions.R .detect_type() so the
# benchmark uses the same types BACE will internally.
# ----------------------------------------------------------------------------
.detect_trait_type <- function(x) {
  if (is.ordered(x)) return("ordinal")
  if (is.factor(x)) {
    k <- nlevels(x)
    if (k <= 1L) return("constant")
    if (k == 2L) return("binary")
    return("categorical")
  }
  if (is.integer(x)) return("count")
  if (is.numeric(x)) return("continuous")
  "unknown"
}

# ----------------------------------------------------------------------------
# Phylogenetic signal (pre-imputation reference values per trait).
# ----------------------------------------------------------------------------
.phylo_signal_cont <- function(tree, values) {
  if (!requireNamespace("phytools", quietly = TRUE))
    return(list(lambda = NA_real_, K = NA_real_))
  x <- values[!is.na(values)]
  if (length(x) < 10L || stats::var(x) == 0)
    return(list(lambda = NA_real_, K = NA_real_))
  tr <- ape::keep.tip(tree, names(x))
  lambda <- tryCatch(
    suppressWarnings(phytools::phylosig(tr, x, method = "lambda")$lambda),
    error = function(e) NA_real_)
  K <- tryCatch(
    suppressWarnings(phytools::phylosig(tr, x, method = "K")),
    error = function(e) NA_real_)
  K_val <- if (is.list(K)) K$K else K
  list(lambda = lambda, K = as.numeric(K_val))
}

.phylo_signal_cat <- function(tree, values, species_names) {
  if (!requireNamespace("caper", quietly = TRUE)) return(NA_real_)
  chr  <- as.character(values)
  keep <- !is.na(chr) & nzchar(chr)
  chr  <- chr[keep]
  sp   <- species_names[keep]
  if (length(chr) < 10L) return(NA_real_)
  classes <- unique(chr)
  if (length(classes) < 2L) return(NA_real_)
  Ds <- vapply(classes, function(cl) {
    bin <- as.integer(chr == cl)
    if (sum(bin) < 2 || sum(bin) > length(bin) - 2) return(NA_real_)
    df <- data.frame(Species = sp, y = bin)
    tr_cl <- ape::keep.tip(tree, sp)
    cd <- tryCatch(
      caper::comparative.data(tr_cl, df, names.col = "Species",
                              na.omit = FALSE, warn.dropped = FALSE),
      error = function(e) NULL)
    if (is.null(cd)) return(NA_real_)
    tryCatch(
      suppressMessages(caper::phylo.d(cd, binvar = y, permut = 100)$DEstimate),
      error = function(e) NA_real_)
  }, numeric(1))
  mean(Ds, na.rm = TRUE)
}

# ----------------------------------------------------------------------------
# Type-aware random masking. Builds a long-format `truth` data.frame
# with one row per hidden cell.
# ----------------------------------------------------------------------------
.apply_mask <- function(masked_df, types, cont_miss_rate, cat_miss_rate) {
  truth <- list()
  for (v in names(types)) {
    type <- types[[v]]
    rate <- if (type %in% c("continuous", "count")) cont_miss_rate
            else                                     cat_miss_rate
    vals <- masked_df[[v]]
    obs  <- which(!is.na(vals))
    if (length(obs) < 10L) next
    n_mask <- floor(length(obs) * rate)
    if (n_mask == 0L) next
    mask_idx <- sample(obs, n_mask)
    is_cat <- type %in% c("ordinal", "categorical", "binary")
    truth[[v]] <- data.frame(
      species_tip = rownames(masked_df)[mask_idx],
      trait       = v,
      type        = type,
      true_value  = if (is_cat) as.character(vals[mask_idx])
                    else        as.numeric(  vals[mask_idx]),
      stringsAsFactors = FALSE
    )
    masked_df[[v]][mask_idx] <- NA
  }
  list(masked = masked_df, truth = do.call(rbind, truth))
}

# ----------------------------------------------------------------------------
# Continuous / count metrics. true_vals is on the comparison scale (log1p
# of raw if is_log else raw). sd_full is marginal sd on the fit scale.
# ----------------------------------------------------------------------------
.metrics_continuous <- function(imp_mat, true_vals, sd_full) {
  per_imp_nrmse <- apply(imp_mat, 2,
    function(iv) sqrt(mean((iv - true_vals)^2)) / sd_full)
  per_imp_mae   <- apply(imp_mat, 2,
    function(iv) mean(abs(iv - true_vals)))
  per_imp_cor   <- apply(imp_mat, 2, function(iv) {
    if (length(unique(true_vals)) > 1 && length(unique(iv)) > 1)
      suppressWarnings(stats::cor(iv, true_vals)) else NA_real_
  })
  lo <- apply(imp_mat, 1, stats::quantile, probs = 0.025,
              na.rm = TRUE, names = FALSE)
  hi <- apply(imp_mat, 1, stats::quantile, probs = 0.975,
              na.rm = TRUE, names = FALSE)
  list(
    nrmse       = mean(per_imp_nrmse, na.rm = TRUE),
    mae_fit     = mean(per_imp_mae,   na.rm = TRUE),
    correlation = mean(per_imp_cor,   na.rm = TRUE),
    coverage95  = mean(true_vals >= lo & true_vals <= hi)
  )
}

# Majority-vote class assignment per hidden cell, breaking ties at random.
.majority_vote <- function(im) {
  apply(im, 1, function(row) {
    row <- row[!is.na(row)]
    if (length(row) == 0L) return(NA_character_)
    tab <- table(row)
    top <- names(tab)[tab == max(tab)]
    if (length(top) == 1L) top else sample(top, 1L)
  })
}

# Categorical / binary metrics: accuracy, balanced_accuracy, multiclass Brier.
.metrics_categorical <- function(im, true_chr, levels_all) {
  vote <- .majority_vote(im)
  accuracy <- mean(vote == true_chr, na.rm = TRUE)
  classes <- unique(true_chr)
  recalls <- vapply(classes, function(cl) {
    is_cl <- true_chr == cl
    if (!any(is_cl)) NA_real_ else mean(vote[is_cl] == cl, na.rm = TRUE)
  }, numeric(1))
  bal_acc <- mean(recalls, na.rm = TRUE)
  prob_mat <- sapply(levels_all, function(cl) rowMeans(im == cl, na.rm = TRUE))
  y_indic  <- sapply(levels_all, function(cl) as.integer(true_chr == cl))
  brier    <- mean(rowSums((prob_mat - y_indic)^2))
  list(accuracy = accuracy, balanced_accuracy = bal_acc, brier = brier,
       vote = vote)
}

# Ordinal: categorical metrics + mean abs level distance (off-by-one is
# better than off-by-two; Hand & Till 2001 Mach Learn 45:171).
.metrics_ordinal <- function(im, true_chr, levels_ordered) {
  m <- .metrics_categorical(im, true_chr, levels_ordered)
  pos_pred  <- match(m$vote,    levels_ordered)
  pos_truth <- match(true_chr,  levels_ordered)
  m$mae_level <- mean(abs(pos_pred - pos_truth), na.rm = TRUE)
  m$vote <- NULL
  m
}

# ----------------------------------------------------------------------------
# Mean-baseline metrics — pigauto-style reference. Predicts column-mean
# for continuous/count traits and modal class for factor types. Single
# constant prediction per trait, so correlation and coverage95 are
# undefined; we emit NA for those and report rmse / mae / accuracy /
# brier just like the BACE rows so the two methods stack cleanly.
# ----------------------------------------------------------------------------
.compute_baseline_metrics <- function(masked, truth_long, types,
                                       log_traits, dataset_name) {
  rows <- lapply(names(types), function(v) {
    type <- types[[v]]
    truth_v <- truth_long[truth_long$trait == v, , drop = FALSE]
    if (nrow(truth_v) == 0L) return(NULL)
    is_log <- v %in% log_traits

    base <- data.frame(
      dataset           = dataset_name,
      method            = "mean_baseline",
      trait             = v,
      type              = type,
      scale             = if (type %in% c("continuous", "count"))
                            (if (is_log) "log" else "raw") else "categorical",
      n_hidden          = nrow(truth_v),
      rmse              = NA_real_,
      nrmse             = NA_real_,
      mae_fit           = NA_real_,
      mae_raw           = NA_real_,
      correlation       = NA_real_,
      coverage95        = NA_real_,
      accuracy          = NA_real_,
      balanced_accuracy = NA_real_,
      brier             = NA_real_,
      mae_level         = NA_real_,
      stringsAsFactors  = FALSE
    )

    if (type %in% c("continuous", "count")) {
      tv_raw  <- as.numeric(truth_v$true_value)
      tv_cmp  <- if (is_log) log1p(tv_raw) else tv_raw
      col_mu  <- mean(masked[[v]], na.rm = TRUE)
      sd_full <- stats::sd(masked[[v]], na.rm = TRUE)
      preds   <- rep(col_mu, length(tv_cmp))
      base$rmse        <- sqrt(mean((preds - tv_cmp)^2))
      base$nrmse       <- base$rmse / sd_full
      base$mae_fit     <- mean(abs(preds - tv_cmp))
      preds_raw        <- if (is_log) expm1(preds) else preds
      base$mae_raw     <- mean(abs(preds_raw - tv_raw))
      # correlation / coverage95 are NA for a constant predictor.
    } else {
      tv_chr <- as.character(truth_v$true_value)
      tab    <- table(masked[[v]], useNA = "no")
      modal  <- names(tab)[which.max(tab)]
      preds  <- rep(modal, length(tv_chr))
      base$accuracy <- mean(preds == tv_chr, na.rm = TRUE)
      classes <- unique(tv_chr)
      recalls <- vapply(classes, function(cl) {
        is_cl <- tv_chr == cl
        if (!any(is_cl)) NA_real_ else mean(preds[is_cl] == cl, na.rm = TRUE)
      }, numeric(1))
      base$balanced_accuracy <- mean(recalls, na.rm = TRUE)
      lvls <- if (is.factor(masked[[v]])) levels(masked[[v]])
              else sort(unique(c(tv_chr, modal)))
      prob_mat <- matrix(0, nrow = length(tv_chr), ncol = length(lvls),
                          dimnames = list(NULL, lvls))
      prob_mat[, modal] <- 1
      y_indic <- sapply(lvls, function(cl) as.integer(tv_chr == cl))
      base$brier <- mean(rowSums((prob_mat - y_indic)^2))
      if (type == "ordinal") {
        pos_pred  <- match(preds,  lvls)
        pos_truth <- match(tv_chr, lvls)
        base$mae_level <- mean(abs(pos_pred - pos_truth), na.rm = TRUE)
      }
    }
    base
  })
  do.call(rbind, rows)
}

# ----------------------------------------------------------------------------
# Top-level: benchmark a dataset.
# ----------------------------------------------------------------------------
benchmark_dataset <- function(
    traits_df, tree, dataset_name,
    log_traits     = character(0),
    trait_subset   = NULL,        # which cols to actually impute (default
                                  # = non-tax columns)
    # 30% MCAR mask matches Shinichi's pigauto cross-dataset bench
    # (2026-05-04 spec). Was 10% historically; the higher rate
    # stresses the imputation pipeline more and matches the
    # comparison baseline.
    cont_miss_rate = 0.30,
    cat_miss_rate  = 0.30,
    subset_n       = 2000L,
    nitt           = 20000, thin = 15, burnin = 4000,
    # n_final=10 for cloud feasibility. Pigauto's canonical
    # N_IMP=20 doubles final-imputation wall time and pushes
    # BACE's MCMCglmm chains past the GHA 5h45m budget for the
    # heavier datasets. Local runs aiming for direct pigauto
    # comparison should override n_final = 20.
    runs           = 5,    n_final = 10,
    max_attempts   = 2,    n_cores = 4L,
    skip_conv      = FALSE,
    seed           = 2026,
    out_dir_root   = "dev/benchmark_results",
    verbose        = TRUE,
    save_imputed   = FALSE) {

  stopifnot(inherits(tree, "phylo"),
            is.data.frame(traits_df),
            !is.null(rownames(traits_df)))

  # Drop higher-tax columns by default — they're known and we don't impute.
  if (is.null(trait_subset)) {
    tax_cols <- intersect(c("order_name", "family_name", "genus_name",
                            "class_name"), colnames(traits_df))
    trait_subset <- setdiff(colnames(traits_df), tax_cols)
  }

  set.seed(seed)

  # Subsample
  n_full <- nrow(traits_df)
  did_subset <- !is.na(subset_n) && subset_n < n_full
  if (did_subset) {
    keep <- sample(rownames(traits_df), subset_n)
    traits_df <- traits_df[keep, , drop = FALSE]
    tree <- ape::keep.tip(tree, keep)
  }
  tree <- .prep_tree(tree)

  # Detect types BEFORE log-transforming (log preserves class)
  types <- vapply(traits_df[, trait_subset, drop = FALSE],
                  .detect_trait_type, character(1))
  types <- types[types %in% c("continuous", "count", "ordinal",
                              "categorical", "binary")]
  if (length(types) == 0L) stop("No imputable traits in trait_subset.")

  # Log-transform requested continuous traits (use log1p to handle zeros)
  for (v in intersect(log_traits, names(types))) {
    if (types[[v]] %in% c("continuous", "count")) {
      traits_df[[v]] <- log1p(traits_df[[v]])
    }
  }

  # Build the input data frame for bace()
  masked <- traits_df[, names(types), drop = FALSE]
  masked$Species <- rownames(traits_df)
  masked <- masked[, c("Species", names(types))]
  rownames(masked) <- masked$Species

  # Apply type-aware masking
  m <- .apply_mask(masked, as.list(types), cont_miss_rate, cat_miss_rate)
  masked     <- m$masked
  truth_long <- m$truth

  # Pre-imputation phylogenetic signal
  if (verbose) cat("Computing phylogenetic signal per trait...\n")
  phylo_signal_df <- do.call(rbind, lapply(names(types), function(v) {
    type <- types[[v]]
    if (type %in% c("continuous", "count")) {
      vals <- stats::setNames(masked[[v]], masked$Species)
      s <- .phylo_signal_cont(tree, vals)
      data.frame(dataset = dataset_name, trait = v, type = type,
                 lambda = s$lambda, K = s$K, D = NA_real_,
                 stringsAsFactors = FALSE)
    } else {
      d <- .phylo_signal_cat(tree, masked[[v]], masked$Species)
      data.frame(dataset = dataset_name, trait = v, type = type,
                 lambda = NA_real_, K = NA_real_, D = d,
                 stringsAsFactors = FALSE)
    }
  }))

  # Run bace()
  fixformulas <- lapply(names(types), function(v) {
    others <- setdiff(names(types), v)
    paste(v, "~", paste(others, collapse = " + "))
  })

  if (verbose) {
    cat(sprintf("\n%s: %d species, %d traits (",
                dataset_name, nrow(masked), length(types)))
    cat(paste(sprintf("%s=%d", names(table(unname(types))),
                       unname(table(unname(types)))), collapse = ", "))
    cat(")\n")
    cat(sprintf("hidden cells: cont=%d  cat/ord/bin=%d\n",
                sum(truth_long$type %in% c("continuous", "count")),
                sum(truth_long$type %in% c("ordinal", "categorical", "binary"))))
  }

  t0 <- Sys.time()
  res <- bace(
    fixformula     = fixformulas,
    ran_phylo_form = "~1|Species",
    phylo          = tree,
    data           = masked,
    nitt           = nitt, thin = thin, burnin = burnin,
    runs           = runs, n_final = n_final,
    species        = FALSE, verbose = verbose,
    skip_conv      = skip_conv, max_attempts = max_attempts,
    n_cores        = n_cores
  )
  runtime_min <- as.numeric(round(difftime(Sys.time(), t0, units = "mins"), 1))

  # Type-aware metrics, one row per imputed trait. Tagged with
  # method = "bace" so it stacks cleanly with the mean_baseline rows.
  metric_rows <- lapply(names(types), function(v) {
    type <- types[[v]]
    truth_v <- truth_long[truth_long$trait == v, ]
    if (nrow(truth_v) == 0L) return(NULL)
    idx <- match(truth_v$species_tip, masked$Species)
    is_log <- v %in% log_traits

    base <- data.frame(
      dataset           = dataset_name,
      method            = "bace",
      trait             = v,
      type              = type,
      scale             = if (type %in% c("continuous", "count"))
                            (if (is_log) "log" else "raw") else "categorical",
      n_hidden          = nrow(truth_v),
      rmse              = NA_real_,
      nrmse             = NA_real_,
      mae_fit           = NA_real_,
      mae_raw           = NA_real_,
      correlation       = NA_real_,
      coverage95        = NA_real_,
      accuracy          = NA_real_,
      balanced_accuracy = NA_real_,
      brier             = NA_real_,
      mae_level         = NA_real_,
      stringsAsFactors  = FALSE
    )

    if (type %in% c("continuous", "count")) {
      imp_mat <- vapply(res$imputed_datasets,
                        function(d) as.numeric(d[[v]][idx]),
                        FUN.VALUE = numeric(nrow(truth_v)))
      if (!is.matrix(imp_mat)) imp_mat <- matrix(imp_mat, nrow = nrow(truth_v))
      tv_raw  <- as.numeric(truth_v$true_value)
      tv_cmp  <- if (is_log) log1p(tv_raw) else tv_raw
      sd_full <- stats::sd(masked[[v]], na.rm = TRUE)
      cm <- .metrics_continuous(imp_mat, tv_cmp, sd_full)
      post_mean      <- rowMeans(imp_mat)
      post_mean_raw  <- if (is_log) expm1(post_mean) else post_mean
      # Raw rmse on the FIT scale (matches pigauto's headline rmse column).
      base$rmse        <- sqrt(mean((post_mean - tv_cmp)^2))
      base$nrmse       <- cm$nrmse
      base$mae_fit     <- cm$mae_fit
      base$mae_raw     <- mean(abs(post_mean_raw - tv_raw))
      base$correlation <- cm$correlation
      base$coverage95  <- cm$coverage95
    } else {
      im <- vapply(res$imputed_datasets,
                   function(d) as.character(d[[v]][idx]),
                   FUN.VALUE = character(nrow(truth_v)))
      if (!is.matrix(im)) im <- matrix(im, nrow = nrow(truth_v))
      lvls   <- if (is.factor(masked[[v]])) levels(masked[[v]])
                else sort(unique(as.character(truth_v$true_value)))
      tv_chr <- as.character(truth_v$true_value)
      cm <- if (type == "ordinal") .metrics_ordinal(im, tv_chr, lvls)
            else                    .metrics_categorical(im, tv_chr, lvls)
      base$accuracy          <- cm$accuracy
      base$balanced_accuracy <- cm$balanced_accuracy
      base$brier             <- cm$brier
      if (!is.null(cm$mae_level)) base$mae_level <- cm$mae_level
    }
    base
  })
  metrics_bace <- do.call(rbind, metric_rows)

  # Mean-baseline reference (pigauto-aligned) -- stacked alongside bace.
  metrics_baseline <- .compute_baseline_metrics(
    masked, truth_long, types, log_traits, dataset_name)

  metrics_df <- rbind(metrics_bace, metrics_baseline)

  # Output
  date_str  <- format(Sys.Date(), "%Y%m%d")
  scope_tag <- if (did_subset) paste0("n", subset_n) else "full"
  out_dir <- file.path(out_dir_root, dataset_name,
                       paste0("run_", scope_tag, "_", date_str))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  utils::write.csv(metrics_df,
                   file.path(out_dir, "summary_metrics.csv"),
                   row.names = FALSE)

  # Pigauto-style tidy long format: (dataset, method, trait, type,
  # metric, value). Each (method, trait) row in the wide frame
  # spreads to N rows, one per metric. Easier to stack across
  # methods / runs / datasets for cross-method comparison.
  metric_cols <- c("rmse", "nrmse", "mae_fit", "mae_raw",
                   "correlation", "coverage95", "accuracy",
                   "balanced_accuracy", "brier", "mae_level")
  metrics_long <- do.call(rbind, lapply(metric_cols, function(m) {
    d <- metrics_df[, c("dataset", "method", "trait", "type",
                        "scale", "n_hidden", m)]
    colnames(d)[ncol(d)] <- "value"
    d$metric <- m
    d <- d[!is.na(d$value), , drop = FALSE]
    d[, c("dataset","method","trait","type","scale","n_hidden",
          "metric","value")]
  }))
  utils::write.csv(metrics_long,
                   file.path(out_dir, "summary_metrics_long.csv"),
                   row.names = FALSE)
  utils::write.csv(phylo_signal_df,
                   file.path(out_dir, "phylo_signal.csv"),
                   row.names = FALSE)
  utils::write.csv(
    data.frame(dataset      = dataset_name,
               n_species    = nrow(masked),
               n_traits     = length(types),
               runtime_min  = runtime_min,
               converged    = res$converged,
               n_attempts   = res$n_attempts,
               nitt = nitt, thin = thin, burnin = burnin,
               runs = runs, n_final = n_final),
    file.path(out_dir, "run_info.csv"), row.names = FALSE)

  if (save_imputed) saveRDS(res, file.path(out_dir, "bace_result.rds"))

  if (verbose) {
    cat(sprintf("\n=== %s metrics ===\n", dataset_name))
    print(metrics_df, digits = 3, row.names = FALSE)
    cat(sprintf("\nRuntime: %.1f min   converged: %s   out: %s\n",
                runtime_min, res$converged, out_dir))
  }

  invisible(list(
    metrics      = metrics_df,
    phylo_signal = phylo_signal_df,
    runtime_min  = runtime_min,
    res          = res,
    out_dir      = out_dir
  ))
}
