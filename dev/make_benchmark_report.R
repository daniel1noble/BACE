# =============================================================================
# make_benchmark_report.R
# =============================================================================
# Aggregate per-dataset benchmark output into a single cross-dataset bundle:
#   all_metrics.csv       long-format per-trait metric table
#   all_phylo_signal.csv  pre-imputation Pagel lambda / Blomberg K / D
#   runtimes.csv          per-dataset runtime + convergence
#   report.md             human-readable summary (renderable on GitHub)
#   report.html           same content as a self-contained HTML page
#
# Inputs: dev/benchmark_results/_artifacts/bench-<dataset>/<dataset>/run_*/
#   summary_metrics.csv, phylo_signal.csv, run_info.csv
#
# This script is invoked by the GitHub Actions "Cross-dataset benchmark"
# workflow's aggregate job. It also works standalone after a local sweep
# where per-dataset bundles live under dev/benchmark_results/<dataset>/.
# =============================================================================

suppressPackageStartupMessages({
  library(knitr)
})

# ---- 0. Locate per-dataset run bundles --------------------------------------

# Two layouts to support:
#   (a) GHA download-artifact  ->  dev/benchmark_results/_artifacts/bench-<ds>/...
#   (b) local sweep            ->  dev/benchmark_results/<ds>/run_*/
search_roots <- c(
  "dev/benchmark_results/_artifacts",
  "dev/benchmark_results"
)
search_roots <- search_roots[dir.exists(search_roots)]

metric_files <- unlist(lapply(search_roots, function(root)
  list.files(root, pattern = "^summary_metrics\\.csv$",
             recursive = TRUE, full.names = TRUE)))
signal_files <- unlist(lapply(search_roots, function(root)
  list.files(root, pattern = "^phylo_signal\\.csv$",
             recursive = TRUE, full.names = TRUE)))
info_files   <- unlist(lapply(search_roots, function(root)
  list.files(root, pattern = "^run_info\\.csv$",
             recursive = TRUE, full.names = TRUE)))

# Avoid double-counting if the cross_dataset/ folder is in search_roots
metric_files <- metric_files[!grepl("/cross_dataset/", metric_files)]
signal_files <- signal_files[!grepl("/cross_dataset/", signal_files)]
info_files   <- info_files  [!grepl("/cross_dataset/", info_files)]

if (length(metric_files) == 0L) {
  stop("No summary_metrics.csv files found under: ",
       paste(search_roots, collapse = ", "))
}

cat(sprintf("Found %d metric files, %d signal files, %d info files.\n",
            length(metric_files), length(signal_files), length(info_files)))

read_safe <- function(f, required_col = "dataset") {
  d <- tryCatch(read.csv(f, stringsAsFactors = FALSE),
                error = function(e) {
                  warning("Failed to read ", f, ": ", conditionMessage(e))
                  NULL
                })
  if (is.null(d)) return(NULL)
  # Filter out pre-engine output bundles (their schema lacks a `dataset`
  # column; rows would silently merge in with NA dataset).
  if (!is.null(required_col) && !(required_col %in% colnames(d))) {
    message("  skipping pre-engine bundle: ", f)
    return(NULL)
  }
  d
}

bind_safe <- function(tables) {
  ok <- !vapply(tables, is.null, logical(1))
  if (!any(ok)) return(NULL)
  tables <- tables[ok]
  # Align columns so old-format bundles or partial outputs do not
  # break rbind. Union of all columns; missing columns filled with NA.
  all_cols <- unique(unlist(lapply(tables, colnames)))
  aligned  <- lapply(tables, function(d) {
    missing <- setdiff(all_cols, colnames(d))
    for (m in missing) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  do.call(rbind, aligned)
}

all_metrics <- bind_safe(lapply(metric_files, read_safe, "dataset"))
all_signal  <- bind_safe(lapply(signal_files, read_safe, "dataset"))
all_info    <- bind_safe(lapply(info_files,   read_safe, "dataset"))

# ---- 1. Output directory ----------------------------------------------------

date_str <- format(Sys.Date(), "%Y%m%d")
out_dir  <- file.path("dev", "benchmark_results", "cross_dataset",
                       paste0("run_", date_str))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!is.null(all_metrics))
  write.csv(all_metrics, file.path(out_dir, "all_metrics.csv"),
            row.names = FALSE)
if (!is.null(all_signal))
  write.csv(all_signal,  file.path(out_dir, "all_phylo_signal.csv"),
            row.names = FALSE)
if (!is.null(all_info))
  write.csv(all_info,    file.path(out_dir, "runtimes.csv"),
            row.names = FALSE)

# ---- 2. Build markdown report -----------------------------------------------

git_sha <- tryCatch(
  substr(system("git rev-parse HEAD", intern = TRUE,
                 ignore.stderr = TRUE), 1, 8),
  error   = function(e) NA_character_,
  warning = function(w) NA_character_)

bace_version <- tryCatch(
  as.character(utils::packageVersion("BACE")),
  error = function(e) NA_character_)

md <- character()
add <- function(...) md <<- c(md, sprintf(...))
nl  <- function() md <<- c(md, "")

add("# BACE cross-dataset benchmark — %s", date_str)
nl()
add("- BACE commit: `%s`", if (is.na(git_sha)) "unknown" else git_sha)
add("- BACE version: %s", if (is.na(bace_version)) "unknown" else bace_version)
add("- Datasets: %s",
    if (!is.null(all_metrics)) paste(sort(unique(all_metrics$dataset)),
                                       collapse = ", ") else "none")
nl()

# ---- Failure callouts -------------------------------------------------------

# Cloud benchmark suite. avonet/pantheria/amphibio/bien are aligned
# with Shinichi's pigauto cross-dataset bench (2026-05-04 spec);
# globtherm and leptraits are BACE-only (not in pigauto). fishbase
# TODO: needs data-raw/make_fishbase.R.
expected_datasets <- c("avonet", "pantheria", "amphibio", "bien",
                       "globtherm", "leptraits")
present_datasets  <- if (!is.null(all_metrics))
  sort(unique(all_metrics$dataset)) else character()
missing_datasets  <- setdiff(expected_datasets, present_datasets)
if (length(missing_datasets)) {
  add("> **Warning:** Missing output for: %s. Check the workflow run logs.",
      paste(missing_datasets, collapse = ", "))
  nl()
}

# ---- Per-dataset summary ----------------------------------------------------

if (!is.null(all_metrics)) {
  add("## Per-dataset summary")
  nl()
  ds_summary <- do.call(rbind, lapply(split(all_metrics, all_metrics$dataset),
    function(d) {
      cont <- d[d$type %in% c("continuous", "count"), , drop = FALSE]
      catg <- d[d$type %in% c("categorical", "binary", "ordinal"),
                , drop = FALSE]
      data.frame(
        dataset       = d$dataset[1],
        n_traits      = nrow(d),
        mean_NRMSE    = if (nrow(cont)) round(mean(cont$nrmse,
                                                    na.rm = TRUE), 3)
                        else NA_real_,
        mean_cor      = if (nrow(cont)) round(mean(cont$correlation,
                                                    na.rm = TRUE), 3)
                        else NA_real_,
        mean_cov95    = if (nrow(cont)) round(mean(cont$coverage95,
                                                    na.rm = TRUE), 3)
                        else NA_real_,
        mean_accuracy = if (nrow(catg)) round(mean(catg$accuracy,
                                                    na.rm = TRUE), 3)
                        else NA_real_,
        mean_bal_acc  = if (nrow(catg)) round(mean(catg$balanced_accuracy,
                                                    na.rm = TRUE), 3)
                        else NA_real_,
        mean_brier    = if (nrow(catg)) round(mean(catg$brier,
                                                    na.rm = TRUE), 3)
                        else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  rownames(ds_summary) <- NULL
  md <- c(md, knitr::kable(ds_summary, format = "markdown"))
  nl()
}

# ---- Continuous trait detail ------------------------------------------------

if (!is.null(all_metrics)) {
  cont <- all_metrics[all_metrics$type %in% c("continuous","count"), ,
                       drop = FALSE]
  if (nrow(cont)) {
    add("## Continuous + count traits")
    nl()
    cols <- intersect(c("dataset","trait","type","scale","n_hidden",
                        "nrmse","mae_fit","mae_raw","correlation",
                        "coverage95"), colnames(cont))
    md <- c(md, knitr::kable(cont[, cols], format = "markdown",
                              digits = 3, row.names = FALSE))
    nl()
  }

  catg <- all_metrics[all_metrics$type %in% c("categorical","binary","ordinal"),
                       , drop = FALSE]
  if (nrow(catg)) {
    add("## Categorical / binary / ordinal traits")
    nl()
    cols <- intersect(c("dataset","trait","type","n_hidden",
                        "accuracy","balanced_accuracy","brier","mae_level"),
                      colnames(catg))
    md <- c(md, knitr::kable(catg[, cols], format = "markdown",
                              digits = 3, row.names = FALSE))
    nl()
  }
}

# ---- Phylogenetic signal context --------------------------------------------

if (!is.null(all_signal)) {
  add("## Pre-imputation phylogenetic signal")
  nl()
  add("Reference values per trait, computed before BACE runs. Pagel λ and Blomberg K for continuous; Fritz-Purvis D (OVR mean) for categorical.")
  nl()
  cols <- intersect(c("dataset","trait","type","lambda","K","D"),
                    colnames(all_signal))
  md <- c(md, knitr::kable(all_signal[, cols], format = "markdown",
                            digits = 3, row.names = FALSE))
  nl()
}

# ---- Run info ---------------------------------------------------------------

if (!is.null(all_info)) {
  add("## Runtime + MCMC config")
  nl()
  cols <- intersect(c("dataset","n_species","n_traits","runtime_min",
                      "converged","n_attempts","nitt","thin","burnin",
                      "runs","n_final"), colnames(all_info))
  md <- c(md, knitr::kable(all_info[, cols], format = "markdown",
                            row.names = FALSE))
  nl()
}

writeLines(md, file.path(out_dir, "report.md"))

# ---- 3. Render HTML self-contained -----------------------------------------

# Wrap the markdown in a minimal self-contained HTML using knit2html-style
# rendering. If rmarkdown + pandoc are available, prefer pandoc for
# better-looking tables; otherwise fall back to a hand-rolled HTML wrapper
# so the report still produces something useful.
md_path   <- file.path(out_dir, "report.md")
html_path <- file.path(out_dir, "report.html")

rendered <- FALSE
if (requireNamespace("rmarkdown", quietly = TRUE) &&
    rmarkdown::pandoc_available()) {
  tryCatch({
    rmarkdown::render(md_path,
                      output_file = "report.html",
                      output_dir  = out_dir,
                      output_format = rmarkdown::html_document(theme = "cosmo"),
                      quiet = TRUE)
    rendered <- TRUE
  }, error = function(e) {
    warning("rmarkdown render failed: ", conditionMessage(e),
            ". Falling back to plain HTML.")
  })
}

if (!rendered) {
  body <- paste(readLines(md_path), collapse = "\n")
  body_html <- gsub("&", "&amp;",  body, fixed = TRUE)
  body_html <- gsub("<", "&lt;",   body_html, fixed = TRUE)
  body_html <- gsub(">", "&gt;",   body_html, fixed = TRUE)
  html <- sprintf(
    "<!doctype html><meta charset=\"utf-8\"><title>BACE benchmark %s</title>%s%s",
    date_str,
    "<style>body{font-family:sans-serif;max-width:1200px;margin:2em auto;padding:0 1em;line-height:1.4} pre{background:#f4f4f4;padding:0.5em;overflow-x:auto} table{border-collapse:collapse} th,td{border:1px solid #ddd;padding:4px 8px}</style>",
    sprintf("<pre>%s</pre>", body_html))
  writeLines(html, html_path)
}

cat(sprintf("\nReport bundle written to: %s\n", out_dir))
cat("Files:\n  ", paste(list.files(out_dir), collapse = "\n  "), "\n")
