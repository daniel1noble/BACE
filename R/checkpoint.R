# Checkpoint/resume support.
#
# Long BACE runs take many hours (22-28 h in the published case studies), so
# losing them to a crashed session is expensive. With `checkpoint_dir` set,
# the convergence phase saves its finished object and each final imputation
# run saves its own file as soon as it completes. Calling the same function
# again with the same arguments and the same directory loads what is already
# there and computes only what is missing.
#
# Two rules make that safe:
#
#  * A "manifest" records the settings the saved work was produced under. If
#    the new call disagrees on anything that would change the results, the
#    resume is REFUSED and the differing settings are named. Results from two
#    different setups are never silently mixed.
#  * The final phase requires `seed`. Each final run draws from its own
#    random-number stream derived from (seed, run number), so a resumed run is
#    identical to one that was never interrupted. Without a seed that
#    guarantee does not hold, so we ask for one rather than quietly producing
#    something different.
#
# The convergence phase does not take a seed (it is unseeded by design), which
# is exactly why its result is SAVED rather than recomputed on resume.

# --- paths -----------------------------------------------------------------

#' @keywords internal
#' @noRd
.ckpt_manifest_file <- function(dir) file.path(dir, "final_manifest.rds")

#' Is this object a completed final imputation run?
#'
#' A file can be readable and still not be one of our runs (a leftover from
#' another analysis, or something the user put there). Accepting the wrong
#' shape would quietly drop an imputation from the pooling instead of failing.
#' @keywords internal
#' @noRd
.ckpt_is_run <- function(x) {
  # "scaling" is required from 1.0.2. A checkpoint written by an earlier build
  # has no z-scoring constants in it, so accepting one would silently leave
  # that run's constants NULL and make its coefficients unrecoverable -- the
  # exact gap the constants were added to close. Such a file fails this check
  # and the run is recomputed, which is the honest outcome.
  is.list(x) && !is.data.frame(x) &&
    all(c("models", "dataset", "prob_preds", "ovr_models", "scaling") %in%
          names(x))
}

#' @keywords internal
#' @noRd
.ckpt_run_file <- function(dir, run) {
  file.path(dir, sprintf("final_run%04d.rds", run))
}

#' Files holding a saved convergence chain
#'
#' The file name carries the number of chained iterations, so the convergence
#' retry ladder (which re-runs with more iterations when the chain has not
#' settled) saves each attempt separately instead of colliding with, or
#' refusing, the previous one. Resuming an interrupted \code{bace()} then
#' reloads the attempts that finished and only runs the ones that did not.
#' @keywords internal
#' @noRd
.ckpt_chain_files <- function(dir, runs, tag = NULL) {
  stem <- if (is.null(tag)) sprintf("chain_runs%04d", as.integer(runs)) else
    sprintf("chain_%s_runs%04d", tag, as.integer(runs))
  list(object   = file.path(dir, paste0(stem, ".rds")),
       manifest = file.path(dir, paste0(stem, "_manifest.rds")))
}

#' Create (if needed) and sanity-check a checkpoint directory
#' @keywords internal
#' @noRd
.ckpt_dir_prepare <- function(dir, context = "bace") {
  if (!is.character(dir) || length(dir) != 1L || is.na(dir) || !nzchar(dir)) {
    stop(context, ": checkpoint_dir must be a single directory path.",
         call. = FALSE)
  }
  if (!dir.exists(dir)) {
    ok <- dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(dir)) {
      stop(context, ": could not create checkpoint_dir '", dir, "'.",
           call. = FALSE)
    }
  }
  probe <- file.path(dir, ".bace_write_probe")
  ok <- tryCatch({ saveRDS(TRUE, probe); TRUE }, error = function(e) FALSE,
                 warning = function(w) FALSE)
  unlink(probe)
  if (!ok) {
    stop(context, ": checkpoint_dir '", dir, "' is not writable.",
         call. = FALSE)
  }
  invisible(normalizePath(dir))
}

# --- atomic write ----------------------------------------------------------

#' Write a checkpoint file atomically
#'
#' Saves to a temporary file in the SAME directory and then renames it, so a
#' crash mid-write can never leave a half-written checkpoint that a later
#' resume would trust. The rename's return value is checked: it returns FALSE
#' (rather than throwing) when it cannot move the file, and an unchecked
#' failure would silently lose the run. Several parallel workers may write
#' here at once, hence the unique temporary name.
#'
#' Compression is off by default: on MCMC output it roughly halves the file
#' but costs about ten times the writing time. Against hours of model fitting
#' that is a fraction of a second per run, so turn it on if disk space is
#' tighter than patience.
#' @keywords internal
#' @noRd
.ckpt_write <- function(object, path, compress = FALSE) {
  tmp <- tempfile(pattern = "bace_ckpt", tmpdir = dirname(path),
                  fileext = ".tmp")
  saveRDS(object, tmp, compress = compress)
  if (!isTRUE(file.rename(tmp, path))) {
    unlink(tmp)
    stop("Could not move the checkpoint into place ('", basename(path),
         "'). The run was NOT saved. Check that the checkpoint directory is ",
         "writable and on the same filesystem as the temporary file.",
         call. = FALSE)
  }
  invisible(path)
}

#' Remove temporary files left behind by an interrupted write
#' @keywords internal
#' @noRd
.ckpt_sweep_stale <- function(dir, verbose = TRUE) {
  stale <- list.files(dir, pattern = "^bace_ckpt.*[.]tmp$", full.names = TRUE)
  if (length(stale) > 0L) {
    unlink(stale)
    if (verbose) {
      message("checkpoint: removed ", length(stale),
              " leftover temporary file(s) from an interrupted write.")
    }
  }
  invisible(length(stale))
}

# --- manifests -------------------------------------------------------------

#' Text form of a formula or formula-like input
#'
#' Formulas are compared as text, never hashed as objects: serialize() follows
#' a formula's environment, so two formulas with identical text but built in
#' different places would look different (and could drag megabytes of
#' unrelated data into the hash).
#' @keywords internal
#' @noRd
.ckpt_formula_text <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !inherits(x, "formula")) {
    return(vapply(x, .ckpt_formula_text, character(1L)))
  }
  paste(deparse(stats::as.formula(x)), collapse = " ")
}

#' Compare two manifests and name the settings that differ
#'
#' `mutable` fields are allowed to differ without refusing the resume: the
#' convergence retry ladder deliberately re-runs with more iterations, so
#' `runs` and `attempt` change between attempts of the same analysis.
#' @return character vector of differing field names (empty = compatible)
#' @keywords internal
#' @noRd
.ckpt_diff <- function(old, new, mutable = character(0)) {
  fields <- setdiff(union(names(old), names(new)), mutable)
  differs <- vapply(fields, function(f) !isTRUE(all.equal(old[[f]], new[[f]])),
                    logical(1L))
  fields[differs]
}

#' Load a manifest, refusing the resume if the settings disagree
#'
#' @param path manifest file
#' @param manifest the manifest for the call being made now
#' @param what human-readable phase name for messages
#' @return TRUE if a compatible manifest already existed (i.e. this is a
#'   resume), FALSE if the manifest was newly written
#' @keywords internal
#' @noRd
.ckpt_manifest_check <- function(path, manifest, what, mutable = character(0),
                                 verbose = TRUE, guard_files = character(0)) {
  if (!file.exists(path)) {
    # Saved work with no manifest beside it cannot be shown to belong to this
    # analysis. Refuse rather than either trusting it or deleting it.
    if (length(guard_files) > 0L) {
      stop("checkpoint: this directory holds ", length(guard_files),
           " saved file(s) but no ", what, " manifest, so they cannot be ",
           "matched to this analysis. Use a different (empty) checkpoint_dir, ",
           "or delete the saved files to start this phase again.",
           call. = FALSE)
    }
    .ckpt_write(manifest, path, compress = TRUE)
    return(FALSE)
  }
  old <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(old)) {
    # Without a readable manifest there is no way to show that saved work
    # belongs to this analysis, so refuse rather than guess.
    if (length(guard_files) > 0L) {
      stop("checkpoint: the ", what, " manifest in this directory cannot be ",
           "read, so the ", length(guard_files), " saved file(s) beside it ",
           "cannot be matched to this analysis. Use a different (empty) ",
           "checkpoint_dir, or delete the saved files to start this phase ",
           "again.", call. = FALSE)
    }
    if (verbose) {
      message("checkpoint: the existing ", what, " manifest could not be read;",
              " starting this phase from scratch.")
    }
    .ckpt_write(manifest, path, compress = TRUE)
    return(FALSE)
  }
  bad <- .ckpt_diff(old, manifest, mutable = mutable)
  if (length(bad) > 0L) {
    hint <- ""
    if ("bace_version" %in% bad) {
      hint <- paste0(
        " The saved work was produced by BACE ",
        if (is.null(old$bace_version)) "(unknown)" else old$bace_version,
        ". To finish an interrupted run, install that version (for example ",
        "with remotes::install_github(\"daniel1noble/BACE@<tag or sha>\")); ",
        "otherwise use a new, empty checkpoint_dir.")
    }
    stop("checkpoint: the saved ", what, " work in this directory was ",
         "produced with different settings, so resuming it could mix ",
         "incompatible results. Differing: ",
         paste(bad, collapse = ", "), ".", hint,
         " Use a different (empty) checkpoint_dir to start fresh.",
         call. = FALSE)
  }
  # refresh mutable fields (e.g. the retry ladder's attempt / runs)
  if (length(mutable) > 0L) {
    for (f in mutable) old[[f]] <- manifest[[f]]
    .ckpt_write(old, path, compress = TRUE)
  }
  TRUE
}

#' Fields every manifest carries
#'
#' The R version is its own field so that a refusal after an R upgrade names
#' the R version rather than blaming the data.
#' @keywords internal
#' @noRd
.ckpt_common <- function() {
  list(bace_version = as.character(utils::packageVersion("BACE")),
       r_version    = paste0(R.version$major, ".", R.version$minor),
       gelman       = if (is.null(getOption("BACE.gelman")))
                        bace_option_defaults()$gelman else
                        getOption("BACE.gelman"))
}

#' Evaluate an expression under a fixed random seed, restoring the caller's state
#'
#' Used for the pooled-posterior subsample: a resumed run skips the chained
#' phase, so the ambient random stream sits in a different place than it would
#' after an uninterrupted run. Seeding the subsample makes it depend on the
#' user's `seed` instead, which is what makes a resumed run match.
#' @keywords internal
#' @noRd
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old)) {
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    } else {
      assign(".Random.seed", old, envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}
