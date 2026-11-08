# bace_validate(): check how well BACE imputes YOUR data.
#
# The idea is simple. Take cells where you DO know the answer, hide a share of
# them, let BACE impute them, and compare the imputations to the hidden truth.
# If BACE recovers values it has never seen, you have grounds to trust the
# values it fills in for the cells that are truly missing.
#
# Design notes live in .agents/plans/next-features-2026-11.md (feature 2).
# zRMSE, its baseline, interval coverage and accuracy (including fractional
# credit for tied votes) match the BACE manuscript exactly, so those numbers
# are directly comparable to the published simulations; balanced accuracy and
# the Brier score are additional metrics the manuscript does not report.

# ---------------------------------------------------------------------------
# Exported helper: per-cell intervals for imputed values
# ---------------------------------------------------------------------------

#' Per-cell intervals for imputed values
#'
#' For every originally missing cell, this gives the mean imputed value and a
#' 95 percent interval built from the spread across the \code{n_final}
#' imputations.
#'
#' The interval is the standard multiple-imputation one (Rubin 1987):
#' mean +/- t * sqrt((1 + 1/M) * v), where M is the number of imputations, v
#' is the variance of the cell's M imputed values, and t is the 97.5 percent
#' point of a t distribution with M - 1 degrees of freedom. For normally
#' spread imputations this interval is exact at any M, and in practice it
#' stays close to its stated level from M of about 20. Taking the 2.5 and
#' 97.5 percent quantiles of the imputed values instead undercovers (about
#' 92 percent at M = 50).
#'
#' Intervals are only meaningful for numeric variables (continuous or count);
#' rows for factor variables are left out.
#'
#' @param object A result from \code{\link{bace}} or
#'   \code{\link{bace_final_imp}}.
#' @param variable Optional: one variable name. Default: all numeric variables
#'   that had missing data.
#' @param level Interval level. Default 0.95.
#' @return A data frame with one row per imputed cell: \code{row} (row number
#'   in the original data), \code{species} (the cluster column's value),
#'   \code{variable}, \code{mean}, \code{sd}, \code{lower}, \code{upper}, and
#'   \code{n_imputations}. If no numeric variable had missing data, a zero-row
#'   data frame is returned with a message. With a single imputation
#'   (\code{n_final = 1}) the spread is undefined: the mean is returned and
#'   the bounds are NA, with a message.
#' @references Rubin, D. B. (1987) Multiple Imputation for Nonresponse in
#'   Surveys. Wiley.
#' @seealso \code{\link{get_imputed_data}}, \code{\link{bace_validate}}
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- sim_bace(response_type = "gaussian", predictor_types = "gaussian",
#'                 phylo_signal = c(0.6, 0.6), missingness = c(0.2, 0),
#'                 n_cases = 40, n_species = 40)
#' fit <- bace(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
#'             phylo = sim$tree, data = sim$data, runs = 2, nitt = 400,
#'             burnin = 100, thin = 2, n_final = 5, skip_conv = TRUE,
#'             verbose = FALSE)
#'
#' # Mean and 95% interval for every imputed cell
#' head(get_imputation_intervals(fit))
#' }
#' @export
get_imputation_intervals <- function(object, variable = NULL, level = 0.95) {
  fin <- if (inherits(object, "bace_complete")) object$final_results else object
  if (!inherits(fin, "bace_final")) {
    stop("`object` must come from bace() or bace_final_imp().", call. = FALSE)
  }
  datasets <- fin$all_datasets
  miss     <- fin$miss_dat
  M        <- length(datasets)
  cluster  <- fin$phylo_ran[["cluster"]]

  empty <- data.frame(row = integer(0), species = character(0),
                      variable = character(0), mean = numeric(0),
                      sd = numeric(0), lower = numeric(0), upper = numeric(0),
                      n_imputations = integer(0))

  vars <- unique(miss$colname)
  # keep numeric variables only (intervals make no sense for classes)
  vars <- vars[vapply(vars, function(v)
    is.numeric(datasets[[1]][[v]]), logical(1))]
  if (!is.null(variable)) {
    if (!variable %in% vars) {
      stop("'", variable, "' is not a numeric variable with missing data.",
           call. = FALSE)
    }
    vars <- variable
  }
  if (length(vars) == 0L) {
    message("No numeric variables had missing data; nothing to interval.")
    return(empty)
  }
  if (M < 2L) {
    message("Only one imputation: the spread between imputations is ",
            "undefined, so bounds are NA. Use n_final >= 20 for intervals.")
  }

  t_mult <- if (M >= 2L) stats::qt(1 - (1 - level) / 2, df = M - 1) else NA_real_
  out <- lapply(vars, function(v) {
    rows <- miss$row[miss$colname == v]
    vals <- vapply(datasets, function(d) as.numeric(d[rows, v]),
                   numeric(length(rows)))
    vals <- matrix(vals, nrow = length(rows))   # cells x imputations
    m  <- rowMeans(vals)
    s  <- if (M >= 2L) apply(vals, 1L, stats::sd) else rep(NA_real_, length(rows))
    hw <- t_mult * sqrt(1 + 1 / M) * s          # half-width
    data.frame(row = rows,
               species = as.character(datasets[[1]][rows, cluster]),
               variable = v, mean = m, sd = s,
               lower = m - hw, upper = m + hw,
               n_imputations = M, row.names = NULL)
  })
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# The main function
# ---------------------------------------------------------------------------

#' Validate BACE imputation on your own data
#'
#' Hides a share of the OBSERVED cells, runs the full BACE pipeline on the
#' masked data, and scores the imputations against the hidden truth. One call
#' answers: "how well does the imputation model predict values it has never
#' seen, in this dataset?"
#'
#' The natural way to use it is on a finished run:
#' \code{bace_validate(fit)} takes the result of \code{\link{bace}} and reuses
#' its data, formulas, tree(s) and model settings, so nothing is re-typed.
#' Note what this does under the hood: a finished run cannot be scored
#' directly -- its imputations were made with every observed cell visible, so
#' there is no held-out truth to compare against. The check therefore always
#' RE-RUNS the pipeline on freshly masked data (one full run per repetition);
#' what the \code{bace_complete} method saves you is the re-typing, not the
#' compute. You can equally run it before any full analysis by giving the
#' data and formulas directly, exactly as in \code{bace()}.
#'
#' Only variables that BACE imputes (the left-hand side of some formula) are
#' masked and scored; a complete covariate that only ever appears as a
#' predictor is left untouched.
#'
#' @details
#' \strong{What the scores mean.}
#' \itemize{
#'   \item \strong{zRMSE} (continuous and count variables): the error of the
#'     mean imputed value, in units of the variable's spread. 0 means the
#'     hidden values were recovered exactly; 1 means no better than filling
#'     every cell with the variable's mean. The mean-filling baseline is
#'     printed next to it.
#'   \item \strong{coverage}: how often the 95 percent interval (see
#'     \code{\link{get_imputation_intervals}}) contained the hidden value.
#'     Should be close to 0.95. For count variables the attainable level is a
#'     little lower because counts are whole numbers; the printed expectation
#'     is simulated for the variable's own rates, and a well-calibrated
#'     imputer typically lands within about 0.01 of it.
#'   \item \strong{accuracy} (factor variables): how often the most frequent
#'     imputed class was the hidden class, next to the baseline of always
#'     guessing the most common class among the cells left observed. Tied
#'     votes get fractional credit (a 2-way tie containing the truth scores
#'     0.5), matching the BACE manuscript. \strong{Brier} is the mean squared
#'     error of the imputed class probabilities (0 is perfect; 2 is a
#'     confidently wrong guess).
#'   \item For \strong{count} variables two reference lines are printed:
#'     the hard floor from counting noise alone (no imputer can score below
#'     it) and the level a well-calibrated imputer should sit near. A zRMSE
#'     near those lines means the imputer is extracting all the signal there
#'     is -- not failing.
#' }
#'
#' \strong{What this check can and cannot tell you.} Masking is random among
#' the observed cells, so the check measures how well the model predicts in
#' the region where you have data. That is necessary, but not sufficient, for
#' trusting the imputations of the truly missing cells: if your data are
#' missing not-at-random (say, small or rare species were never measured),
#' the truly missing cells live outside the tested region and no masking
#' check can see that. The approach is inspired by Amelia's
#' \code{overimpute()} (Honaker, King & Blackwell 2011); note Amelia hides
#' one cell at a time, while \code{bace_validate()} hides a share of them at
#' once -- a deliberately harder test, so do not expect identical numbers
#' from the two tools.
#'
#' \strong{Cost.} Each repetition runs the complete BACE pipeline once (and
#' the convergence retry can multiply that -- see \code{max_attempts} and
#' \code{skip_conv}). For a quick first look, use lighter MCMC settings and
#' remember the scores then reflect those settings. With a
#' \code{bace_complete} input the stored settings are reused; to choose
#' lighter ones, call \code{bace_validate()} with the data and formulas
#' directly.
#'
#' @param x Either a finished run (a \code{bace_complete} object from
#'   \code{\link{bace}}) -- its data, formulas, tree(s), number of chained
#'   iterations (\code{runs}) and model settings are reused -- or the
#'   \code{fixformula} (a single string, a formula, or a list of formulas),
#'   with the data and tree given alongside, exactly as in \code{bace()}.
#' @param ran_phylo_form,phylo,data,nitt,thin,burnin,final_nitt,final_thin,final_burnin,runs,n_final,species,verbose,max_attempts,skip_conv,n_cores,nitt_cat_mult,ovr_categorical,priors,sweeps
#'   Passed to \code{\link{bace}}, with the same meaning and defaults. For a
#'   \code{bace_complete} input these are taken from the finished run; only
#'   \code{verbose}, \code{n_cores}, \code{max_attempts} and \code{skip_conv}
#'   remain free (they are run-time choices, not model settings), and any
#'   other named argument is refused with a warning rather than silently
#'   ignored.
#' @param ... Passed between methods.
#' @param prop Share of each imputed variable's OBSERVED cells to hide.
#'   Default 0.1. Variables with fewer than 10 observed cells are skipped
#'   with a warning; masking always leaves at least 2 observed cells, and for
#'   factors it never empties an observed class.
#' @param reps How many times to repeat the whole check with a fresh random
#'   mask. Default 1. With 2 or more, the summary reports each score as a
#'   mean with a standard error of that mean across repetitions -- the honest
#'   uncertainty to quote.
#' @param scale Optional named character vector for continuous variables
#'   recorded as proportions, e.g. \code{c(prop_trait = "logit")}: the
#'   variable is scored on the logit scale (values first clipped to the range
#'   0.0001 to 0.9999), matching the BACE manuscript. Default: all variables
#'   scored on their raw scale. Unknown names or values are warned about.
#' @param seed Optional integer (any value \code{set.seed} accepts). Fixes
#'   the masks exactly (each repetition's mask is derived from the seed and
#'   the repetition number) and makes the final-phase draws reproducible
#'   given the converged data. The scores still vary a little between
#'   identical calls, because the convergence phase that precedes the final
#'   draws is deliberately unseeded (see \code{\link{bace}}). If NULL, mask
#'   seeds are drawn at random and recorded in the result's
#'   \code{$settings}.
#' @return An object of class \code{bace_validation}: \code{$table} (one row
#'   per variable per repetition, all scores and baselines),
#'   \code{$summary} (scores averaged over repetitions, with a between-rep
#'   standard error when reps > 1), \code{$masked_cells} (which cells were
#'   hidden, with their true values), \code{$settings}, and \code{$n_final}.
#'   Print it for a plain-language report.
#' @references Honaker, J., King, G. & Blackwell, M. (2011) Amelia II: A
#'   program for missing data. Journal of Statistical Software 45(7), 1-47.
#' @seealso \code{\link{bace}}, \code{\link{get_imputation_intervals}}
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- sim_bace(response_type = "gaussian", predictor_types = "gaussian",
#'                 phylo_signal = c(0.6, 0.6), missingness = c(0.2, 0),
#'                 n_cases = 40, n_species = 40)
#'
#' # A (deliberately fast) run, then validate it: the finished object's
#' # data, tree and settings are reused -- nothing re-typed.
#' fit <- bace(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
#'             phylo = sim$tree, data = sim$data, runs = 2, nitt = 400,
#'             burnin = 100, thin = 2, n_final = 5, skip_conv = TRUE,
#'             verbose = FALSE)
#' v <- bace_validate(fit, prop = 0.2, seed = 1, verbose = FALSE)
#' v
#'
#' # The same check straight from data and formulas (no prior run needed):
#' v2 <- bace_validate("y ~ x1", ran_phylo_form = "~ 1 | Species",
#'                     phylo = sim$tree, data = sim$data, prop = 0.2,
#'                     runs = 2, nitt = 400, burnin = 100, thin = 2,
#'                     n_final = 5, skip_conv = TRUE, verbose = FALSE,
#'                     seed = 1)
#' }
#' @export
bace_validate <- function(x, ...) UseMethod("bace_validate")

#' @rdname bace_validate
#' @export
bace_validate.bace_complete <- function(x, prop = 0.1, reps = 1, scale = NULL,
                                        seed = NULL, verbose = TRUE,
                                        n_cores = 1L, max_attempts = 3,
                                        skip_conv = FALSE, ...) {
  fin <- x$final_results
  if (is.null(fin$settings) || is.null(fin$trees_used)) {
    stop("This result was made by a BACE version that did not store its ",
         "settings, so they cannot be reused. Call bace_validate() with the ",
         "data and formulas directly instead (see ?bace_validate).",
         call. = FALSE)
  }
  # Model settings come from the finished run; a named nitt = ... here would
  # be silently ignored, so refuse it loudly (mirrors bace()'s dots warning).
  dots <- list(...)
  if (length(dots) > 0L) {
    nm <- names(dots); if (is.null(nm)) nm <- rep("", length(dots))
    nm[nm == ""] <- "<unnamed>"
    warning("bace_validate() on a finished run takes its model settings ",
            "from that run; argument(s) ignored: ",
            paste(nm, collapse = ", "),
            ". To choose settings yourself, call bace_validate() with the ",
            "data and formulas directly.", call. = FALSE)
  }
  st    <- fin$settings
  trees <- fin$trees_used
  # The chained and final phases can have different MCMC settings. bace()
  # records both on the result; fin$settings holds only what the FINAL runs
  # used, so fall back to it for objects made before bace() stored its own.
  bs    <- if (is.list(x$settings)) x$settings else list()
  ch_nitt   <- if (!is.null(bs$nitt))   bs$nitt   else st$nitt
  ch_thin   <- if (!is.null(bs$thin))   bs$thin   else st$thin
  ch_burnin <- if (!is.null(bs$burnin)) bs$burnin else st$burnin
  .bace_validate_run(
    # the per-variable formulas the run actually used, as strings
    fixformula = lapply(fin$formulas,
                        function(f) paste(deparse(f), collapse = " ")),
    ran_phylo_form = st$ran_phylo_form,
    phylo = if (length(trees) == 1L) trees[[1L]] else trees,
    # the raw data as the run received it (original NAs intact)
    data = x$initial_results$data[["Initial_Data"]],
    prop = prop, reps = reps, scale = scale, seed = seed,
    nitt = ch_nitt, thin = ch_thin, burnin = ch_burnin,
    final_nitt = st$nitt, final_thin = st$thin, final_burnin = st$burnin,
    runs = length(x$initial_results$data) - 1L,  # chained iterations as run
    n_final = fin$n_final, species = st$species, verbose = verbose,
    max_attempts = max_attempts, skip_conv = skip_conv, n_cores = n_cores,
    nitt_cat_mult = st$nitt_cat_mult, ovr_categorical = st$ovr_categorical,
    priors = st$priors,
    sweeps = if (!is.null(fin$sweeps_requested)) fin$sweeps_requested else 3L,
    expected_types = fin$types)
}

#' @rdname bace_validate
#' @export
bace_validate.character <- function(x, ran_phylo_form, phylo, data,
                                    prop = 0.1, reps = 1, scale = NULL,
                                    seed = NULL, nitt = 6000, thin = 5,
                                    burnin = 1000, runs = 10, n_final = 50,
                                    species = FALSE, verbose = TRUE,
                                    max_attempts = 3, skip_conv = FALSE,
                                    n_cores = 1L, nitt_cat_mult = 1L,
                                    ovr_categorical = FALSE, priors = NULL,
                                    sweeps = 3L, final_nitt = NULL,
                                    final_thin = NULL, final_burnin = NULL,
                                    ...) {
  if (length(x) > 1L) {
    stop("Supply several formulas as a list(), not a character vector.",
         call. = FALSE)
  }
  .bace_validate_run(fixformula = x, ran_phylo_form = ran_phylo_form,
                     phylo = phylo, data = data, prop = prop, reps = reps,
                     scale = scale, seed = seed, nitt = nitt, thin = thin,
                     burnin = burnin, runs = runs, n_final = n_final,
                     species = species, verbose = verbose,
                     max_attempts = max_attempts, skip_conv = skip_conv,
                     n_cores = n_cores, nitt_cat_mult = nitt_cat_mult,
                     ovr_categorical = ovr_categorical, priors = priors,
                     sweeps = sweeps, final_nitt = final_nitt,
                     final_thin = final_thin, final_burnin = final_burnin)
}

#' @rdname bace_validate
#' @export
bace_validate.formula <- function(x, ran_phylo_form, phylo, data, ...) {
  bace_validate.character(paste(deparse(x), collapse = " "),
                          ran_phylo_form = ran_phylo_form, phylo = phylo,
                          data = data, ...)
}

#' @rdname bace_validate
#' @export
bace_validate.list <- function(x, ran_phylo_form, phylo, data,
                               prop = 0.1, reps = 1, scale = NULL,
                               seed = NULL, nitt = 6000, thin = 5,
                               burnin = 1000, runs = 10, n_final = 50,
                               species = FALSE, verbose = TRUE,
                               max_attempts = 3, skip_conv = FALSE,
                               n_cores = 1L, nitt_cat_mult = 1L,
                               ovr_categorical = FALSE, priors = NULL,
                               sweeps = 3L, final_nitt = NULL,
                               final_thin = NULL, final_burnin = NULL,
                               ...) {
  .bace_validate_run(fixformula = x, ran_phylo_form = ran_phylo_form,
                     phylo = phylo, data = data, prop = prop, reps = reps,
                     scale = scale, seed = seed, nitt = nitt, thin = thin,
                     burnin = burnin, runs = runs, n_final = n_final,
                     species = species, verbose = verbose,
                     max_attempts = max_attempts, skip_conv = skip_conv,
                     n_cores = n_cores, nitt_cat_mult = nitt_cat_mult,
                     ovr_categorical = ovr_categorical, priors = priors,
                     sweeps = sweeps, final_nitt = final_nitt,
                     final_thin = final_thin, final_burnin = final_burnin)
}

#' Per-repetition mask seeds, overflow-safe
#'
#' Derived through the RNG (not seed arithmetic): seed * 1000 overflows for
#' date-style seeds like 20261107. Any seed set.seed() accepts works here,
#' and the caller's random state is restored afterwards.
#' @keywords internal
#' @noRd
.derive_mask_seeds <- function(seed, reps) {
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  # Put the caller's random state back exactly as it was. If there was none
  # (a fresh session), remove the `.Random.seed` that set.seed() created here:
  # leaving it behind would silently continue this function's fixed stream into
  # the user's own random draws.
  on.exit({
    if (is.null(old_seed)) {
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  sample.int(.Machine$integer.max, reps)
}

#' The workhorse behind every bace_validate() method: mask, run, score.
#' @keywords internal
#' @noRd
.bace_validate_run <- function(fixformula, ran_phylo_form, phylo, data,
                          prop = 0.1, reps = 1, scale = NULL, seed = NULL,
                          nitt = 6000, thin = 5, burnin = 1000, runs = 10,
                          n_final = 50, species = FALSE, verbose = TRUE,
                          max_attempts = 3, skip_conv = FALSE, n_cores = 1L,
                          nitt_cat_mult = 1L, ovr_categorical = FALSE,
                          priors = NULL, sweeps = 3L,
                          final_nitt = NULL, final_thin = NULL,
                          final_burnin = NULL, expected_types = NULL) {

  if (!is.numeric(prop) || length(prop) != 1L || prop <= 0 || prop >= 1) {
    stop("`prop` must be a single number between 0 and 1.", call. = FALSE)
  }
  reps <- as.integer(reps)
  if (is.na(reps) || reps < 1L) stop("`reps` must be a whole number >= 1.",
                                     call. = FALSE)

  # Variables to mask: only the ones BACE imputes, i.e. the LEFT-hand side of
  # some formula. A complete covariate that only appears as a predictor has
  # no imputation formula -- masking it would make the pipeline error.
  if (is.list(fixformula)) {
    vars <- unique(unlist(lapply(fixformula,
                                 function(f) all.vars(stats::as.formula(f)[[2]]))))
  } else {
    # a single formula string: BACE builds one formula per named variable,
    # so every named variable is imputed and can be masked
    vars <- as.vector(unlist(.get_variables(fixformula)))
  }
  cluster <- .get_variables(ran_phylo_form, fix = FALSE)[["cluster"]]
  vars <- setdiff(vars, cluster)

  # scale sanity: a typo must not silently score a variable on the wrong scale
  if (!is.null(scale)) {
    bad_names <- setdiff(names(scale), vars)
    if (length(bad_names) > 0L) {
      warning("`scale` names not among the validated variables (ignored): ",
              paste(bad_names, collapse = ", "), call. = FALSE)
    }
    bad_vals <- setdiff(unique(unname(scale)), c("logit", "raw"))
    if (length(bad_vals) > 0L) {
      warning("`scale` values must be \"logit\" or \"raw\"; ignoring: ",
              paste(bad_vals, collapse = ", "), call. = FALSE)
    }
  }

  # One mask seed per repetition, fixed BEFORE any model runs so masks and
  # MCMC draws can never get tangled up. Each repetition's pipeline also gets
  # its own seed (the mask seed), so repetitions are independent and
  # reproducible. With seed = NULL the mask seeds are drawn from the session
  # RNG and recorded, so the masks remain reportable afterwards.
  mask_seeds <- .derive_mask_seeds(seed, reps)

  per_rep_tables <- vector("list", reps)
  per_rep_cells  <- vector("list", reps)

  for (r in seq_len(reps)) {
    if (verbose && reps > 1L) {
      cat("\n=== bace_validate: repetition", r, "of", reps, "===\n")
    }

    masked <- .mask_observed_cells(data, vars = vars,
                                   prop = prop, mask_seed = mask_seeds[r])
    if (nrow(masked$cells) == 0L) {
      stop("No cells could be masked (too few observed values?).",
           call. = FALSE)
    }

    fit <- bace(fixformula = fixformula, ran_phylo_form = ran_phylo_form,
                phylo = phylo, data = masked$data,
                nitt = nitt, thin = thin, burnin = burnin, runs = runs,
                n_final = n_final, species = species, verbose = verbose,
                max_attempts = max_attempts, skip_conv = skip_conv,
                n_cores = n_cores, nitt_cat_mult = nitt_cat_mult,
                ovr_categorical = ovr_categorical, priors = priors,
                sweeps = sweeps, seed = mask_seeds[r],
                final_nitt = final_nitt, final_thin = final_thin,
                final_burnin = final_burnin)

    # Masking can, in corner cases, flip a variable's detected type (e.g.
    # hiding an integer variable's only negative value makes it look like a
    # count). Say so rather than score a different family in silence.
    if (!is.null(expected_types)) {
      got <- fit$final_results$types
      for (v in intersect(names(expected_types), names(got))) {
        if (!identical(expected_types[[v]], got[[v]])) {
          warning("Masking changed the detected type of '", v, "' from ",
                  expected_types[[v]], " to ", got[[v]],
                  "; its scores reflect the re-detected type.",
                  call. = FALSE)
        }
      }
    }

    per_rep_tables[[r]] <- .score_validation(
      truth_data   = data,
      masked_cells = masked$cells,
      masked_data  = masked$data,
      datasets     = fit$imputed_datasets,
      prob_preds   = fit$prob_preds,
      types        = fit$final_results$types,
      scale        = scale)
    per_rep_tables[[r]]$rep <- r
    masked$cells$rep <- r
    per_rep_cells[[r]] <- masked$cells
  }

  table_all <- do.call(rbind, per_rep_tables)

  out <- list(
    table        = table_all,
    summary      = .summarise_validation(table_all),
    masked_cells = do.call(rbind, per_rep_cells),
    settings     = list(prop = prop, reps = reps, seed = seed,
                        mask_seeds = mask_seeds, scale = scale,
                        n_final = n_final, nitt = nitt, burnin = burnin,
                        thin = thin, sweeps = sweeps),
    n_final      = n_final
  )
  class(out) <- "bace_validation"
  if (verbose) print(out)
  invisible(out)
}

# ---------------------------------------------------------------------------
# Masking: hide a share of the observed cells
# ---------------------------------------------------------------------------

#' Hide a share of the observed cells, safely
#'
#' Picks cells to hide completely at random among the OBSERVED cells of each
#' variable, with three protections: a variable with fewer than 10 observed
#' cells is skipped (too few to score; a warning says so); at least 2
#' observed cells always remain (the error scale needs them); and for factors
#' no observed class is ever emptied (losing a class would change the model
#' itself, not just the data).
#'
#' @return list(data = the data with the chosen cells set to NA,
#'   cells = data.frame(row, variable, truth) recording what was hidden;
#'   `truth` is stored as character so mixed types fit one column).
#' @keywords internal
#' @noRd
.mask_observed_cells <- function(data, vars, prop, mask_seed) {
  # own RNG sandbox: restore the caller's random state afterwards
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  # Put the caller's random state back exactly as it was. If there was none
  # (a fresh session), remove the `.Random.seed` that set.seed() created here:
  # leaving it behind would silently continue this function's fixed stream into
  # the user's own random draws.
  on.exit({
    if (is.null(old_seed)) {
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  set.seed(mask_seed)

  out_data <- data
  cells <- list()

  for (v in vars) {
    obs <- which(!is.na(data[[v]]))
    if (length(obs) < 10L) {
      warning("bace_validate: '", v, "' has only ", length(obs),
              " observed cells (need at least 10) -- skipped.", call. = FALSE)
      next
    }
    if (is.factor(data[[v]])) {
      # Sample within each observed class, always leaving at least one
      # observed member per class.
      take <- unlist(lapply(
        split(obs, droplevels(data[[v]][obs])),
        function(ii) {
          k <- min(length(ii) - 1L, round(prop * length(ii)))
          if (k > 0L) sample(ii, k) else integer(0)
        }), use.names = FALSE)
    } else {
      # never mask everything: at least 2 observed cells must remain so the
      # variable's spread (the zRMSE denominator) exists
      k <- min(length(obs) - 2L, max(1L, round(prop * length(obs))))
      take <- if (k > 0L) sample(obs, k) else integer(0)
    }
    if (length(take) == 0L) {
      warning("bace_validate: '", v, "' -- no cells could be masked without ",
              "emptying a class; skipped.", call. = FALSE)
      next
    }
    out_data[take, v] <- NA
    cells[[v]] <- data.frame(row = sort(take), variable = v,
                             truth = as.character(data[[v]][sort(take)]),
                             stringsAsFactors = FALSE)
  }

  list(data = out_data,
       cells = if (length(cells)) do.call(rbind, c(cells, list(make.row.names = FALSE)))
               else data.frame(row = integer(0), variable = character(0),
                               truth = character(0)))
}

# ---------------------------------------------------------------------------
# Scoring: compare imputations to the hidden truth
# ---------------------------------------------------------------------------

#' Score imputations against hidden truth
#'
#' This is the arithmetic core of bace_validate(), kept separate so it can be
#' tested on hand-made inputs with known answers. One row per variable.
#'
#' @param truth_data the ORIGINAL data (nothing hidden)
#' @param masked_cells data.frame(row, variable, truth) from the masker
#' @param masked_data the data the pipeline actually saw (hidden cells = NA)
#' @param datasets list of M completed datasets
#' @param prob_preds per-run class-probability matrices (factors only)
#' @param types named list of variable types from the pipeline
#' @param scale optional named vector marking variables to score on the
#'   logit scale
#' @keywords internal
#' @noRd
.score_validation <- function(truth_data, masked_cells, masked_data,
                              datasets, prob_preds, types, scale = NULL) {
  M <- length(datasets)
  rows_out <- list()

  for (v in unique(masked_cells$variable)) {
    rows_v <- masked_cells$row[masked_cells$variable == v]
    type_v <- types[[v]]
    if (is.null(type_v)) next

    if (type_v %in% c("gaussian", "poisson")) {
      rows_out[[v]] <- .score_numeric_var(
        v, type_v, rows_v, truth_data, masked_data, datasets, M,
        use_logit = !is.null(scale) && identical(unname(scale[v]), "logit"))
    } else {
      rows_out[[v]] <- .score_class_var(
        v, type_v, rows_v, truth_data, masked_data, datasets, prob_preds, M)
    }
  }
  do.call(rbind, c(rows_out, list(make.row.names = FALSE)))
}

#' Score one continuous or count variable
#'
#' zRMSE exactly as in the BACE manuscript: the root-mean-square error of the
#' ACROSS-IMPUTATION MEAN at the hidden cells, divided by the standard
#' deviation of the variable's REMAINING observed cells. 0 = perfect
#' recovery; 1 = no better than filling in the observed mean. Coverage uses
#' the per-cell Rubin t-interval (same construction as
#' get_imputation_intervals()).
#' @keywords internal
#' @noRd
.score_numeric_var <- function(v, type_v, rows_v, truth_data, masked_data,
                               datasets, M, use_logit = FALSE) {
  # value of variable v at the hidden cells, in every completed dataset
  imp <- vapply(datasets, function(d) as.numeric(d[rows_v, v]),
                numeric(length(rows_v)))
  imp <- matrix(imp, nrow = length(rows_v))     # cells x imputations
  truth   <- as.numeric(truth_data[rows_v, v])
  obs_rest <- as.numeric(masked_data[[v]][!is.na(masked_data[[v]])])

  # optional logit scoring for proportion-like variables (ms convention:
  # clip to the range 0.0001 to 0.9999, then log-odds; the baseline SD
  # moves to the same scale)
  if (use_logit) {
    squash <- function(x) stats::qlogis(pmin(pmax(x, 1e-4), 1 - 1e-4))
    imp <- squash(imp); truth <- squash(truth); obs_rest <- squash(obs_rest)
  }

  ybar  <- rowMeans(imp)
  s     <- apply(imp, 1L, stats::sd)
  s_obs <- stats::sd(obs_rest)

  zrmse      <- sqrt(mean((ybar - truth)^2)) / s_obs
  zrmse_base <- sqrt(mean((mean(obs_rest) - truth)^2)) / s_obs
  r_pearson  <- if (length(truth) >= 3L && stats::sd(ybar) > 0 &&
                    stats::sd(truth) > 0) stats::cor(ybar, truth) else NA_real_

  # 95% Rubin t-interval per cell; count how often it caught the truth
  hw <- stats::qt(0.975, df = M - 1) * sqrt(1 + 1 / M) * s
  covered  <- mean(truth >= ybar - hw & truth <= ybar + hw)
  n_zero_w <- sum(s == 0)

  # Count variables: the two reference lines and a fair coverage expectation
  # (whole numbers make exactly 0.95 unattainable; always simulated from the
  # variable's own rates -- a mean-rate shortcut misreports mixtures).
  floor_abs <- floor_cal <- NA_real_
  cov_expect <- 0.95
  if (type_v == "poisson" && !use_logit) {
    floor_abs  <- sqrt(mean(ybar)) / s_obs          # counting noise alone
    floor_cal  <- sqrt((1 + 1 / M) * mean(s^2)) / s_obs  # calibrated level
    cov_expect <- .expected_coverage_counts(ybar, M)
  }

  data.frame(variable = v, type = type_v,
             scored_on = if (use_logit) "logit" else "raw",
             n_masked = length(rows_v),
             zrmse = zrmse, zrmse_baseline = zrmse_base, r = r_pearson,
             coverage = covered, coverage_expected = cov_expect,
             n_zero_width = n_zero_w,
             floor_counting = floor_abs, floor_calibrated = floor_cal,
             accuracy = NA_real_, accuracy_baseline = NA_real_,
             balanced_accuracy = NA_real_, brier = NA_real_,
             stringsAsFactors = FALSE)
}

#' Fractional-credit class vote (the manuscript's tie rule)
#'
#' Per cell: the classes with the most votes form the winning set; the cell
#' scores 1/|set| if the truth is in it, else 0 (a 2-way tie containing the
#' truth scores 0.5). This matches the BACE manuscript; breaking ties to the
#' alphabetically first class would allocate credit by label order.
#' @keywords internal
#' @noRd
.vote_credit <- function(votes_per_cell, truth) {
  vapply(seq_along(truth), function(i) {
    tab <- table(votes_per_cell[[i]])
    winners <- names(tab)[tab == max(tab)]
    if (truth[i] %in% winners) 1 / length(winners) else 0
  }, numeric(1))
}

#' Score one factor variable
#'
#' Accuracy: how often the most frequent imputed class was the hidden class,
#' with fractional credit for ties (see .vote_credit). Its baseline: always
#' guessing the most common class among the cells left observed (same tie
#' rule). Balanced accuracy: the same, averaged over the true classes, so
#' rare classes count as much as common ones. Brier: the mean squared error
#' of the imputed class probabilities (0 = perfect; 2 = confidently wrong).
#' @keywords internal
#' @noRd
.score_class_var <- function(v, type_v, rows_v, truth_data, masked_data,
                             datasets, prob_preds, M) {
  truth <- as.character(truth_data[rows_v, v])

  # all M imputed classes per cell ("votes")
  imp <- vapply(datasets, function(d) as.character(d[rows_v, v]),
                character(length(rows_v)))
  imp <- matrix(imp, nrow = length(rows_v))
  votes <- lapply(seq_len(nrow(imp)), function(i) imp[i, ])

  credit <- .vote_credit(votes, truth)
  acc <- mean(credit)
  per_class_acc <- vapply(unique(truth),
                          function(k) mean(credit[truth == k]), numeric(1))
  bal_acc <- mean(per_class_acc)

  # baseline: guess the most common class among the remaining observed cells,
  # with the same fractional tie rule
  obs_rest <- as.character(masked_data[[v]][!is.na(masked_data[[v]])])
  obs_tab  <- table(obs_rest)
  modal_set <- names(obs_tab)[obs_tab == max(obs_tab)]
  acc_base <- mean(vapply(truth, function(tt)
    if (tt %in% modal_set) 1 / length(modal_set) else 0, numeric(1)))

  # Brier score from the class probabilities, averaged over the M runs.
  # The probability matrices cover ALL missing rows of v (original NAs plus
  # the hidden cells), in the row order the data has; pick out the hidden
  # cells' rows by position.
  brier <- NA_real_
  pm_list <- lapply(prob_preds, function(pp) pp[[v]])
  pm_list <- pm_list[!vapply(pm_list, is.null, logical(1))]
  if (length(pm_list) > 0L) {
    na_rows <- which(is.na(masked_data[[v]]))
    pos <- match(rows_v, na_rows)
    if (!anyNA(pos) && all(dim(pm_list[[1]])[1] == length(na_rows))) {
      p_avg <- Reduce(`+`, pm_list) / length(pm_list)
      p_hid <- p_avg[pos, , drop = FALSE]
      classes <- colnames(p_hid)
      truth_ind <- outer(truth, classes, `==`) * 1
      brier <- mean(rowSums((p_hid - truth_ind)^2))
    }
  }

  data.frame(variable = v, type = type_v, scored_on = "class",
             n_masked = length(rows_v),
             zrmse = NA_real_, zrmse_baseline = NA_real_, r = NA_real_,
             coverage = NA_real_, coverage_expected = NA_real_,
             n_zero_width = NA_integer_,
             floor_counting = NA_real_, floor_calibrated = NA_real_,
             accuracy = acc, accuracy_baseline = acc_base,
             balanced_accuracy = bal_acc, brier = brier,
             stringsAsFactors = FALSE)
}

#' Fair coverage expectation for counts
#'
#' Counts are whole numbers, so a 95 percent interval cannot always hit
#' exactly 95 percent. This simulates the attainable level for the variable's
#' OWN rates: pretend the imputer is exactly right (truth and imputations
#' drawn from the same Poisson rates) and measure the coverage the t-interval
#' then achieves. Run under a fixed internal seed (it is a reference value --
#' two identical validations must print the same expectation), with the
#' caller's RNG state restored. A genuinely calibrated imputer with parameter
#' uncertainty typically lands within about 0.01 of this value.
#' @keywords internal
#' @noRd
.expected_coverage_counts <- function(rates, M, n_total = 20000L) {
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  # Put the caller's random state back exactly as it was. If there was none
  # (a fresh session), remove the `.Random.seed` that set.seed() created here:
  # leaving it behind would silently continue this function's fixed stream into
  # the user's own random draws.
  on.exit({
    if (is.null(old_seed)) {
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  set.seed(20261107L)

  rates <- pmax(rates, 1e-6)
  n_sim <- max(400L, ceiling(n_total / length(rates)))
  t_mult <- stats::qt(0.975, df = M - 1)
  covered <- vapply(rates, function(r) {
    truth <- stats::rpois(n_sim, r)
    draws <- matrix(stats::rpois(n_sim * M, r), nrow = n_sim)
    m  <- rowMeans(draws)
    hw <- t_mult * sqrt(1 + 1 / M) * apply(draws, 1L, stats::sd)
    mean(truth >= m - hw & truth <= m + hw)
  }, numeric(1))
  mean(covered)
}

#' Average the per-repetition scores
#'
#' One row per variable: the mean of each score over repetitions and, when
#' there is more than one repetition, the standard error of that mean. NAs
#' from degenerate repetitions are dropped per metric (the count used is in
#' reps; a metric absent in every repetition stays NA).
#' @keywords internal
#' @noRd
.summarise_validation <- function(table_all) {
  num_cols <- c("zrmse", "zrmse_baseline", "r", "coverage",
                "coverage_expected", "floor_counting", "floor_calibrated",
                "accuracy", "accuracy_baseline", "balanced_accuracy", "brier")
  out <- lapply(split(table_all, table_all$variable), function(g) {
    means <- vapply(num_cols, function(cc) {
      x <- g[[cc]][!is.na(g[[cc]])]
      if (length(x)) mean(x) else NA_real_
    }, numeric(1))
    ses <- vapply(num_cols, function(cc) {
      x <- g[[cc]][!is.na(g[[cc]])]
      if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
    }, numeric(1))
    names(ses) <- paste0(num_cols, "_se")
    cbind(data.frame(variable = g$variable[1], type = g$type[1],
                     scored_on = g$scored_on[1],
                     n_masked = round(mean(g$n_masked)), reps = nrow(g),
                     stringsAsFactors = FALSE),
          as.data.frame(as.list(means)), as.data.frame(as.list(ses)))
  })
  do.call(rbind, c(out, list(make.row.names = FALSE)))
}

# ---------------------------------------------------------------------------
# Printing: a plain-language report
# ---------------------------------------------------------------------------

#' @export
print.bace_validation <- function(x, ...) {
  s <- x$summary
  reps <- x$settings$reps
  cat("\n=== BACE validation: hide known cells, impute, compare ===\n\n")
  cat("Hidden per variable: ~", round(100 * x$settings$prop),
      "% of observed cells",
      if (reps > 1) paste0(" x ", reps, " repetitions"),
      "; ", x$n_final, " imputations per repetition.\n", sep = "")
  if (reps > 1L) {
    cat("Values are mean (+/- 1 standard error of the mean) over ",
        reps, " repetitions.\n", sep = "")
  }
  cat("\n")

  fmt <- function(val, se = NA) {
    if (is.na(val)) return("--")
    if (!is.na(se)) sprintf("%.2f (+/- %.2f)", val, se) else sprintf("%.2f", val)
  }

  for (i in seq_len(nrow(s))) {
    g <- s[i, ]
    cat(g$variable, " (", g$type,
        if (g$scored_on == "logit") ", scored on log-odds scale" else "",
        ", ", g$n_masked, " hidden cells):\n", sep = "")
    if (g$scored_on != "class") {
      if (is.na(g$zrmse)) {
        cat("  error (zRMSE): -- (not computable: too few observed cells",
            "remained)\n")
      } else {
        cat("  error (zRMSE): ", fmt(g$zrmse, g$zrmse_se),
            "   [0 = perfect; mean-filling scores ",
            fmt(g$zrmse_baseline), "]\n", sep = "")
      }
      if (!is.na(g$floor_counting)) {
        cat("  count references: counting-noise floor ",
            fmt(g$floor_counting), "; a calibrated imputer sits near ",
            fmt(g$floor_calibrated),
            " -- a zRMSE near these lines is a GOOD result\n", sep = "")
      }
      if (!is.na(g$r)) cat("  correlation with truth: ", fmt(g$r, g$r_se),
                           "\n", sep = "")
      if (!is.na(g$coverage)) {
        cat("  95% interval coverage: ", fmt(g$coverage, g$coverage_se),
            "   [expect about ", fmt(g$coverage_expected),
            if (!is.na(g$floor_counting))
              "; calibrated imputers land within ~0.01 of this" else "",
            "]\n", sep = "")
      }
      nz <- x$table$n_zero_width[x$table$variable == g$variable]
      if (any(nz > 0, na.rm = TRUE)) {
        cat("  note: ", sum(nz, na.rm = TRUE), " cell-instance(s) across ",
            reps, " repetition(s) had all imputations identical ",
            "(zero-width intervals)\n", sep = "")
      }
    } else {
      cat("  accuracy: ", fmt(g$accuracy, g$accuracy_se),
          "   [guessing the commonest observed class scores ",
          fmt(g$accuracy_baseline), "; ties get fractional credit]\n",
          sep = "")
      cat("  balanced accuracy (rare classes count equally): ",
          fmt(g$balanced_accuracy, g$balanced_accuracy_se), "\n", sep = "")
      if (!is.na(g$brier)) {
        cat("  probability error (Brier; 0 = perfect, 2 = confidently",
            "wrong): ", fmt(g$brier, g$brier_se), "\n", sep = "")
      }
    }
    cat("\n")
  }

  if (reps == 1L) {
    cat("Single repetition: treat the numbers as estimates; rerun with",
        "reps >= 2\nfor an honest between-repetition uncertainty.\n")
  }
  cat("Keep in mind: this checks prediction where data were OBSERVED.",
      "It cannot\ndetect bias from data that are missing not-at-random",
      "(e.g. unmeasured rare\nspecies).\n")
  invisible(x)
}
