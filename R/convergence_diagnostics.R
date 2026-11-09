#' @title assess_convergence
#' @description Check a BACE convergence chain for signs that the chained
#'   equations have NOT settled. The check watches, for every variable with
#'   missing data, two series across the chained runs: the mean and the
#'   standard deviation of that variable's imputed values (the two quantities
#'   recommended for monitoring by van Buuren & Groothuis-Oudshoorn 2011 and
#'   van Buuren 2018, ch. 6.5.2). The first half of the chain is treated as
#'   warm-up; all movement is measured against the variable's own imputation
#'   spread, never against its level, so the verdict does not depend on where
#'   or on what scale a variable lives.
#'
#'   Three named pathologies are flagged:
#'   \itemize{
#'     \item \strong{drift}: the late chain is still travelling — the mean
#'       moved by more than \code{drift_threshold} spread units across the
#'       retained window, with a supporting trend. The one pathology that
#'       more chained runs can cure.
#'     \item \strong{oscillation}: the mean flip-flops between runs, loudly —
#'       strongly negative lag-1 autocorrelation AND steps well above the
#'       Monte-Carlo noise floor. Usually two variables feeding back on each
#'       other; more runs will not fix it (van Buuren 2018, ch. 6.5.2).
#'     \item \strong{stuck}: the imputation spread is collapsing towards zero
#'       run after run — the textbook pathology in van Buuren &
#'       Groothuis-Oudshoorn (2011, p. 37), where the sampler locks up. The
#'       pre-1.0.1 check called this "converged", which was exactly backwards.
#'   }
#'
#'   \strong{What a pass means.} "Converged" here means \emph{no pathology
#'   detected} — a one-directional diagnostic, not a certificate. There is no
#'   agreed numerical test of chained-equations convergence ("There is no
#'   clear-cut method for determining when the MICE algorithm has converged",
#'   van Buuren 2018, ch. 6.5.2), and short chains cannot support one: with
#'   the default \code{runs = 10} the check sees ten points per variable.
#'   Empirically, valid inference tends to arrive after 5--10 chained
#'   iterations, well before formal diagnostics look settled (Oberman, van
#'   Buuren & Vink 2021, arXiv:2110.11951). Always inspect the trace plots
#'   (\code{plot(<result>)}) rather than relying on the flag alone.
#'
#'   Chained samplers also alternate a little by nature (each run conditions
#'   on the previous values of the other variables), and a proportion watched
#'   over few imputed cells moves in quantised jumps. Mild versions of both
#'   are therefore reported as notes for the eye, not failures.
#' @param bace_object An object of class 'bace' returned from bace_imp
#' @param method "stability" (the check described here). "summary" is accepted
#'   as a deprecated alias and warns. The pre-1.0.1 methods ("summary.acf",
#'   "summary.percentage", "summary.trend", "summary.geweke", "energy",
#'   "wasserstein", "all") are defunct: all of them measured a statistic
#'   against its own level or mean, so their verdicts tracked the variable's
#'   scale rather than its stability, failing settled chains and passing
#'   drifting ones. See NEWS for the full account.
#' @param variables Character vector of variable names to assess. If NULL
#'   (default), all variables with missing data are assessed. A warning is
#'   issued for variables without missing data, which are ignored.
#' @param use_all_data Logical. If TRUE, monitor whole columns rather than
#'   only the imputed cells. Default FALSE (imputed cells only, the
#'   recommended practice — observed cells never change, so including them
#'   only dilutes the signal).
#' @param min_iterations Minimum number of chained runs before any assessment
#'   is attempted. Default is 3.
#' @param drift_threshold Net late-chain movement of the mean, in units of
#'   the variable's imputation spread, above which (with a supporting trend)
#'   the variable is flagged as drifting. Default 0.75. Calibrated so that
#'   settled chains pass >= 99 per cent of the time and a steady drift of half a
#'   spread per run is caught >= 99 per cent of the time at \code{runs = 10}; at
#'   \code{runs = 5} the drift catch rate is ~75 per cent, which is as much as five
#'   points can support.
#' @param stuck_ratio If the late-chain imputation spread falls below this
#'   fraction of the early-chain spread, the variable is flagged as stuck.
#'   Default 0.25.
#' @return A list of class \code{bace_convergence}:
#'   \itemize{
#'     \item converged: TRUE when no pathology was detected in any variable
#'       (NA if nothing could be assessed)
#'     \item verdicts: per-variable, one of "ok", "drift", "oscillation",
#'       "stuck", or "unassessable" (fewer than 3 finite points)
#'     \item notes: per-variable advisory notes ("" when none)
#'     \item summary_stats: data frame of the mean series (iteration x
#'       variable) that the check watched
#'     \item sd_stats: data frame of the SD series
#'     \item method_results$summary: per-variable detail (verdict, the
#'       statistics behind it, and the thresholds used)
#'   }
#' @references
#'   van Buuren S & Groothuis-Oudshoorn K (2011). mice: Multivariate
#'   Imputation by Chained Equations in R. Journal of Statistical Software,
#'   45(3). \cr
#'   van Buuren S (2018). Flexible Imputation of Missing Data, 2nd ed.,
#'   ch. 6.5.2. CRC Press. \cr
#'   Oberman HI, van Buuren S & Vink G (2021). Missing the Point:
#'   Non-Convergence in Iterative Imputation Algorithms. arXiv:2110.11951.
#' @importFrom stats cor cor.test var sd median approx
#' @examples
#' \donttest{
#' # Simulate a small phylogenetic dataset with missing values, run a short
#' # bace_imp() convergence chain, then check it for pathologies
#' set.seed(1)
#' sim <- sim_bace(response_type = "gaussian", predictor_types = "gaussian",
#'                 phylo_signal = c(0.6, 0.6), missingness = c(0.2, 0),
#'                 n_cases = 40, n_species = 40)
#' fit <- bace_imp(fixformula = "y ~ x1", ran_phylo_form = "~ 1 | Species",
#'                 phylo = sim$tree, data = sim$data,
#'                 runs = 3, nitt = 600, burnin = 100, thin = 2, verbose = FALSE)
#' conv <- assess_convergence(fit)
#' conv$converged
#' }
#' @export
assess_convergence <- function(bace_object,
                               method = c("stability", "summary"),
                               variables = NULL,
                               use_all_data = FALSE,
                               min_iterations = 3,
                               drift_threshold = 0.75,
                               stuck_ratio = 0.25) {

  # Check inputs
  if (!inherits(bace_object, "bace")) {
    stop("Input must be an object of class 'bace' from bace_imp function")
  }

  # The pre-1.0.1 methods are defunct, with an explanation rather than a
  # match.arg error: their verdicts measured scale, not stability.
  defunct <- c("summary.acf", "summary.percentage", "summary.trend",
               "summary.geweke", "energy", "wasserstein", "all")
  if (length(method) == 1 && method %in% defunct) {
    stop("method = \"", method, "\" is defunct. The pre-1.0.1 convergence ",
         "criteria compared each statistic to its own level or mean, so ",
         "their verdicts tracked the variable's scale rather than its ",
         "stability (failing settled chains and passing drifting ones). ",
         "Use the default method = \"stability\" and inspect plot() of the ",
         "result. See NEWS for details.")
  }
  method <- match.arg(method)
  if (method == "summary") {
    warning("method = \"summary\" is deprecated; using method = \"stability\" ",
            "(the repaired check). See NEWS.")
    method <- "stability"
  }

  if (!is.numeric(min_iterations) || length(min_iterations) != 1 ||
      !is.finite(min_iterations) || min_iterations < 3) {
    stop("min_iterations must be a single number >= 3: fewer than three ",
         "chained runs cannot support any stability statement.")
  }

  # Extract imputed datasets (element 1 is the un-imputed starting frame)
  imputed_data <- bace_object$data[-1]
  n_iter <- length(imputed_data)

  if (n_iter < min_iterations) {
    warning(paste("Only", n_iter, "iterations available. Minimum of",
                  min_iterations, "recommended for convergence assessment"))
    return(list(converged = NA,
                message = "Insufficient iterations for convergence assessment"))
  }

  # Get variable types and missing data info
  types <- bace_object$types
  miss_dat <- bace_object$miss_dat

  # Filter variables if specified
  vars_with_missing <- unique(miss_dat$colname)

  if (!is.null(variables)) {
    vars_without_missing <- setdiff(variables, vars_with_missing)

    if (length(vars_without_missing) > 0) {
      warning(paste("The following variables have no missing data and will be ignored:",
                    paste(vars_without_missing, collapse = ", ")))
    }

    variables <- intersect(variables, vars_with_missing)

    if (length(variables) == 0) {
      stop("None of the specified variables have missing data")
    }
  }

  summary_result <- .assess_summary_convergence(
    imputed_data, types, miss_dat,
    variables = variables,
    use_all_data = use_all_data,
    drift_threshold = drift_threshold,
    stuck_ratio = stuck_ratio)

  results <- list(
    converged = summary_result$converged,
    verdicts = summary_result$verdicts,
    notes = summary_result$notes,
    method_results = list(summary = summary_result),
    summary_stats = summary_result$stats_df,
    sd_stats = summary_result$sd_df,
    diagnostics = list(n_iter = n_iter,
                       warmup = summary_result$warmup,
                       thresholds = summary_result$thresholds)
  )

  results$types <- bace_object$types[setdiff(names(results$summary_stats), "iteration")]

  class(results) <- c("bace_convergence", "list")
  return(results)
}


#' @title .conv_assess_series
#' @description The stability check for ONE monitored series: the shared core
#'   used by the convergence assessment and by the diagnostic plots, so the
#'   two can never disagree. Scale-free and location-free throughout: all
#'   movement is measured in units of the imputation spread, never the level.
#' @param m Numeric vector: the monitored statistic at each chained run — the
#'   mean of the imputed values (\code{kind = "numeric"}) or the proportion of
#'   imputed cells in one level (\code{kind = "proportion"}).
#' @param s Numeric vector: for \code{kind = "numeric"}, the SD of those
#'   imputed values at each run — an independent measurement of how uncertain
#'   the imputations are, so it can be watched for collapse or explosion. For
#'   \code{kind = "proportion"} it is ignored: sqrt(p(1-p)) is an arithmetic
#'   function of \code{m}, not a second measurement, so watching it would
#'   re-read the mean and mislabel an extreme proportion as a collapse.
#' @param n_cells How many cells the statistic at each run is computed over.
#' @param kind "numeric" or "proportion". Sets the yardstick: a numeric series
#'   is standardised by its own imputation spread, a proportion is already on
#'   a 0--1 scale and is measured in proportion units.
#' @param drift_threshold,stuck_ratio See \code{assess_convergence}.
#' @param osc_acf,osc_amp,osc_floor,wander,trend_alpha Fixed internals.
#'   Oscillation needs lag-1 autocorrelation < \code{osc_acf} AND a mean step
#'   above BOTH \code{osc_amp} noise floors and the absolute floor
#'   \code{osc_floor}. The absolute floor matters: the Monte-Carlo noise floor
#'   sqrt(2/n_cells) shrinks with the number of imputed cells, but real
#'   between-run movement is refit noise that does not, so without it the
#'   amplitude gate stops screening on large datasets and a lag-1 correlation
#'   over a handful of pairs decides the verdict on its own (false
#'   oscillation ~25% at n_cells = 1000). \code{wander} catches a chain that
#'   keeps taking large steps without going anywhere. Calibrated against
#'   constructed chains and the real harvested ones; see
#'   \code{.agents/investigation-convergence-2026-11.md}.
#' @return list(verdict, note, stats)
#' @noRd
.conv_assess_series <- function(m, s, n_cells,
                                kind = c("numeric", "proportion"),
                                drift_threshold = 0.75,
                                stuck_ratio = 0.25,
                                osc_acf = -0.60,
                                osc_amp = 3.0,
                                osc_floor = 0.60,
                                wander = 1.00,
                                trend_alpha = 0.10) {

  kind <- match.arg(kind)
  if (missing(s) || is.null(s)) s <- rep(NA_real_, length(m))
  if (!is.finite(n_cells) || n_cells < 1) n_cells <- 1

  ## A series that cannot be assessed must say so, never pass silently.
  ## Non-finite runs (a run whose imputations are all NA) drop out pairwise.
  finite <- is.finite(m) & (is.finite(s) | kind == "proportion")
  base_stats <- c(net = NA_real_, step = NA_real_, acf1 = NA_real_,
                  sd_ratio = NA_real_, n_kept = NA_real_, ref = NA_real_)
  if (sum(finite) < 3) {
    base_stats["n_kept"] <- sum(finite)
    return(list(verdict = "unassessable",
                note = "fewer than 3 finite points in the monitored series",
                stats = base_stats))
  }
  m <- m[finite]
  s <- s[finite]

  runs <- length(m)
  keep <- max(1, floor(runs / 2)):runs        # discard the warm-up half
  mk   <- m[keep]
  sk   <- s[keep]
  note <- character(0)

  add_note <- function(x) note <<- c(note, x)

  out <- function(verdict, net = NA_real_, step = NA_real_, acf1 = NA_real_,
                  sd_ratio = NA_real_, ref = NA_real_) {
    list(verdict = verdict,
         note = if (length(note)) paste(note, collapse = "; ") else "",
         stats = c(net = net, step = step, acf1 = acf1, sd_ratio = sd_ratio,
                   n_kept = length(mk), ref = ref))
  }

  ## 0. constant chain: identical imputations every run. Legitimate certainty
  ##    (e.g. a binary trait the model is sure about), not a pathology — but
  ##    worth a note, since a frozen chain also looks like this. The tolerance
  ##    is RELATIVE so the branch cannot swallow a drift that happens to live
  ##    at a tiny absolute scale.
  ## The tolerance is RELATIVE to the series' own scale (with a true-zero
  ## special case), so this branch cannot swallow a drift that happens to
  ## live at a tiny absolute scale — an absolute 1e-12 floor did exactly
  ## that, calling a clean 0 -> 1e-13 ramp "identical imputations".
  m_scale <- max(abs(m), na.rm = TRUE)
  tol     <- if (m_scale > 0) 1e-12 * m_scale else 0
  s_const <- if (kind == "proportion") TRUE else all(s <= tol, na.rm = TRUE)
  if (stats::sd(m) <= tol && s_const) {
    add_note("identical imputations in every run; inspect if unexpected")
    return(out("ok", net = 0, step = 0, sd_ratio = 1))
  }

  ## 1. the imputation spread, watched in BOTH directions — but only when it
  ##    is a real measurement. For a proportion, sqrt(p(1-p)) is a function of
  ##    the mean, so reading it as a spread turns "the model became confident"
  ##    into a false collapse (and a drifting proportion into a false stuck).
  sd_ratio <- NA_real_
  if (kind == "numeric") {
    first_s <- stats::median(s[seq_len(ceiling(runs / 2))], na.rm = TRUE)
    late_s  <- mean(utils::tail(sk, 2), na.rm = TRUE)
    if (is.finite(first_s) && first_s > 0 && is.finite(late_s)) {
      sd_ratio <- late_s / first_s
      if (sd_ratio < stuck_ratio)
        return(out("stuck", sd_ratio = sd_ratio))
      if (sd_ratio > 1 / stuck_ratio)
        return(out("diverging", sd_ratio = sd_ratio))
    }
  }

  ## 2. the yardsticks.
  ##    numeric    -> standardise by the variable's own late imputation spread
  ##    proportion -> already on 0-1; measure in proportion units
  ##    ref is the Monte-Carlo noise floor of the statistic, with an absolute
  ##    floor (see @param osc_floor): it shrinks with n_cells, real refit
  ##    movement does not.
  if (kind == "numeric") {
    spread <- mean(sk, na.rm = TRUE)
    if (!is.finite(spread) || spread <= 0) spread <- stats::sd(mk)
    if (!is.finite(spread) || spread <= 0) spread <- .Machine$double.eps
    ref <- sqrt(2 / max(n_cells, 1))
    amp_gate <- max(osc_amp * ref, osc_floor)
  } else {
    ## A proportion is already on 0-1, so it carries its own scale: movement
    ## is in proportion units and the largest possible step is 1 (every cell
    ## changing level). The gates are therefore smaller AND capped below 1,
    ## or they could never fire when few cells make the sampling noise of the
    ## proportion itself large.
    spread <- 1
    pbar <- mean(mk, na.rm = TRUE)
    ref  <- max(sqrt(2 * pbar * (1 - pbar) / max(n_cells, 1)),
                1 / max(n_cells, 1))
    ## A proportion over FEW cells is intrinsically jumpy (sampling SD
    ## sqrt(p(1-p)/n)), so the drift bar scales with that noise: tight when
    ## there is enough data to see a drift, honest when there is not.
    drift_threshold <- max(min(drift_threshold, 0.20), 2 * ref)
    amp_gate <- min(max(2 * ref, 0.15), 0.80)
  }
  steps     <- abs(diff(mk)) / spread
  step_mean <- mean(steps, na.rm = TRUE)

  ## 3. flip-flop: strong negative lag-1 autocorrelation with LOUD amplitude.
  ##    Chained equations alternate a little by nature, and a proportion over
  ##    few cells moves in quantised jumps, so the amplitude must clear an
  ##    absolute floor as well as the noise floor.
  acf1 <- NA_real_
  if (length(mk) >= 5 && stats::sd(mk) > 0) {
    acf1 <- suppressWarnings(stats::cor(mk[-length(mk)], mk[-1]))
    if (is.finite(acf1) && acf1 < osc_acf && step_mean > amp_gate)
      return(out("oscillation", step = step_mean, acf1 = acf1,
                 sd_ratio = sd_ratio, ref = ref))
    if (is.finite(acf1) && acf1 < -0.7 && step_mean > ref)
      add_note("alternating structure; inspect the trace")
  } else if (length(mk) < 5) {
    add_note("chain too short to test for flip-flopping (needs runs >= 7)")
  }

  ## 4. still moving? Either a net journey (with trend support, or so large
  ##    that no test is needed) or sustained large steps going nowhere.
  net <- NA_real_
  k   <- max(1, floor(length(mk) / 3))
  net <- abs(mean(utils::tail(mk, k)) - mean(utils::head(mk, k))) / spread
  trend_p <- if (length(mk) >= 4) {
    suppressWarnings(stats::cor.test(seq_along(mk), mk,
                                     method = "kendall")$p.value)
  } else {
    NA_real_
  }
  big_net <- net > 3 * drift_threshold          # too large to need a test
  trended <- net > drift_threshold &&
             is.finite(trend_p) && trend_p < trend_alpha
  if (big_net || trended || step_mean > wander)
    return(out("drift", net = net, step = step_mean, acf1 = acf1,
               sd_ratio = sd_ratio, ref = ref))
  if (length(mk) < 4)
    add_note("chain too short to test for drift (needs runs >= 5)")

  out("ok", net = net, step = step_mean, acf1 = acf1,
      sd_ratio = sd_ratio, ref = ref)
}


#' @title .assess_summary_convergence
#' @description Internal: build the monitored series (mean and SD of the
#'   imputed values per chained run, per variable) and run the stability
#'   check on each. See \code{assess_convergence} for semantics.
#'
#'   Categorical and ordinal variables are monitored on EVERY level, not just
#'   the first, and take their variable's verdict from the worst level. The
#'   first level alone is blind to movement among the others: a chain flipping
#'   all its imputed cells between levels 2 and 3 leaves level 1's proportion
#'   flat, which reads as a perfectly settled chain.
#' @param imputed_data List of imputed datasets across iterations
#' @param types Named list of variable types
#' @param miss_dat Data frame with missing data information
#' @param variables Character vector of variables to check, or NULL for all
#' @param use_all_data Logical. Use all data (TRUE) or only imputed values (FALSE)
#' @param drift_threshold,stuck_ratio See \code{assess_convergence}
#' @return List with convergence assessment
#' @noRd
.assess_summary_convergence <- function(imputed_data, types, miss_dat,
                                        variables = NULL,
                                        use_all_data = FALSE,
                                        drift_threshold = 0.75,
                                        stuck_ratio = 0.25) {

  n_iter <- length(imputed_data)

  # Get unique variables with missing data
  vars_to_check <- unique(miss_dat$colname)

  # Filter to specified variables if provided
  if (!is.null(variables)) {
    vars_to_check <- intersect(vars_to_check, variables)
  }

  ## Nothing to assess is a documented outcome, not a crash.
  if (length(vars_to_check) == 0) {
    empty <- data.frame(iteration = seq_len(n_iter))
    return(list(converged = NA, var_converged = logical(0),
                verdicts = character(0), notes = character(0),
                assessments = list(), stats_df = empty, sd_df = empty,
                warmup = max(1, floor(n_iter / 2)) - 1,
                thresholds = c(drift_threshold = drift_threshold,
                               stuck_ratio = stuck_ratio,
                               osc_acf = -0.60, osc_amp = 3.0,
                               osc_floor = 0.50, wander = 1.0,
                               trend_alpha = 0.10),
                types = list()))
  }

  mean_list <- list()   # one entry per variable: the series shown to the user
  sd_list   <- list()
  series    <- list()   # one entry per variable: list of series to ASSESS

  for (var in vars_to_check) {
    var_type <- types[[var]]
    if (is.null(var_type) ||
        !var_type %in% c("gaussian", "poisson", "categorical", "threshold")) {
      # The pre-1.0.1 code fell through silently here, leaving an all-zero
      # series that auto-passed every criterion. Unknown must be loud.
      stop("assess_convergence: variable '", var, "' has unrecognised type '",
           if (is.null(var_type)) "NULL" else var_type,
           "'; cannot monitor it.")
    }

    miss_idx <- miss_dat[miss_dat$colname == var, "row"]

    vals_by_run <- lapply(seq_len(n_iter), function(i) {
      if (use_all_data) imputed_data[[i]][[var]]
      else imputed_data[[i]][[var]][miss_idx]
    })
    n_cells <- length(vals_by_run[[1]])

    if (var_type %in% c("gaussian", "poisson")) {
      ## mean and SD of the imputed values: two independent measurements.
      mu <- vapply(vals_by_run, function(v) mean(v, na.rm = TRUE), numeric(1))
      sg <- vapply(vals_by_run, function(v) {
        z <- stats::sd(v, na.rm = TRUE)       # na.rm: one NA cell must not
        if (is.na(z)) 0 else z                # masquerade as a collapsed spread
      }, numeric(1))
      mean_list[[var]] <- mu
      sd_list[[var]]   <- sg
      series[[var]] <- list(list(label = var, kind = "numeric",
                                 m = mu, s = sg, n_cells = n_cells))
    } else {
      ## one proportion series per level; the variable takes the worst verdict.
      lev <- levels(vals_by_run[[1]])
      if (is.null(lev)) lev <- sort(unique(unlist(lapply(vals_by_run,
                                    function(v) as.character(v)))))
      props <- lapply(lev, function(L)
        vapply(vals_by_run,
               function(v) mean(as.character(v) == L, na.rm = TRUE),
               numeric(1)))
      names(props) <- lev
      ## shown to the user: the first level, as before (and its binomial sd)
      mean_list[[var]] <- props[[1]]
      sd_list[[var]]   <- sqrt(props[[1]] * (1 - props[[1]]))
      series[[var]] <- lapply(lev, function(L)
        list(label = paste0(var, ":", L), kind = "proportion",
             m = props[[L]], s = NULL, n_cells = n_cells))
    }
  }

  # Data frames of the monitored series (kept for plots and for the user)
  stats_df <- data.frame(iteration = seq_len(n_iter),
                         do.call(cbind, mean_list), check.names = FALSE)
  names(stats_df)[-1] <- names(mean_list)
  sd_df <- data.frame(iteration = seq_len(n_iter),
                      do.call(cbind, sd_list), check.names = FALSE)
  names(sd_df)[-1] <- names(sd_list)

  ## Assess every series; a variable's verdict is its worst level's.
  rank_of <- c(ok = 0, unassessable = 1, drift = 2, oscillation = 3,
               stuck = 4, diverging = 4)
  assessments <- lapply(vars_to_check, function(var) {
    parts <- lapply(series[[var]], function(z)
      .conv_assess_series(z$m, z$s, z$n_cells, kind = z$kind,
                          drift_threshold = drift_threshold,
                          stuck_ratio = stuck_ratio))
    worst <- which.max(vapply(parts, function(p)
      unname(rank_of[p$verdict]), numeric(1)))
    a <- parts[[worst]]
    ## name the offending level when it is not the first
    if (length(parts) > 1 && a$verdict != "ok") {
      lbl <- series[[var]][[worst]]$label
      a$note <- paste(c(a$note[nzchar(a$note)],
                        paste0("level ", sub("^[^:]*:", "", lbl))),
                      collapse = "; ")
    }
    if (length(parts) > 1 && a$verdict == "ok") {
      nts <- unique(unlist(lapply(parts, function(p) p$note)))
      nts <- nts[nzchar(nts)]
      a$note <- if (length(nts)) paste(nts, collapse = "; ") else ""
    }
    a
  })
  names(assessments) <- vars_to_check

  verdicts <- vapply(assessments, function(a) a$verdict, character(1))
  notes    <- vapply(assessments, function(a) a$note, character(1))
  var_converged <- verdicts == "ok"

  # One pathology anywhere is a failure, and it is named. (The pre-1.0.1 rule,
  # mean(var_converged) > 0.7, demanded 100% agreement for 1-3 variables but
  # only 71% for 7 — an accident of arithmetic, not a design.) A variable that
  # could not be assessed is NOT a pathology: if nothing could be assessed the
  # verdict is NA, as documented.
  overall_converged <- if (length(verdicts) == 0 ||
                           all(verdicts == "unassessable")) NA
                       else all(verdicts[verdicts != "unassessable"] == "ok")

  return(list(
    converged = overall_converged,
    var_converged = var_converged,
    verdicts = verdicts,
    notes = notes,
    assessments = assessments,
    stats_df = stats_df,
    sd_df = sd_df,
    warmup = max(1, floor(n_iter / 2)) - 1,
    thresholds = c(drift_threshold = drift_threshold,
                   stuck_ratio = stuck_ratio,
                   osc_acf = -0.60, osc_amp = 3.0, osc_floor = 0.50,
                   wander = 1.0, trend_alpha = 0.10),
    types = types[vars_to_check]
  ))
}


#' @title print.bace_convergence
#' @description Print method for bace_convergence objects
#' @param x Object of class bace_convergence
#' @param ... Additional arguments
#' @return Returns \code{x} invisibly; called for its side effect of printing
#'   the convergence assessment to the console.
#' @export
print.bace_convergence <- function(x, ...) {
  cat("\n=== BACE Chained-Equations Stability Check ===\n\n")

  if (is.na(x$converged)) {
    cat("Nothing could be assessed.\n")
    if (!is.null(x$message)) cat(x$message, "\n")
    if (!is.null(x$verdicts) && length(x$verdicts))
      cat("Variables seen:", paste(names(x$verdicts), collapse = ", "), "\n")
    return(invisible(x))
  }

  cat("Pathologies detected:", ifelse(x$converged, "none", "YES"), "\n")
  n_iter <- x$diagnostics$n_iter
  if (!is.null(n_iter)) {
    cat("Chained runs:", n_iter,
        sprintf("(first %d treated as warm-up)\n", x$diagnostics$warmup))
  }
  cat("\n")

  det <- x$method_results$summary$assessments
  if (!is.null(det)) {
    tab <- data.frame(
      variable = names(det),
      verdict  = vapply(det, function(a) a$verdict, character(1)),
      net_drift = round(vapply(det, function(a)
        unname(a$stats["net"]), numeric(1)), 2),
      step     = round(vapply(det, function(a)
        unname(a$stats["step"]), numeric(1)), 2),
      acf1     = round(vapply(det, function(a)
        unname(a$stats["acf1"]), numeric(1)), 2),
      note     = vapply(det, function(a)
        if (nzchar(a$note)) a$note else "", character(1)),
      row.names = NULL, check.names = FALSE)
    print(tab, row.names = FALSE)
    cat("\n(net_drift and step are in units of the variable's imputation",
        "spread;\n under settled chains, step has a noise floor of about",
        "sqrt(2/n_missing).)\n\n")
  }

  if (!x$converged) {
    v <- x$verdicts[!x$verdicts %in% c("ok", "unassessable")]
    for (var in names(v)) {
      remedy <- switch(v[[var]],
        drift = "still travelling; more chained runs (larger 'runs') can settle it",
        oscillation = paste0("flip-flops between runs; more runs will NOT fix this - ",
                             "look for variables feeding back on each other ",
                             "(derived or collinear predictors)"),
        stuck = paste0("imputation spread is collapsing; check this variable's ",
                       "model (priors, separation, near-constant values)"),
        diverging = paste0("imputation spread is growing without bound; check ",
                           "this variable's model (priors, scaling, separation)"),
        "could not be assessed; inspect the imputations for this variable")
      cat(" - ", var, ": ", remedy, "\n", sep = "")
    }
    cat("\n")
  }

  cat("A pass means NO PATHOLOGY DETECTED, not a certificate of convergence:\n")
  cat("no clear-cut test of chained-equations convergence exists (van Buuren\n")
  cat("2018, ch. 6.5.2). Inspect the traces: plot(<this object>).\n")
  invisible(x)
}
