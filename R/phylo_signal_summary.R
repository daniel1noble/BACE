#' @title phylo_signal_summary
#' @description Diagnostic summary of phylogenetic signal for every variable
#'   in a comparative dataset, fit as univariate phylogenetic mixed models with
#'   MCMCglmm. Returns posterior-mean phylogenetic heritability H² (plus 95%
#'   HPD), and optionally classical metrics (Pagel's lambda, Blomberg's K,
#'   Fritz-Purvis D). Designed as a cheap pre-flight check before committing
#'   to a full `bace()` imputation run.
#'
#' @section Caveat (important):
#' Phylogenetic signal is *necessary but not sufficient* for good imputation.
#' A high-H² trait can still impute poorly under MNAR missingness; a low-H²
#' trait can still impute well if tightly coupled to an observed predictor.
#' The table sets expectations; it does not predict imputation performance.
#'
#' @section Base models:
#' For each variable `v`, the function fits one univariate MCMCglmm with
#' `v ~ 1` as the fixed effect, and either `~Species` (single phylogenetic
#' random effect, `species = FALSE`) or `~Species + Species2` (dual random
#' effects: phylogenetic + non-phylogenetic species effect, `species = TRUE`).
#' The random-effect structure mirrors what BACE itself fits per variable, so
#' H² reported here is the same quantity BACE sees internally.
#'
#' Priors follow `.make_prior()`: gaussian/poisson use
#' `R = list(V = 1, nu = 2)` and Gelman 2006 parameter-expanded G; threshold
#' (binary / ordered) uses `R = list(V = 1, fix = 1)`; multinomial uses
#' `R = list(V = (I + J)/K, fix = 1)` (Hadfield's IJ parameterization).
#'
#' @section MCMC defaults:
#' If `nitt`/`burnin`/`thin` are `NULL`, type-specific defaults are used,
#' balanced to deliver 1500 post-burn draws with ESS > 1000 on every variance
#' component under typical data:
#'
#' | Type | nitt | burnin | thin |
#' |------|------|--------|------|
#' | gaussian | 20000 | 5000 | 10 |
#' | poisson | 30000 | 7500 | 15 |
#' | binary / ordered | 50000 | 12500 | 25 |
#' | multinomial | 80000 | 20000 | 40 |
#'
#' `quick = TRUE` halves nitt and burnin for a rough first-pass exploration
#' (post-burn draws ~750, ESS target ~500). Expect `low_ess` flags to trigger
#' more often in quick mode.
#'
#' @section Tree-size warnings:
#' Defaults are calibrated for 50–500 species. Smaller trees make phylogenetic
#' variance weakly identifiable (H² pulled toward prior mean ~0.5; HPD very
#' wide); larger trees have slower per-iteration MCMC. The function emits
#' `small_tree` (n_species < 30) and `large_tree` (n_species > 1000) warnings
#' automatically, both attached to the returned object's `$warnings` slot.
#'
#' @section Convergence diagnostics:
#' Every fit is checked with `coda::effectiveSize()` and `coda::geweke.diag()`
#' on all non-fixed variance components. Rows with `min_ess < min_ess`
#' threshold are flagged `low_ess`; rows with any `|z| > 2` are flagged
#' `geweke_fail`; either triggers `unreliable` in the `flag` column. Do not
#' interpret H² cells with either flag as data-driven estimates.
#'
#' @section H² formulas:
#' With V_A = phylo variance, V_S = non-phylo species variance (dual RE
#' only), V_R = residual variance, the latent-scale H² formulas are:
#'
#' * gaussian / poisson single RE: V_A / (V_A + V_R)
#' * gaussian / poisson dual RE:   V_A / (V_A + V_S + V_R)
#' * threshold / ordered single:   V_A / (V_A + 1)
#' * threshold / ordered dual:     V_A / (V_A + V_S + 1)
#' * multinomial single RE:        tr(G_phylo) / (tr(G_phylo) + tr(R))
#' * multinomial dual RE:          tr(G_phylo) / (tr(G_phylo) + tr(G_species) + tr(R))
#'
#' References: de Villemereuil et al. 2016 *Genetics* 204:1281; Nakagawa &
#' Schielzeth 2013 *Methods Ecol Evol* 4:133; Hadfield & Nakagawa 2010
#' *J Evol Biol* 23:494; Lynch 1991 *Evolution* 45:1065.
#'
#' @section Multinomial: H² versus lambda_nominal:
#' For `categorical` traits the table reports BOTH `H2_*` and
#' `lambda_nominal_*`. They measure related but DIFFERENT quantities and
#' typically disagree by 0.1-0.2 — that gap is expected, not a bug.
#'
#' * `H2_mean` is trace-based on MCMCglmm's working scale: residuals
#'   take their MCMCglmm-reported diagonal `2/(J+1)` per non-baseline
#'   level, with no rescaling. Faithful to the model's own variance
#'   decomposition; useful for within-categorical reporting.
#' * `lambda_nominal_mean` applies the Hadfield (2010 §3.7) / Amemiya
#'   (1981) c² correction to put each level on the standard probit
#'   scale (residual = 1 per level), then averages
#'   `lambda_k = G'_kk / (G'_kk + 1)` with
#'   `G'_kk = G_phylo[k,k] / (1 + c²·2/(J+1))`. This is the
#'   Pagel-comparable statistic — same scale and interpretation as
#'   Pagel's lambda for continuous traits, so it can be compared
#'   across trait types in a multi-trait phylo_signal table.
#'
#' Use `lambda_nominal_mean` for cross-trait comparisons; use
#' `H2_mean` when you specifically want MCMCglmm-native variance
#' partitioning. Per-level (one value per non-baseline category)
#' detail is on `attr(model, "lambda_nominal_per_level")` after a
#' run with `keep_models = TRUE`. Reference: Ayumi (2024)
#' multinomial GLMM tutorial,
#' \url{https://ayumi-495.github.io/multinomial-GLMM-tutorial/#nominal}.
#'
#' @param data data.frame with rows keyed by the species column.
#' @param tree phylogenetic tree of class `phylo` (ape package).
#' @param species_col character; column name in `data` matching `tree$tip.label`.
#' @param variables character vector; which columns to summarise. Default:
#'   all columns except `species_col`.
#' @param species logical; `FALSE` (default) = single phylogenetic random
#'   effect; `TRUE` = dual (phylo + non-phylo species). Dual requires
#'   within-species replication.
#' @param methods character; which metrics to compute. Default `"auto"` picks
#'   per-type (H² always; lambda/K for gaussian; D for binary). Or pass
#'   `c("H2","lambda","K","D")`.
#' @param ovr_categorical logical; if `TRUE`, multinomial variables are fit as
#'   J binary threshold models (one-vs-rest) instead of a multinomial probit.
#'   Matches BACE's own OVR path.
#' @param nitt,burnin,thin integer; MCMC settings. `NULL` (default) uses
#'   type-specific defaults (see Details).
#' @param quick logical; if `TRUE`, halves nitt and burnin for a rough
#'   first-pass exploration. Default `FALSE`.
#' @param prior list; MCMCglmm prior. `NULL` (default) uses `.make_prior()`.
#' @param n_sim integer; simulations for Blomberg's K p-value. Default 999.
#' @param min_ess numeric; minimum acceptable effective sample size on any
#'   variance component. Rows below this get `low_ess` flag. Default 1000.
#' @param keep_models logical; if `TRUE` (default), return the fitted MCMCglmm
#'   objects on `$models`. For `ovr_categorical = TRUE` runs, the per-level
#'   binary fits are returned as a named list keyed by factor level. Set to
#'   `FALSE` for memory-constrained runs (e.g. many variables x large trees).
#' @param verbose logical; print progress messages. Default `TRUE`.
#'
#' @return An S3 `phylo_signal` object with slots `$table`, `$models`,
#'   `$warnings`, `$call`, `$species`, `$n_sim`. Print with
#'   `print(x)` for a formatted table.
#'
#' @importFrom stats as.formula complete.cases setNames
#' @importFrom coda effectiveSize geweke.diag HPDinterval as.mcmc
#' @export
phylo_signal_summary <- function(
    data,
    tree,
    species_col     = "Species",
    variables       = NULL,
    species         = FALSE,
    methods         = "auto",
    ovr_categorical = FALSE,
    nitt            = NULL,
    burnin          = NULL,
    thin            = NULL,
    quick           = FALSE,
    prior           = NULL,
    n_sim           = 999,
    min_ess         = 1000,
    keep_models     = TRUE,
    verbose         = TRUE) {

  # ----- Validate inputs -----
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be a phylo object (ape package).")
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }
  if (!species_col %in% colnames(data)) {
    stop("`species_col` ('", species_col, "') not found in data.")
  }
  if (is.null(variables)) {
    variables <- setdiff(colnames(data), species_col)
  }
  missing_vars <- setdiff(variables, colnames(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  if (species) {
    counts <- table(data[[species_col]])
    n_replicated <- sum(counts > 1)
    if (n_replicated < 2) {
      stop("`species = TRUE` requires within-species replicates (at least 2 ",
           "species with > 1 row). Found ", n_replicated, ". ",
           "Use `species = FALSE` for one-row-per-species data.")
    }
  }

  # ----- Tree-size warnings -----
  warnings_list <- list()
  n_species <- length(tree$tip.label)
  if (n_species < 30) {
    w <- paste0(
      "Tree has ", n_species, " species (<30). Phylogenetic variance is ",
      "weakly identified; H2 posterior is prior-sensitive. Interpret with ",
      "caution, and consider the HPD width more than the posterior mean."
    )
    warning(w, call. = FALSE)
    warnings_list$small_tree <- w
  }
  if (n_species > 1000) {
    scale_factor <- sqrt(n_species / 500)
    w <- sprintf(
      paste0("Tree has %d species (>1000). Default MCMC may underrun; ",
             "consider increasing nitt by sqrt(n_species/500) ~ %.2fx ",
             "or using quick=TRUE for a rough first pass."),
      n_species, scale_factor
    )
    warning(w, call. = FALSE)
    warnings_list$large_tree <- w
  }
  if (quick) {
    w <- "quick = TRUE: nitt and burnin halved. Expect more low_ess flags."
    warning(w, call. = FALSE)
    warnings_list$quick_mode <- w
  }

  # ----- Detect types per variable -----
  var_types <- .detect_types_for_signal(data, variables)
  if (verbose) {
    message("Variable types detected: ",
            paste(sprintf("%s=%s", variables, unlist(var_types[variables])),
                  collapse = ", "))
  }

  # ----- Phylo inverse matrix (computed once) -----
  A <- MCMCglmm::inverseA(tree, nodes = "ALL")$Ainv

  # ----- Fit per-variable models -----
  rows   <- vector("list", length(variables))
  models <- if (keep_models) vector("list", length(variables)) else NULL
  names(rows) <- variables
  if (!is.null(models)) names(models) <- variables

  for (i in seq_along(variables)) {
    v <- variables[[i]]
    type_v <- var_types[[v]]
    if (verbose) message(sprintf("[%d/%d] Fitting '%s' (type=%s)...",
                                 i, length(variables), v, type_v))

    res <- tryCatch(
      .fit_signal_one(
        v = v, type = type_v, data = data, tree = tree, A = A,
        species_col = species_col, species = species,
        nitt = nitt, burnin = burnin, thin = thin, quick = quick,
        prior = prior, min_ess = min_ess, methods = methods,
        ovr_categorical = ovr_categorical, n_sim = n_sim,
        keep_models = keep_models
      ),
      error = function(e) {
        warning(sprintf("Fit for '%s' failed: %s", v, conditionMessage(e)),
                call. = FALSE)
        NULL
      }
    )

    if (is.null(res)) {
      rows[[v]] <- .empty_signal_row(v, type_v, species)
      rows[[v]]$flag <- "fit_error"
    } else {
      rows[[v]] <- res$row
      if (keep_models) models[[v]] <- res$model
    }
  }

  tbl <- do.call(rbind, rows)
  rownames(tbl) <- NULL

  out <- list(
    table    = tbl,
    models   = models,
    warnings = warnings_list,
    call     = match.call(),
    species  = species,
    n_sim    = n_sim
  )
  class(out) <- "phylo_signal"
  out
}


# =============================================================================
# Internal helpers
# =============================================================================

#' Per-variable type detection wrapper for phylo_signal_summary.
#' Reuses .get_type()'s logic via its underlying summary.
#' @keywords internal
#' @noRd
.detect_types_for_signal <- function(data, variables) {
  # Build the per-variable summary .get_type expects
  sum_tbl <- do.call(rbind, lapply(variables, function(v) {
    x <- data[[v]]
    data.frame(
      variable   = v,
      is_numeric = is.numeric(x),
      is_integer = is.integer(x) ||
                   (is.numeric(x) && all(!is.na(x) & x == floor(x))),
      is_factor  = is.factor(x),
      is_ordered = is.ordered(x),
      n_levels   = if (is.factor(x)) nlevels(x) else NA_integer_,
      stringsAsFactors = FALSE
    )
  }))
  out <- stats::setNames(
    lapply(variables, function(v) .get_type(sum_tbl, v, data)),
    variables
  )
  # Fallback: if .get_type returned NULL (unclassifiable), mark as gaussian
  # for numeric, else skip with warning.
  for (v in variables) {
    if (is.null(out[[v]])) {
      if (is.numeric(data[[v]])) out[[v]] <- "gaussian"
      else stop("Could not classify variable '", v,
                "'. Check its class (numeric/integer/factor/ordered).")
    }
  }
  out
}


#' Type-specific MCMC defaults, with `quick` halving.
#' @keywords internal
#' @noRd
.resolve_mcmc_defaults <- function(type, nitt, burnin, thin, quick) {
  defaults <- switch(
    type,
    gaussian    = list(nitt = 20000, burnin = 5000,  thin = 10),
    poisson     = list(nitt = 30000, burnin = 7500,  thin = 15),
    threshold   = list(nitt = 50000, burnin = 12500, thin = 25),
    ordinal     = list(nitt = 50000, burnin = 12500, thin = 25),
    categorical = list(nitt = 80000, burnin = 20000, thin = 40),
    stop("Unknown type for MCMC defaults: '", type, "'")
  )
  if (quick) {
    defaults$nitt   <- defaults$nitt   / 2L
    defaults$burnin <- defaults$burnin / 2L
  }
  if (!is.null(nitt))   defaults$nitt   <- nitt
  if (!is.null(burnin)) defaults$burnin <- burnin
  if (!is.null(thin))   defaults$thin   <- thin
  defaults
}


#' Fit one signal model and return a single-row data.frame.
#' @keywords internal
#' @noRd
.fit_signal_one <- function(v, type, data, tree, A, species_col, species,
                            nitt, burnin, thin, quick, prior, min_ess,
                            methods, ovr_categorical, n_sim,
                            keep_models = TRUE) {

  # --- Subset to rows with observed v and species in tree ---
  d <- data[!is.na(data[[v]]), c(v, species_col), drop = FALSE]
  d <- d[d[[species_col]] %in% tree$tip.label, , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No observed rows for variable '", v,
         "' after filtering to species in tree.")
  }
  # Ensure species column is a factor with levels from tree
  d[[species_col]] <- factor(d[[species_col]], levels = tree$tip.label)
  if (species) {
    d[[paste0(species_col, "2")]] <- d[[species_col]]
  }
  # z-scale gaussian responses: matches BACE's .data_prep() convention and
  # keeps priors (V = 1, nu = 2) well-matched to the data scale. H² is a
  # scale-invariant ratio so this does not bias the quantity we report.
  # Poisson counts are NOT scaled (breaks the integer-count family).
  if (type == "gaussian") {
    d[[v]] <- as.numeric(scale(as.numeric(d[[v]])))
  }

  # --- MCMC settings ---
  mcmc <- .resolve_mcmc_defaults(type, nitt, burnin, thin, quick)

  # --- Prior ---
  n_rand <- if (species) 2L else 1L
  n_levels <- if (type %in% c("threshold", "ordinal", "categorical") &&
                  is.factor(d[[v]])) nlevels(d[[v]]) else NULL
  pr <- prior
  if (is.null(pr)) {
    pr <- .make_prior(
      n_rand = n_rand, type = type,
      n_levels = n_levels, par_expand = TRUE
    )
  }

  # --- Formulas ---
  fixform <- stats::as.formula(paste(v, "~ 1"))
  if (species) {
    # .model_fit() expects a named list for dual RE with elements
    # 'phylo' (gets the ginverse) and 'species' (identity ginverse).
    randform <- list(
      phylo   = stats::as.formula(paste0("~ ", species_col)),
      species = stats::as.formula(paste0("~ ", species_col))
    )
  } else {
    randform <- stats::as.formula(paste0("~ ", species_col))
  }

  # --- OVR multinomial shortcut ---
  if (type == "categorical" && ovr_categorical) {
    return(.fit_signal_ovr(
      v = v, d = d, type = type, tree = tree, A = A,
      species_col = species_col, species = species,
      mcmc = mcmc, n_levels = n_levels, min_ess = min_ess,
      methods = methods, n_sim = n_sim, data_full = data,
      keep_models = keep_models
    ))
  }

  # --- Fit via .model_fit() (handles trait-expansion for categorical,
  # ginverse for dual RE, and slice sampling automatically) ---
  model <- .model_fit(
    data       = d,
    tree       = tree,
    fixformula = fixform,
    randformula = randform,
    type       = type,
    prior      = pr,
    nitt       = mcmc$nitt,
    burnin     = mcmc$burnin,
    thin       = mcmc$thin
  )

  # --- Compute H² per posterior draw ---
  h2_stats <- .compute_H2_stats(model, type, species, species_col, n_levels)

  # --- Per-level multinomial Pagel-style lambda (categorical only) ---
  lam_nom <- if (type == "categorical") {
    .compute_lambda_nominal(model, species_col, n_levels)
  } else {
    list(avg_mean = NA_real_, avg_lo = NA_real_, avg_hi = NA_real_,
         per_level = NULL)
  }

  # --- Compute R_species (dual RE only) ---
  rs_stats <- if (species) {
    .compute_R_species_stats(model, type, species_col, n_levels)
  } else {
    list(mean = NA_real_, lo = NA_real_, hi = NA_real_)
  }

  # --- Convergence diagnostics ---
  diag_cols <- .get_vc_diag_cols(model, type, species, species_col)
  ess_vals  <- if (length(diag_cols) > 0) {
    coda::effectiveSize(coda::as.mcmc(as.matrix(model$VCV)[, diag_cols,
                                                           drop = FALSE]))
  } else NA_real_
  min_ess_val <- if (all(is.na(ess_vals))) NA_real_ else min(ess_vals, na.rm = TRUE)
  geweke_fail <- .geweke_fail(model, diag_cols)

  # --- Classical metrics ---
  classical <- .compute_classical_metrics(
    v = v, type = type, data = data, tree = tree,
    species_col = species_col, methods = methods, n_sim = n_sim
  )

  # --- Flag ---
  n_obs <- nrow(d)
  flag_vec <- character(0)
  if (!is.na(min_ess_val) && min_ess_val < min_ess) flag_vec <- c(flag_vec, "low_ess")
  if (isTRUE(geweke_fail)) flag_vec <- c(flag_vec, "geweke_fail")
  if (length(flag_vec) > 0) flag_vec <- c(flag_vec, "unreliable")
  if (n_obs < 20) flag_vec <- c(flag_vec, "low_n")
  flag_str <- if (length(flag_vec) == 0) "" else paste(unique(flag_vec), collapse = ";")

  # --- Row ---
  row <- data.frame(
    variable       = v,
    type           = type,
    n_obs          = n_obs,
    H2_mean        = h2_stats$mean,
    H2_lo          = h2_stats$lo,
    H2_hi          = h2_stats$hi,
    lambda_nominal_mean = lam_nom$avg_mean,
    lambda_nominal_lo   = lam_nom$avg_lo,
    lambda_nominal_hi   = lam_nom$avg_hi,
    R_species_mean = rs_stats$mean,
    R_species_lo   = rs_stats$lo,
    R_species_hi   = rs_stats$hi,
    lambda         = classical$lambda,
    lambda_p       = classical$lambda_p,
    K              = classical$K,
    K_p            = classical$K_p,
    D              = classical$D,
    D_random_p     = classical$D_random_p,
    D_BM_p         = classical$D_BM_p,
    min_ess        = min_ess_val,
    interpretation = .interpret_H2(h2_stats$mean),
    flag           = flag_str,
    stringsAsFactors = FALSE
  )

  # Stash per-level lambda detail on the model object for later access
  if (!is.null(lam_nom$per_level)) {
    attr(model, "lambda_nominal_per_level") <- lam_nom$per_level
  }
  list(row = row, model = model)
}


#' Fit J binary threshold models for a multinomial variable (OVR path).
#' @keywords internal
#' @noRd
.fit_signal_ovr <- function(v, d, type, tree, A, species_col, species,
                             mcmc, n_levels, min_ess, methods, n_sim,
                             data_full, keep_models = TRUE) {
  lvls <- levels(d[[v]])
  per_level_h2 <- numeric(length(lvls))
  per_level_ess <- numeric(length(lvls))
  per_level_gwk <- logical(length(lvls))
  per_level_fits <- if (keep_models) vector("list", length(lvls)) else NULL
  if (!is.null(per_level_fits)) names(per_level_fits) <- lvls
  n_obs_total <- nrow(d)

  # Build threshold prior once
  pr <- .make_prior(
    n_rand = if (species) 2L else 1L,
    type   = "threshold",
    n_levels = 2, par_expand = TRUE
  )
  randform <- if (species) {
    stats::as.formula(paste0("~ ", species_col, " + ", species_col, "2"))
  } else {
    stats::as.formula(paste0("~ ", species_col))
  }
  ginv <- stats::setNames(list(A), species_col)

  for (j in seq_along(lvls)) {
    lvl <- lvls[[j]]
    d_j <- d
    d_j[[paste0(v, "_ovr")]] <- factor(ifelse(d[[v]] == lvl, 1L, 0L),
                                       levels = c(0L, 1L))
    fixform_j <- stats::as.formula(paste0(v, "_ovr ~ 1"))

    m_j <- MCMCglmm::MCMCglmm(
      fixed    = fixform_j,
      random   = randform,
      data     = d_j,
      family   = "threshold",
      ginverse = ginv,
      prior    = pr,
      nitt     = mcmc$nitt,
      burnin   = mcmc$burnin,
      thin     = mcmc$thin,
      slice    = TRUE,
      verbose  = FALSE
    )
    h2_s <- .compute_H2_stats(m_j, "threshold", species, species_col, NULL)
    per_level_h2[j] <- h2_s$mean

    diag_cols <- .get_vc_diag_cols(m_j, "threshold", species, species_col)
    ess <- coda::effectiveSize(coda::as.mcmc(
      as.matrix(m_j$VCV)[, diag_cols, drop = FALSE]))
    per_level_ess[j] <- min(ess, na.rm = TRUE)
    per_level_gwk[j] <- .geweke_fail(m_j, diag_cols)
    if (keep_models) per_level_fits[[j]] <- m_j
  }

  h2_mean <- mean(per_level_h2)
  min_ess_val <- min(per_level_ess, na.rm = TRUE)

  flag_vec <- character(0)
  if (!is.na(min_ess_val) && min_ess_val < min_ess) flag_vec <- c(flag_vec, "low_ess")
  if (any(per_level_gwk)) flag_vec <- c(flag_vec, "geweke_fail")
  if (length(flag_vec) > 0) flag_vec <- c(flag_vec, "unreliable")
  if (n_obs_total < 20) flag_vec <- c(flag_vec, "low_n")
  flag_str <- if (length(flag_vec) == 0) "" else paste(unique(flag_vec), collapse = ";")

  classical <- .compute_classical_metrics(
    v = v, type = "categorical", data = data_full, tree = tree,
    species_col = species_col, methods = methods, n_sim = n_sim
  )

  row <- data.frame(
    variable       = v,
    type           = "categorical (OVR)",
    n_obs          = n_obs_total,
    H2_mean        = h2_mean,
    H2_lo          = NA_real_,  # OVR aggregation: per-level HPDs not pooled
    H2_hi          = NA_real_,
    lambda_nominal_mean = NA_real_,  # OVR path uses per-level threshold fits
    lambda_nominal_lo   = NA_real_,  # not the multinomial probit; lambda_nominal
    lambda_nominal_hi   = NA_real_,  # is defined only on the joint G_phylo
    R_species_mean = NA_real_,
    R_species_lo   = NA_real_,
    R_species_hi   = NA_real_,
    lambda         = classical$lambda,
    lambda_p       = classical$lambda_p,
    K              = classical$K,
    K_p            = classical$K_p,
    D              = classical$D,
    D_random_p     = classical$D_random_p,
    D_BM_p         = classical$D_BM_p,
    min_ess        = min_ess_val,
    interpretation = .interpret_H2(h2_mean),
    flag           = flag_str,
    stringsAsFactors = FALSE
  )
  list(row = row, model = per_level_fits)
}


#' Posterior H² draws and summaries, depending on type and RE structure.
#' @keywords internal
#' @noRd
.compute_H2_stats <- function(model, type, species, species_col, n_levels) {
  vcv <- as.matrix(model$VCV)
  cols <- colnames(vcv)
  phylo_col  <- species_col
  sp2_col    <- paste0(species_col, "2")
  units_col  <- "units"

  if (type %in% c("gaussian", "poisson")) {
    V_A <- vcv[, phylo_col]
    V_R <- vcv[, units_col]
    V_S <- if (species) vcv[, sp2_col] else 0
    h2 <- V_A / (V_A + V_S + V_R)

  } else if (type %in% c("threshold", "ordinal")) {
    V_A <- vcv[, phylo_col]
    V_S <- if (species) vcv[, sp2_col] else 0
    # Residual is fixed at 1 (probit identification)
    h2 <- V_A / (V_A + V_S + 1)

  } else if (type == "categorical") {
    # Multinomial: VCV columns for G_phylo are (K-1)^2 entries named like
    # "traitY.1:traitY.1.Species" etc. Sum the diagonal entries
    # (positions where trait_i == trait_j) for the phylo block.
    J <- n_levels - 1L
    phylo_diag <- .get_multinom_diag_cols(cols, phylo_col, J)
    tr_Gphylo <- rowSums(vcv[, phylo_diag, drop = FALSE])

    if (species) {
      sp2_diag  <- .get_multinom_diag_cols(cols, sp2_col, J)
      tr_Gsp    <- rowSums(vcv[, sp2_diag, drop = FALSE])
    } else {
      tr_Gsp <- 0
    }
    units_diag <- .get_multinom_diag_cols(cols, "units", J)
    tr_R <- if (length(units_diag) > 0) {
      rowSums(vcv[, units_diag, drop = FALSE])
    } else {
      # Fallback: residual fixed per IJ parameterization
      (J * J) / n_levels
    }

    denom <- tr_Gphylo + tr_Gsp + tr_R
    h2 <- tr_Gphylo / denom

  } else {
    stop("Unknown type in .compute_H2_stats: ", type)
  }

  hpd <- tryCatch(
    coda::HPDinterval(coda::as.mcmc(h2)),
    error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1,
                               dimnames = list(NULL, c("lower", "upper")))
  )
  list(mean = mean(h2), lo = hpd[1, "lower"], hi = hpd[1, "upper"])
}


#' Per-level Pagel-style lambda for multinomial categorical traits.
#'
#' For each non-baseline category k, computes
#'   lambda_k = (G_phylo[k,k] / (1 + c_k)) / (G_phylo[k,k]/(1 + c_k) + 1)
#' where c_k = c2 * IJ[k,k] = c2 * 2/(J+1) is the Amemiya c² correction
#' (Amemiya 1981) and J = K-1 is the number of non-baseline levels,
#' c2 = (16*sqrt(3)/(15*pi))^2 ≈ 0.3458.
#'
#' Each lambda_k has the same Pagel interpretation as for continuous
#' traits: fraction of LATENT-scale variance for category k that is
#' phylogenetic. Useful when the per-category interpretation matters
#' (some categories may have strong phylogenetic signal, others weak).
#' Differs from the trace-based H² returned alongside, which is one
#' summary number across the whole multinomial.
#'
#' Returns the average across non-baseline levels plus per-level
#' details. Requires MCMCglmm `family = "categorical"` (multinomial
#' probit with R fixed at I/(J+1) + 1/(J+1)*1*1').
#'
#' Reference: Ayumi (2024) multinomial GLMM tutorial,
#' \url{https://ayumi-495.github.io/multinomial-GLMM-tutorial/#nominal};
#' Hadfield (2010) §3.7; Amemiya (1981) Econometrica 49:1483.
#' @keywords internal
#' @noRd
.compute_lambda_nominal <- function(model, species_col, n_levels) {
  J <- n_levels - 1L
  if (J < 1L) {
    return(list(avg_mean = NA_real_, avg_lo = NA_real_, avg_hi = NA_real_,
                per_level = NULL))
  }
  vcv <- as.matrix(model$VCV)
  cols <- colnames(vcv)
  phylo_diag <- .get_multinom_diag_cols(cols, species_col, J)
  if (length(phylo_diag) == 0L) {
    return(list(avg_mean = NA_real_, avg_lo = NA_real_, avg_hi = NA_real_,
                per_level = NULL))
  }

  # Amemiya c² with the IJ-diagonal correction. IJ = (1/(J+1))(I + 11')
  # so diag(IJ) = 2/(J+1) for every k.
  c2      <- (16 * sqrt(3) / (15 * pi))^2
  c2_corr <- c2 * 2 / (J + 1)

  per_level <- vector("list", length(phylo_diag))
  lambda_avg <- numeric(nrow(vcv))
  for (k in seq_along(phylo_diag)) {
    G_kk    <- vcv[, phylo_diag[k]]
    G_corr  <- G_kk / (1 + c2_corr)
    lam_k   <- G_corr / (G_corr + 1)
    hpd_k <- tryCatch(
      coda::HPDinterval(coda::as.mcmc(lam_k)),
      error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1,
                                  dimnames = list(NULL, c("lower", "upper")))
    )
    per_level[[k]] <- list(
      column = phylo_diag[k],
      mean   = mean(lam_k),
      lo     = hpd_k[1, "lower"],
      hi     = hpd_k[1, "upper"]
    )
    lambda_avg <- lambda_avg + lam_k
  }
  lambda_avg <- lambda_avg / length(phylo_diag)
  hpd_avg <- tryCatch(
    coda::HPDinterval(coda::as.mcmc(lambda_avg)),
    error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1,
                                dimnames = list(NULL, c("lower", "upper")))
  )
  list(
    avg_mean  = mean(lambda_avg),
    avg_lo    = hpd_avg[1, "lower"],
    avg_hi    = hpd_avg[1, "upper"],
    per_level = per_level
  )
}


#' Find the diagonal column names of a (K-1) x (K-1) VCV block for a given RE.
#' MCMCglmm names them like `traitY.1:traitY.1.<RE>` where the two trait
#' indices must match.
#' @keywords internal
#' @noRd
.get_multinom_diag_cols <- function(cols, re_name, J) {
  # Match "traitXXX.i:traitXXX.i.<re_name>" pattern
  pat <- sprintf("trait[^.]+\\.(\\d+):trait[^.]+\\.\\1\\.%s$",
                 gsub("([^[:alnum:]])", "\\\\\\1", re_name))
  diag_idx <- grep(pat, cols, perl = TRUE)
  if (length(diag_idx) == 0) {
    # Fallback: any column ending in .<re_name> — caller may aggregate more
    # than just diagonals, which overestimates. Emit a warning.
    fallback <- grep(paste0("\\.", gsub("([^[:alnum:]])", "\\\\\\1", re_name),
                            "$"), cols)
    if (length(fallback) == 0) return(character(0))
    # Take only the first J columns as a heuristic diagonal guess.
    return(cols[fallback[seq_len(min(J, length(fallback)))]])
  }
  cols[diag_idx]
}


#' Posterior R_species (non-phylo species variance fraction) for dual RE.
#' @keywords internal
#' @noRd
.compute_R_species_stats <- function(model, type, species_col, n_levels) {
  vcv <- as.matrix(model$VCV)
  cols <- colnames(vcv)
  phylo_col <- species_col
  sp2_col   <- paste0(species_col, "2")

  if (type %in% c("gaussian", "poisson")) {
    V_A <- vcv[, phylo_col]
    V_S <- vcv[, sp2_col]
    V_R <- vcv[, "units"]
    rs <- V_S / (V_A + V_S + V_R)

  } else if (type %in% c("threshold", "ordinal")) {
    V_A <- vcv[, phylo_col]
    V_S <- vcv[, sp2_col]
    rs <- V_S / (V_A + V_S + 1)

  } else if (type == "categorical") {
    J <- n_levels - 1L
    tr_Gphylo <- rowSums(vcv[, .get_multinom_diag_cols(cols, phylo_col, J),
                             drop = FALSE])
    tr_Gsp    <- rowSums(vcv[, .get_multinom_diag_cols(cols, sp2_col, J),
                             drop = FALSE])
    units_diag <- .get_multinom_diag_cols(cols, "units", J)
    tr_R <- if (length(units_diag) > 0) {
      rowSums(vcv[, units_diag, drop = FALSE])
    } else {
      (J * J) / n_levels
    }
    rs <- tr_Gsp / (tr_Gphylo + tr_Gsp + tr_R)

  } else {
    rs <- NA_real_
  }

  hpd <- tryCatch(
    coda::HPDinterval(coda::as.mcmc(rs)),
    error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1,
                               dimnames = list(NULL, c("lower", "upper")))
  )
  list(mean = mean(rs), lo = hpd[1, "lower"], hi = hpd[1, "upper"])
}


#' Get VCV diagonal column names for ESS/Geweke checks (excluding fixed).
#' @keywords internal
#' @noRd
.get_vc_diag_cols <- function(model, type, species, species_col) {
  vcv_cols <- colnames(as.matrix(model$VCV))
  phylo_col <- species_col
  sp2_col   <- paste0(species_col, "2")

  if (type %in% c("gaussian", "poisson")) {
    cols <- c(phylo_col, if (species) sp2_col, "units")
  } else if (type %in% c("threshold", "ordinal")) {
    # Residual fixed at 1: drop it.
    cols <- c(phylo_col, if (species) sp2_col)
  } else if (type == "categorical") {
    # Diagonal cols of G blocks; drop R (fixed IJ).
    J <- NA
    # Determine J from the actual VCV structure
    J_guess <- length(.get_multinom_diag_cols(vcv_cols, phylo_col, J = 10L))
    cols <- c(.get_multinom_diag_cols(vcv_cols, phylo_col, J_guess),
              if (species) .get_multinom_diag_cols(vcv_cols, sp2_col, J_guess))
  } else {
    cols <- character(0)
  }
  cols[cols %in% vcv_cols]
}


#' Geweke diagnostic on given VCV columns; returns TRUE if any |z| > 2.
#' @keywords internal
#' @noRd
.geweke_fail <- function(model, diag_cols) {
  if (length(diag_cols) == 0) return(NA)
  mc <- coda::as.mcmc(as.matrix(model$VCV)[, diag_cols, drop = FALSE])
  gw <- tryCatch(coda::geweke.diag(mc), error = function(e) NULL)
  if (is.null(gw)) return(NA)
  z <- if (is.list(gw)) gw$z else gw
  any(abs(z) > 2, na.rm = TRUE)
}


#' Classical phylogenetic signal metrics (lambda, K, D) with soft package deps.
#' @keywords internal
#' @noRd
.compute_classical_metrics <- function(v, type, data, tree, species_col,
                                        methods, n_sim) {
  out <- list(
    lambda = NA_real_, lambda_p = NA_real_,
    K = NA_real_, K_p = NA_real_,
    D = NA_real_, D_random_p = NA_real_, D_BM_p = NA_real_
  )
  auto <- identical(methods, "auto")
  want <- function(m) auto || (m %in% methods)

  # Trait vector aligned to tree tip labels (one value per tip;
  # for multi-row-per-species data, use the first non-NA value).
  trait <- .tip_vector(data, v, species_col, tree$tip.label, type)
  if (all(is.na(trait))) return(out)

  # --- lambda + K for continuous ---
  if (type %in% c("gaussian", "poisson") && (want("lambda") || want("K"))) {
    if (requireNamespace("phytools", quietly = TRUE)) {
      has <- !is.na(trait)
      if (sum(has) >= 5) {
        tr_sub <- ape::keep.tip(tree, tree$tip.label[has])
        trait_sub <- trait[has]
        names(trait_sub) <- tree$tip.label[has]
        if (want("lambda")) {
          lam <- tryCatch(
            phytools::phylosig(tr_sub, trait_sub,
                                method = "lambda", test = TRUE),
            error = function(e) NULL
          )
          if (!is.null(lam)) {
            out$lambda   <- as.numeric(lam$lambda)
            out$lambda_p <- as.numeric(lam$P)
          }
        }
        if (want("K")) {
          kk <- tryCatch(
            phytools::phylosig(tr_sub, trait_sub,
                                method = "K", test = TRUE,
                                nsim = n_sim),
            error = function(e) NULL
          )
          if (!is.null(kk)) {
            out$K   <- as.numeric(kk$K)
            out$K_p <- as.numeric(kk$P)
          }
        }
      }
    }
  }

  # --- D for binary ---
  if (type == "threshold" && want("D")) {
    # D only meaningful for binary; check n_levels.
    is_bin <- is.factor(data[[v]]) && nlevels(data[[v]]) == 2
    if (is_bin && requireNamespace("caper", quietly = TRUE)) {
      df_d <- data.frame(
        species = tree$tip.label,
        val     = trait,
        stringsAsFactors = FALSE
      )
      df_d <- df_d[!is.na(df_d$val), , drop = FALSE]
      if (nrow(df_d) >= 5) {
        tr_sub <- ape::keep.tip(tree, df_d$species)
        df_d$val <- as.integer(as.character(df_d$val))
        # caper::phylo.d uses NSE for names.col and binvar. Go through
        # do.call + as.name to avoid the global-variable R CMD check note.
        res <- tryCatch(
          do.call(caper::phylo.d, list(
            data      = df_d,
            phy       = tr_sub,
            names.col = as.name("species"),
            binvar    = as.name("val"),
            permut    = n_sim
          )),
          error = function(e) NULL
        )
        if (!is.null(res)) {
          out$D          <- as.numeric(res$DEstimate)
          out$D_random_p <- as.numeric(res$Pval1)
          out$D_BM_p     <- as.numeric(res$Pval0)
        }
      }
    }
  }
  out
}


#' Collapse a data column to one value per tree tip (first non-NA).
#' @keywords internal
#' @noRd
.tip_vector <- function(data, v, species_col, tip_labels, type) {
  idx <- match(as.character(data[[species_col]]), tip_labels)
  out <- rep(NA, length(tip_labels))
  for (k in seq_along(tip_labels)) {
    rows <- which(idx == k)
    if (length(rows) == 0) next
    vals <- data[[v]][rows]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) next
    if (type %in% c("gaussian", "poisson")) {
      out[k] <- mean(as.numeric(vals))
    } else {
      # factor: first observed
      out[k] <- as.character(vals[[1]])
    }
  }
  out
}


#' Interpretation band text.
#' @keywords internal
#' @noRd
.interpret_H2 <- function(h2) {
  if (is.na(h2)) return(NA_character_)
  if (h2 < 0.2) "low"
  else if (h2 < 0.5) "moderate"
  else "high"
}


#' Empty row used when a fit errors out.
#' @keywords internal
#' @noRd
.empty_signal_row <- function(v, type, species) {
  data.frame(
    variable = v, type = type, n_obs = NA_integer_,
    H2_mean = NA_real_, H2_lo = NA_real_, H2_hi = NA_real_,
    lambda_nominal_mean = NA_real_,
    lambda_nominal_lo   = NA_real_,
    lambda_nominal_hi   = NA_real_,
    R_species_mean = NA_real_, R_species_lo = NA_real_,
    R_species_hi = NA_real_,
    lambda = NA_real_, lambda_p = NA_real_, K = NA_real_, K_p = NA_real_,
    D = NA_real_, D_random_p = NA_real_, D_BM_p = NA_real_,
    min_ess = NA_real_, interpretation = NA_character_, flag = "",
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# S3 print
# =============================================================================

#' @export
print.phylo_signal <- function(x, ...) {
  cat("=============================================================\n")
  cat("  Phylogenetic signal summary\n")
  cat("=============================================================\n")
  cat("  Random-effect structure: ",
      if (isTRUE(x$species)) "dual (phylo + non-phylo species)" else "single (phylo only)",
      "\n", sep = "")
  cat("  Variables                :", nrow(x$table), "\n")
  if (length(x$warnings) > 0) {
    cat("\n--- Warnings ---\n")
    for (nm in names(x$warnings)) {
      cat("  [", nm, "] ", x$warnings[[nm]], "\n", sep = "")
    }
  }

  cat("\n--- Table ---\n")
  tbl <- x$table

  # Build a human-friendly display sub-table. Hide empty classical columns.
  disp_cols <- c("variable", "type", "n_obs",
                 "H2_mean", "H2_lo", "H2_hi")
  if (isTRUE(x$species)) {
    disp_cols <- c(disp_cols, "R_species_mean")
  }
  for (mcol in c("lambda", "K", "D")) {
    if (any(!is.na(tbl[[mcol]]))) disp_cols <- c(disp_cols, mcol)
  }
  disp_cols <- c(disp_cols, "min_ess", "interpretation", "flag")

  disp <- tbl[, disp_cols, drop = FALSE]

  # Mark unreliable rows with asterisk on H2_mean
  unrel <- grepl("unreliable", tbl$flag)
  if (any(unrel)) {
    disp$H2_mean[unrel] <- paste0(
      formatC(disp$H2_mean[unrel], format = "f", digits = 3), "*")
  }

  print(disp, row.names = FALSE, digits = 3)

  cat("\nInterpretation bands: H2 < 0.2 = low, 0.2-0.5 = moderate, >0.5 = high.\n")
  cat("Flag meanings:\n")
  cat("  low_ess     -- effective sample size below threshold (default 1000)\n")
  cat("  geweke_fail -- Geweke diagnostic |z| > 2 on any variance component\n")
  cat("  unreliable  -- H2 cell should NOT be trusted (low_ess or geweke_fail)\n")
  cat("  low_n       -- n_obs < 20; posterior is data-sparse\n")
  cat("  fit_error   -- MCMCglmm fit failed (check warnings)\n")
  cat("\n")
  cat("Note: phylogenetic signal is necessary but NOT sufficient for good\n")
  cat("imputation. A high-H2 trait can still impute poorly under MNAR missingness;\n")
  cat("a low-H2 trait can still impute well if coupled to an observed predictor.\n")
  invisible(x)
}
