# =================================================================================================
# STEP 2 of 4 -- one simulated dataset from one scenario's settings.
#
# Two halves:
#   simulate_traits()   calls BACE::sim_bace() to make the traits and the tree, with the truth kept
#   apply_missingness() hides values in one of three ways: at random, driven, or whole clades
#
# sim_bace() can only hide values at random, so the second half is ours. We ask it for NO
# missingness and mask the complete data ourselves, so all three mechanisms go through one code
# path and differ only in which cells they choose.
# =================================================================================================


# ---- the six traits, and the two helper variables ----------------------------------------------
#
#   cat3     3-level category   <- the RESPONSE. It has to be, see note below.
#   x1, x2   continuous         <- the downstream analysis is the PGLS slope of x2 on x1
#   cnt      count
#   bin      binary
#   ord      4-level ordinal
#
#   driver   continuous, ALWAYS observed. Causes the gaps under "mar". Related to x1 at a fixed
#            0.35 in every regime, so that "MAR" really is missing-at-random-given-an-observed-
#            variable, rather than missing at random on a variable unrelated to anything.
#   aux      continuous, ALWAYS observed. The auxiliary -- think body size. Present only in the
#            good_covariate regime, related to x1 at the scenario's aux_beta.
#
# WHY THE CATEGORY IS THE RESPONSE: originally forced, now a free choice we keep. Before BACE
# 1.0.1, sim_bace() dropped all random effects from a multinomial PREDICTOR with no dependencies,
# so a categorical predictor came out with essentially no phylogenetic signal whatever you asked
# for (BACE issue #16, fixed in 1.0.1 -- liability SD across species went 0.000 -> 1.000). The
# constraint no longer binds, but the response path was always the correct one and the generated
# data are identical either way, so there is nothing to gain by changing it.
TRAITS <- c("cat3", "x1", "x2", "cnt", "bin", "ord")
DRIVER_REL <- 0.35


#' Make one dataset: tree, traits, and the complete truth.
#'
#' @param scenario one row of make_scenarios()
#' @param seed     replicate seed
#' @return list(truth, tree, aux, realised) -- no missingness yet, that is the next function
simulate_traits <- function(scenario, seed) {

  ## STEP 1. Who is the response and who are the predictors.
  ##         sim_bace() puts the response first and warns unless the first predictor is gaussian.
  predictors <- c("x1", "x2", "cnt", "bin", "ord", "driver")
  types      <- c("gaussian", "gaussian", "poisson", "binary", "threshold4", "gaussian")
  if (scenario$aux_beta > 0) {
    predictors <- c(predictors, "aux")
    types      <- c(types, "gaussian")
  }
  n_pred <- length(predictors)

  ## STEP 2. How the variables are related. sim_bace() has no correlation matrix -- variables are
  ##         related through regression paths -- so beta_matrix[i, j] is the effect of predictor j
  ##         on predictor i. Everything hangs off x1, which is independent.
  ##         These are REGRESSION COEFFICIENTS, not correlations. The realised correlation is
  ##         measured in step 5 and reported, never assumed.
  B <- matrix(0, n_pred, n_pred, dimnames = list(predictors, predictors))
  for (v in c("x2", "cnt", "bin", "ord")) B[v, "x1"] <- scenario$trait_beta
  B["driver", "x1"] <- DRIVER_REL
  if (scenario$aux_beta > 0) B["aux", "x1"] <- scenario$aux_beta

  ## STEP 3. Phylogenetic signal, one value per variable, response first.
  ##
  ##         This is phylogenetic heritability on the LATENT scale, which is the scale BACE's
  ##         threshold models work on, so it is the scale that matters for imputation. It is
  ##         delivered as specified and equally for every trait type -- measured over 5 datasets
  ##         at n = 300, var(u_phylo) / total latent variance came out:
  ##
  ##            asked 0.3 -> gaussian 0.22 | poisson 0.23 | binary 0.24 | ordinal 0.26
  ##            asked 0.7 -> gaussian 0.60 | poisson 0.61 | binary 0.63 | ordinal 0.66
  ##            asked 0.9 -> gaussian 0.84 | poisson 0.88 | binary 0.87 | ordinal 0.90
  ##
  ##         (The shortfall is sample realisation, not bias: a single Brownian draw has few
  ##         effective degrees of freedom, so its sample variance sits below the population value.)
  ##
  ##         Measured on the OBSERVED values, a binary trait looks far less phylogenetic than a
  ##         continuous one at the same setting. That is what dichotomising a measurement does to
  ##         it, not a flaw in the design -- but it is worth reporting, so step 5 measures it.
  signal <- rep(scenario$phylo_signal, 1 + n_pred)

  ## STEP 3b. Where to cut the ordinal. The hidden value behind `ord` is
  ##
  ##            trait_beta * x1  +  its own species, phylogenetic and residual effects
  ##
  ##         x1 has population variance 2 / (1 - signal) (sim_bace's gaussian defaults: species
  ##         and residual variance 1 each, phylogenetic variance set from the signal), and the
  ##         ordinal's own effects total 0.2 / (1 - signal) (non-gaussian defaults 0.1 + 0.1).
  ##         So the hidden value is Normal with variance (2 * trait_beta^2 + 0.2) / (1 - signal),
  ##         and cutting it at the quartiles of THAT distribution gives four levels of 25 species
  ##         each in expectation at every signal level. These are population constants: they do
  ##         not look at the dataset, so hidden and observed cells are draws from the same
  ##         population and nothing ties one to the other (a sample-quantile cut would make the
  ##         hidden cells' level counts the exact complement of the observed cells', which pushes
  ##         any cutpoint-estimating imputer the wrong way -- measured 2026-11-14). The price is
  ##         that at strong signal a whole tree can shift together: at signal 0.9 and n = 100 about
  ##         1 dataset in 20 has a level with no species in truth, and 1 in 10-20 has a level
  ##         unobserved after masking; at n = 500 none in 10 seeds. Both are recorded per dataset
  ##         below (level counts), so the paper can report them rather than discover them.
  ord_sd  <- sqrt((2 * scenario$trait_beta^2 + 0.2) / (1 - scenario$phylo_signal))
  ord_thr <- stats::qnorm(c(0.25, 0.5, 0.75)) * ord_sd

  ## STEP 4. Simulate. One row per species, and NO missingness -- we add that ourselves next.
  set.seed(seed)
  sim <- suppressMessages(BACE::sim_bace(
    response_type   = "multinomial3",
    predictor_types = types,
    var_names       = c("cat3", predictors),
    beta_matrix     = B,
    beta_resp       = setNames(as.list(rep(0.4, n_pred)), predictors),
    phylo_signal    = signal,
    n_cases         = scenario$n_species,
    n_species       = scenario$n_species,
    missingness     = NULL,
    ordinal_thresholds = list(ord = ord_thr)))

  truth <- sim$complete_data
  rownames(truth) <- truth$Species

  ## STEP 5. Measure what we actually got, rather than trusting the settings. Both of these go into
  ##         the results, because both are weaker than the nominal value asks for and the paper
  ##         should report them.
  realised <- list(
    lambda_x1  = realised_lambda(truth$x1,  sim$tree, truth$Species),
    lambda_cnt = realised_lambda(truth$cnt, sim$tree, truth$Species),
    lambda_bin = realised_lambda(as.integer(truth$bin), sim$tree, truth$Species),
    cor_x1_x2  = stats::cor(truth$x1, truth$x2),
    cor_x1_aux = if (scenario$aux_beta > 0) stats::cor(truth$x1, truth$aux) else NA_real_,
    ## how many species each level of the discrete traits actually has, before masking. An
    ## ordinal level with no species (possible at strong signal, see step 3b) makes that dataset
    ## a 3-level problem; the count is what lets the analysis say so.
    levels_truth = list(ord  = as.integer(table(factor(truth$ord,  levels = 1:4))),
                        cat3 = as.integer(table(factor(truth$cat3, levels = c("A", "B", "C")))),
                        bin  = as.integer(table(factor(truth$bin,  levels = c("0", "1"))))),
    ord_thresholds = ord_thr)

  list(truth = truth, tree = sim$tree, realised = realised)
}


#' Pagel's lambda actually realised by a trait, on the scale it is observed.
#' Returns NA rather than stopping if the fit fails -- one trait's diagnostic is not worth
#' losing a dataset over.
realised_lambda <- function(x, tree, species) {
  v <- tapply(as.numeric(x), species, mean)[tree$tip.label]
  tryCatch(phytools::phylosig(tree, v, method = "lambda")$lambda,
           error = function(e) NA_real_)
}


#' Hide values. The three mechanisms differ only in which cells they pick.
#'
#' `driver` and `aux` are never hidden: one causes the gaps, the other is the fully observed
#' auxiliary. Hiding either would defeat its purpose.
#'
#' @param truth     complete data from simulate_traits()
#' @param tree      the phylogeny, needed by "clade"
#' @param mechanism "mcar", "mar" or "clade"
#' @param frac      target fraction of trait values hidden
#' @param seed      replicate seed (offset so masking does not reuse the simulation's stream)
#' @return list(data = truth with NAs, mask = logical matrix, realised_frac)
apply_missingness <- function(truth, tree, mechanism, frac, seed) {

  ## STEP 1. Which columns may be hidden. Species is an identifier; driver and aux must stay.
  targets <- setdiff(names(truth), c("Species", "driver", "aux"))
  n <- nrow(truth)
  mask <- matrix(FALSE, n, length(targets), dimnames = list(rownames(truth), targets))

  set.seed(seed + 10000L)

  if (mechanism == "mcar") {
    ## Every cell equally likely, each column independently.
    for (v in targets) mask[sample(n, round(frac * n)), v] <- TRUE

  } else if (mechanism == "mar") {
    ## Gaps depend on `driver`, which is fully observed. Species with a high driver value are
    ## likelier to be missing. The intercept is solved so the expected fraction comes out at
    ## `frac` -- otherwise the three mechanisms would not be comparable.
    z <- as.numeric(scale(truth$driver))
    a <- stats::uniroot(function(a) mean(stats::plogis(a + log(3) * z)) - frac, c(-30, 30))$root
    p <- stats::plogis(a + log(3) * z)
    for (v in targets) mask[stats::runif(n) < p, v] <- TRUE

  } else if (mechanism == "clade") {
    ## Whole clades lose everything at once, because whole groups of species go unstudied. This is
    ## the hardest case for phylogenetic imputation -- the close relatives you would borrow from
    ## are exactly the ones that are gone.
    nodes <- (n + 2L):(n + tree$Nnode)
    lo <- max(1L, ceiling(0.05 * n)); hi <- max(lo, ceiling(0.15 * n))
    tries <- 0L
    while (mean(mask) < frac && tries < 2000L) {
      tries <- tries + 1L
      cl <- tryCatch(ape::extract.clade(tree, sample(nodes, 1)), error = function(e) NULL)
      if (is.null(cl) || length(cl$tip.label) < lo || length(cl$tip.label) > hi) next
      mask[rownames(mask) %in% cl$tip.label, ] <- TRUE
    }
    ## The last clade usually overshoots the target; give back the excess at random so the three
    ## mechanisms hide the same amount and differ only in WHERE.
    over <- sum(mask) - round(frac * length(mask))
    if (over > 0) mask[sample(which(mask), over)] <- FALSE

  } else {
    stop("unknown missingness mechanism: ", mechanism)
  }

  ## STEP 2. Punch the holes.
  data <- truth
  for (v in targets) data[mask[, v], v] <- NA

  list(data = data, mask = mask, realised_frac = mean(mask))
}


#' Both halves together: one scenario row + one seed -> one ready-to-impute dataset.
simulate_dataset <- function(scenario, seed) {
  sim <- simulate_traits(scenario, seed)
  miss <- apply_missingness(sim$truth, sim$tree, scenario$mechanism, scenario$miss_frac, seed)
  ## and the level counts BACE will actually see. A level observed 0 times cannot be imputed by
  ## any method and makes MCMCglmm's cutpoint sampler crawl (the hour-long runs of 2026-11-13).
  miss$levels_observed <- list(
    ord  = as.integer(table(factor(miss$data$ord,  levels = 1:4))),
    cat3 = as.integer(table(factor(miss$data$cat3, levels = c("A", "B", "C")))),
    bin  = as.integer(table(factor(miss$data$bin,  levels = c("0", "1")))))
  c(sim, miss)      # truth, tree, realised, data, mask, realised_frac, levels_observed
}
