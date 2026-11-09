# =================================================================================================
# STEP 1 of 4 -- the list of scenarios.
#
# A "scenario" is one combination of simulation settings. This file turns the design in
# .agents/plans/simulation-plan-2026-11.md into a plain data frame, one row per scenario.
# Everything downstream reads a row of this table and nothing else, so if you want to know what
# the campaign will run, print this table.
# =================================================================================================


#' Build the scenario grid.
#'
#' @param n_species    tree sizes to run
#' @param phylo_signal phylogenetic signal levels (Pagel's lambda, roughly)
#' @param mechanism    how the gaps are made
#' @param miss_frac    fraction of values hidden in each trait
#' @return a data frame, one row per scenario, with an `id` you can use as a filename
make_scenarios <- function(n_species    = c(100, 500),
                           phylo_signal = c(0.3, 0.7, 0.9),
                           mechanism    = c("mcar", "mar", "clade"),
                           miss_frac    = 0.30) {

  ## The three INFORMATION REGIMES. This is ONE axis with three levels, not two crossed axes.
  ## Each says how much a missing trait can borrow from somewhere other than the phylogeny:
  ##
  ##   floor           traits unrelated, no auxiliary    -- nothing to borrow but the tree
  ##   typical         traits related, no auxiliary      -- the ordinary situation
  ##   good_covariate  traits related, PLUS one fully    -- you have a body-size-like variable
  ##                   observed variable related 0.8        measured for every species
  ##
  ## Claim 2 of the manuscript is the difference between the last two at low signal.
  regimes <- data.frame(
    regime     = c("floor", "typical", "good_covariate"),
    trait_rel  = c(0,       0.5,       0.5),   # how related the traits are to each other
    aux_rel    = c(0,       0,         0.8)    # how related the auxiliary is to its target; 0 = none
  )

  ## Every combination of the four axes.
  grid <- expand.grid(n_species    = n_species,
                      phylo_signal = phylo_signal,
                      regime       = regimes$regime,
                      mechanism    = mechanism,
                      stringsAsFactors = FALSE)

  ## Attach each regime's two numbers.
  grid <- merge(grid, regimes, by = "regime")

  grid$miss_frac <- miss_frac

  ## A short, unique, filename-safe name for each scenario. EVERY setting that changes the result
  ## is in it: if two scenarios shared a name, the runner would skip the second as "already done"
  ## and it would never run.
  grid$id <- paste0("n",   grid$n_species,
                    "_ps", grid$phylo_signal,
                    "_",   grid$regime,
                    "_",   grid$mechanism)

  ## Put the columns in a readable order and sort so the table reads top to bottom.
  grid <- grid[order(grid$n_species, grid$phylo_signal, grid$regime, grid$mechanism),
               c("id", "n_species", "phylo_signal", "regime", "trait_rel", "aux_rel",
                 "mechanism", "miss_frac")]
  rownames(grid) <- NULL

  stopifnot(!anyDuplicated(grid$id))
  grid
}


#' Expand the scenario grid into one row per DATASET: every scenario x every replicate seed.
#' This is the full list of jobs the campaign has to run.
#'
#' @param scenarios a make_scenarios() table
#' @param n_reps    replicate datasets per scenario
expand_jobs <- function(scenarios, n_reps = 1000) {
  jobs <- scenarios[rep(seq_len(nrow(scenarios)), each = n_reps), ]
  jobs$seed <- rep(seq_len(n_reps), times = nrow(scenarios))
  jobs$job  <- paste0(jobs$id, "_s", jobs$seed)
  rownames(jobs) <- NULL
  stopifnot(!anyDuplicated(jobs$job))
  jobs
}
