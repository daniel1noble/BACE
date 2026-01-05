#' @title simBACE: simulator of mock data for BACE
#' @description simBACE is a function to simulate mock data for BACE.
#'  The function simulates data with defined phylogenetic structure and
#'  with multiple dependencies on additional covariates. The data can be
#'  gaussian, binary, poissonian, nomial (unordered) or ordinal.
#' @return Names list with a dataframe of simulated data and underlying tree

simBACE <- function(
    predictor_types,
    response_type = "gaussian",
    var_names = NULL,
    beta_matrix = NULL,
    beta_resp = NULL,
    intercept_p = rep(0, length(predictor_types)),
    intercept_r = 0,
    n_cases = 200,
    n_species = 75,
    missingness = rep(0.25, length(predictor_types) + 1),
    phylo_signal = rep(0.5, length(predictor_types) + 1),
    sigmas = list(
        sigma_s = rep(1, length(predictor_types) + 1),
        sigma_e = rep(1.2, length(predictor_types) + 1)
    ),
    birth = 0.8,
    death = 0.4,
    control_par = NULL) {
    # preparing the variables
    if (is.null(var_names)) var_names <- var_name_gen(length(predictor_types))
    mydesign <- design_size(predictor_types, beta_matrix) # predictor dimensions

    # produce tree and VCV phylo matrix
    taxa <- sim_tree(n_species, birth, death, n_cases)
    cor_phylo <- ape::vcv(taxa$tree, corr = TRUE)

    # define VC variances
    sigma2_s <- sigmas$sigma_s^2
    sigma2_e <- sigmas$sigma_e^2
    sigma2_p <- (phylo_signal * (sigma2_s + sigma2_e)) / (1 - phylo_signal)

    # sample random effects coefficients
    u_s <- lapply(
        sigma2_s,
        function(x) {
            stats::rnorm(length(ape::Ntip(taxa$tree)),
                mean = 0, sd = sqrt(x)
            )
        }
    )
    u_p <- lapply(
        sigma2_p,
        function(x) {
            MASS::mvrnorm(
                n = 1,
                mu = rep(0, length(ape::Ntip(taxa$tree)),
                    Sigma = cor_phylo
                )
            )
        }
    )
    e <- lapply(
        sigma2_e,
        function(x) {
            stats::rnorm(n_cases,
                mean = 0, sd = sqrt(x)
            )
        }
    )

    # form Z matrix for species
    # Z-matrix for species
    Z <- stats::model.matrix(
        ~ 0 + factor(species),
        data.frame(species = taxa$cases)
    )

    # create output predictor set
    covars <- as.data.frame(matrix(0,
        nrow = n_cases,
        ncol = length(predictor_types),
        dimnames = list(NULL, var_names[-1])
    ))

    # simulate the independent predictor
    # must be gaussian and not dependent
    # on other variables (rowSums(beta_matrix)[i]==0)

    covars[, mydesign$sim_order[1]] <- intercept_p[sim_order[1]] +
        Z %*% u_s[[sim_order[1]]] +
        Z %*% u_p[[sim_order[1]]] +
        e[[sim_order[1]]]

    # simulate the remaining predictors

    for (vari in mydesign$sim_order[-1]) {
        if (predictor_types[vari] == "gaussian" ||
            predictor_types[vari] == "binary") {
            # CODE
        } else if (substring(predictor_types[vari], 1, 11) == "categorical") {
            # CODE
        } else {
            stop("Unknown predictor type: ", predictor_types[vari])
        }

        # simulate the response
    }
}
