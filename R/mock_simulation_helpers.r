#' @title Multinomial liability to category
#' @description Used to convert (multivariate) liability scale odds-ratios
#'  to a categorical scale
#' @param liability A matrix of liabilities (N-1 columns when N categories)
#' @param n_cases The number of cases (rows of data)
#' @param categories A vector of category names to label final response levels
#' @return A character vector of category names

mnom_liab2cat <- function(liability, n_cases, categories) {
    # turn stcked liabilities to N-1 columns
    liab_wide <- matrix(liability, nrow = n_cases)

    # I = identity matrix, J = unit matrix
    # IJ defines the final N-simplex space used to constrain the probabilities
    IJ <- (1 / 3) * (diag(2) + matrix(1, 2, 2))
    c2 <- (16 * sqrt(3) / (15 * pi))^2

    # form contrast matrix for data categories
    delta <- diag(length(categories) - 1)
    delta <- rbind(rep(-1, length(categories) - 1), delta)

    # helper matrix to project liabilities to probabilities
    Delta <- MASS::ginv(delta %*% t(delta)) %*% delta

    # project liabilities to (constrained) log-odds
    # assuming zero residual variance
    projection <- t(apply(liab_wide, 1, function(x) {
        Delta %*% (x / sqrt(1 + c2 * diag(IJ)))
    }))

    # convert log-odds to probabilities
    prob <- exp(projection) / rowSums(exp(projection))

    return(
        apply(prob, 1, function(x) {
            sample(categories, 1, prob = x)
        })
    )
}


#' @title Design size & order
#' @description Function to calculate the design size and order of the
#'   predictors based on the predictor types and the dependence matrix
#' @param predictor_types A vector of strings indicating the type of each
#'   predictor
#' @param dep_matrix A matrix indicating the dependence between predictors
#' @return A list with the design size and the order of the predictors

design_size <- function(predictor_types, dep_matrix) {
    if (!(nrow(dep_matrix) == ncol(dep_matrix) &&
        nrow(dep_matrix) == length(predictor_types))) {
        stop("Dependence matrix and predictor types do not match")
    }

    ns <- numeric(length(predictor_types))

    for (i in seq_along(predictor_types)) {
        type <- predictor_types[i]
        if (type == "gaussian") {
            n <- 1
        } else if (type == "binary") {
            n <- 1
        } else if (substring(type, 1, 11) == "categorical") {
            n <- as.numeric(substring(type, 12, 12)) - 1
        } else {
            stop("Unknown predictor type: ", type)
        }
        ns[i] <- n
    }

    if (any(rowSums(dep_matrix) == 0)) {
        sim_order <- which(rowSums(dep_matrix) == 0)
    } else {
        stop("Min one predictor must be independent")
    }

    while (!all(seq_along(predictor_types) %in% sim_order)) {
        # print(sim_order)
        next_ix <- apply(dep_matrix, 1, function(x, ix) {
            if (any(x[ix] == 1) && all(x[-ix] == 0)) {
                return(1)
            } else {
                return(0)
            }
        }, ix = sim_order)
        next_ix[sim_order] <- 0
        sim_order <- c(sim_order, which(as.logical(next_ix)))
    }

    if (predictor_types[sim_order[1]] != "gaussian") {
        stop("First independent predictor must be gaussian")
    }

    return(list(ns = ns, sim_order = sim_order))
}


#' @title Beta generator
#' @description This code generates regression coefficients
#'   in correct format and in correct numbers
#'   based on the supplied design properties
#' @param ns A vector of numbers of beta coefs needed from
#'   design_size function (depends on the type of variable)
#' @param beta_matrix A matrix of beta coefficients conforming to
#'  the dependence matrix (reg coefs instead of 1's)
#' @param vari The index of the predictor to generate betas for
#' @param def_b The default beta coefficient to use if no beta_matrix
#' @param no_beta A boolean indicating whether to use beta_matrix or
#'  just the default beta
#' @return A list of beta coefs

beta_generator <- function(
    beta_matrix, ns, vari,
    def_b = 0.5, no_beta = TRUE) {
    vars <- which(beta_matrix[vari, ] != 0)
    betas <- numeric()

    if (no_beta) {
        betas <- rep(def_b, sum(ns[vars]))
    } else {
        raw_betas <- beta_matrix[vari, ]
        for (i in vars) {
            if (ns[i] == 1) {
                betas <- c(betas, raw_betas[i])
            } else {
                betas <- c(
                    betas,
                    seq(
                        from = -1 * (raw_betas[i] - 1),
                        to = raw_betas[i],
                        length.out = ns[i]
                    )
                )
            }
        }
    }
    return(betas)
}
