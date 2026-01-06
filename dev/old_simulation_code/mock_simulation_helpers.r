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
            if (any(x[ix] != 0) && all(x[-ix] == 0)) {
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
    beta_row, ns,
    def_b = 0.5, no_beta = TRUE) {
    vars <- which(beta_row != 0)
    betas <- numeric()

    if (no_beta) {
        betas <- rep(def_b, sum(ns[vars]))
    } else {
        for (i in vars) {
            if (ns[i] == 1) {
                betas <- c(betas, beta_row[i])
            } else {
                betas <- c(
                    betas,
                    seq(
                        from = -1 * (beta_row[i] - 1),
                        to = beta_row[i],
                        length.out = ns[i]
                    )
                )
            }
        }
    }
    return(betas)
}


sim_tree <- function(n_species, birth, death, n_cases, strl = 5) {
    # generate names for species
    mock_names <- sort(stringi::stri_rand_strings(n_species, strl))
    your_cases <- sample(mock_names, n_cases, replace = TRUE)
    your_species <- sort(unique(your_cases))

    # simulate tree
    mytree <- ape::rphylo(
        n = length(your_species),
        birth = birth, death = death, T0 = 100,
        fossils = FALSE
    )
    mytree$tip.label <- your_species

    return(list(tree = mytree, cases = your_cases))
}



var_name_gen <- function(n) {
    var_names <- paste0("x", 1:n)
    var_names <- c("y", var_names)
    return(var_names)
}

beta_extractor <- function(betas) {
    blist <- lapply(betas, function(x) {
        if (length(x) == 1) {
            c(b = x, n = 1)
        } else {
            cbind(b = x, n = seq_along(x))
        }
    })

    blist <- do.call(rbind, blist)
    if (blist[1, "n"] != 1) stop("Check your betas!")
    dim2 <- numeric(nrow(blist))
    dim3 <- numeric(nrow(blist))

    n <- blist[1, "n"]
    dim2[1] <- 1
    dim3[1] <- n
    for (i in 2:nrow(blist)) {
        if (blist[i, "n"] == 1) {
            n <- n + 1
            dim2[i] <- n
            dim3[i] <- blist[i, "n"]
        } else if (blist[i, "n"] > 1) {
            dim2[i] <- n
            dim3[i] <- blist[i, "n"]
        } else {
            stop("Beta coefs do not make sense.")
        }
    }

    blist <- cbind(blist, dim2, dim3)
    return(blist)
}

array_3_fill <- function(myarray, indexes, values) {
    for (i in seq_along(values)) {
        myarray[
            indexes[i, 1],
            indexes[i, 2],
            indexes[i, 3]
        ] <- values[i]
    }

    return(myarray)
}


beta_constructor <- function(mydesign, beta_coefs) {
    ns <- mydesign$ns

    if (!(all(sapply(beta_coefs, length) == length(ns)))) {
        stop("Ill-formed beta definition list!")
    }

    beta_array <- array(dim = c(
        length(ns),
        length(ns),
        max(ns)^2
    ))

    for (i in seq_along(ns)) {
        betas <- beta_coefs[[i]]
        indexes <- beta_extractor(betas)

        beta_array <- array_3_fill(
            beta_array,
            cbind(i, indexes[, "dim2"], indexes[, "dim3"]),
            indexes[, "b"]
        )
    }

    dimnames(
        beta_array,
        d
    )
    return(beta_array)
}

