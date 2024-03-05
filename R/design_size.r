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

    return(list(ns = ns, sim_order = sim_order))
}