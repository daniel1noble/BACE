#' @importFrom stringi stri_rand_strings

# functions and helpers ------------------------------------------
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

design_size(predictor_types, dep_matrix)

# entry conditions -----------------------------------------------
response_type <- "gaussian"
predictor_types <- c(
    "gaussian", "gaussian",
    "binary", "categorical3"
)
var_names("resp", "x1", "x2", "x3", "x4")

dep_matrix <- rbind(
    c(0, 1, 0, 1),
    c(0, 0, 0, 0),
    c(1, 0, 0, 0),
    c(0, 1, 0, 0)
)

phylo_signal_resp <- 0.8
phylo_signal_vars <- c(0.5, 0.3, 0.2, 0.1)

n_cases <- 200
n_species <- 75
missingness <- c(0.5, 0.1, 0.2, 0, 0.1)

mock_names <- stringi::stri_rand_strings(n_species, 5)

# parameters for tree diversification
birth <- 0.8
death <- 0.4

# simulate the tree -----------------------------------------------
mytree <- ape::rphylo(
    n = n_species, birth = birth, death = death, T0 = 100,
    fossils = FALSE
)
mytree$tip.label <- mock_names
plot(mytree, cex = 0.5)
cor_phylo <- ape::vcv(mytree, corr = TRUE)

your_species <- sample(mock_names, n_cases, replace = TRUE)
