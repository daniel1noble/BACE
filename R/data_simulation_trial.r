#' @importFrom stringi stri_rand_strings

response_type <- "gaussian"
predictor_types <- c(
    "gaussian", "gaussian",
    "binary", "categorical"
)

n_cases <- 200
n_species <- 75
missingness <- c(0.5, 0.1, 0.2, 0, 0.1)

mock_names <- stringi::stri_rand_strings(n_species, 5)


