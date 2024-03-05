#' @title Beta generator
#' @description This code generates default regression coefficients
#'   based on the supplied design properties
#' @param ns A vector of numbers of beta coefs needed from
#'   design_size function
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
