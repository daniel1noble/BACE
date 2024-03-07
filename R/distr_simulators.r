#' @title Data sampling function for gaussian, binary, Poisson
#' @description This function samples data from a given distribution
#'  for a defined linear predictor
#' @return A data vector on the response scale or liability scale
#' @param covars A dataframe with the covariate columns
#' @param beta_matrix A matrix of regression coefficients (mice-style)
#' @param beta0 A numeric value for the intercept
#' @param mydesign A list with the design size and the order of the predictors
#' @param u_s A vector of random effect coefs for species
#' @param u_p A vector of random effect coefs for phylogeny
#' @param e A vector of random residuals
#' @param Z A design model matrix for species/phylogeny
#' @param family A string indicating the distribution of the focal variable
#' @param formula A formula for the model (not implemented)
#' @param Liab A logical indicating whether to return the liability scale
#' @param event A logical indicating whether to return the response event scale
#' @param vari A numeric indicating the index of the focal variable
#' @export

simGaussBernPois <- function(
    covars,
    beta_matrix,
    beta0,
    mydesign,
    u_s, u_p, e, Z,
    vari,
    family = "gaussian",
    formula = NULL,
    Liab = FALSE, event = FALSE) {

    which_x <- which(beta_matrix[vari, ] != 0)

    # generate model formula
    xi_formula <- paste(
        "~1",
        paste(names(covars)[which_x], collapse = " + "),
        sep = " + "
    )

    # generate model matrix
    X <- stats::model.matrix(
        stats::as.formula(xi_formula),
        covars[, which_x, drop = FALSE]
    )

    # generate beta coefficients
    if (is.matrix(beta_matrix)) {
        beta <- c(beta0, beta_generator(
            beta_matrix = beta_matrix,
            ns = mydesign$ns,
            vari = vari,
            no_beta = FALSE
        ))
    } else {
        beta <- c(beta0, beta_matrix)
    }

    cat("Variable mode: ", family, "\n")
    cat("Generating formula: ", xi_formula, "\n")
    cat("Beta coefficients: ", beta, "\n")

    # generate liability
    liab_i <- X %*% beta +
        Z %*% u_s[[vari]] +
        Z %*% u_p[[vari]] +
        e[[vari]]
    if (Liab || family == "gaussian") {
        return(liab_i)
    } else if (family == "binary") {
        resp_i <- stats::plogis(liab_i)
        if (event) {
            return(stats::rbinom(n_cases, 1, resp_i))
        } else {
            return(resp_i)
        }
    } else if (family == "poisson") {
        resp_i <- exp(liab_i)
        if (event) {
            return(stats::rpois(n_cases, resp_i))
        } else {
            return(resp_i)
        }
    } else {
        stop("Family not recognized")
    }
}


#' @title Data sampling function for categorical variables (unordered)
#' @description This function samples data from a multinomial logit process
#' @return A data vector on the response scale or (matrix) on liability scale
#' @param covars A dataframe with the covariate columns
#' @param beta_matrix A matrix of regression coefficients (mice-style)
#' @param beta0 A numeric value for the intercept
#' @param mydesign A list with the design size and the order of the predictors
#' @param u_s A vector of random effect coefs for species
#' @param u_p A vector of random effect coefs for phylogeny
#' @param e A vector of random residuals
#' @param Z A design model matrix for species/phylogeny
#' @param vari A numeric indicating the index of the focal variable
#' @param catlabels A vector of labels for the categories
#' @param family A string indicating the distribution of the focal variable
#' @param formula A formula for the model (not implemented)
#' @param Liab A logical indicating whether to return the liability scale
#' @param event A logical indicating whether to return the response event scale
#' @export


simMnomial <- function(
    covars,
    beta_matrix,
    beta0, # 2-element vector
    mydesign,
    u_s, u_p, e, Z,
    vari,
    family = "categorical",
    formula = NULL,
    Liab = FALSE, event = FALSE,
    catlabels = LETTERS[1:(mydesign$ns[vari] + 1)]) {

    which_x <- which(beta_matrix[vari, ] != 0)

    if (family != "categorical") {
        stop("Family not recognized")
    }

    # generate model formula
    xi_formula <- paste(
        "~trait-1",
        paste(names(covars)[which_x], collapse = " + "),
        sep = " + "
    )

    covars_mv <- rbind(covars, covars)
    covars_mv$trait <- factor(rep(paste0("L", 1:mydesign$ns[vari]),
        each = nrow(covars)
    ))

    # generate model matrix
    X <- stats::model.matrix(
        stats::as.formula(xi_formula),
        covars_mv[, c(which_x, ncol(covars_mv))]
    )

    # generate beta coefficients
    beta <- c(beta0, beta_generator(
        beta_matrix = beta_matrix,
        ns = mydesign$ns,
        vari = vari,
        no_beta = FALSE
    ))

    # generate liability
    liab_i <- X %*% beta +
        rbind(Z, Z) %*% u_s[[vari]] +
        rbind(Z, Z) %*% u_p[[vari]] +
        sample(e[[vari]], nrow(covars_mv), replace = TRUE)

    if (Liab) {
        return(liab_i)
    } else if (event) {
        resp_i <- mnom_liab2cat(
            liab_i, nrow(covars),
            categories = catlabels
        )
    } else {
        cat("Warning: even probabilities not scaled to unity sum!")
        return(stats::plogis(liab_i))
    }
}