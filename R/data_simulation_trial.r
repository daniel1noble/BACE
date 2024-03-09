#' @importFrom stringi stri_rand_strings

# functions and helpers ------------------------------------------
source("./R/design_size.r")
source("./R/beta_generator.r")
source("./R/mnom_liab2cat.r")

# fine-tuning parameters ------------------------------------------
control_par <- list(
    mnomial_cor = 0.25
)

# entry conditions -----------------------------------------------
response_type <- "gaussian"
predictor_types <- c(
    "gaussian", "gaussian",
    "binary", "categorical3"
)
var_names <- c("resp", "x1", "x2", "x3", "x4")

dep_matrix <- rbind(
    c(0, 1, 0, 1),
    c(0, 0, 0, 0),
    c(1, 0, 0, 0),
    c(0, 1, 0, 0)
)

beta_matrix <- rbind(
    c(0, 0.5, 0, 4),
    c(0, 0, 0, 0),
    c(0.1, 0, 0, 0),
    c(0, 0.15, 0, 0)
)

beta_resp <- c(1, 2, 0.5, 0.25)

intercept_p <- c(7, 2, 0, -0.05)
intercept_r <- 10
resp_betas <- c(0.5, 0.25, 0.75, 0.1, -0.1)

phylo_signal_resp <- 0.8
phylo_signal_vars <- c(0.5, 0.3, 0.2, 0.1)
phylo_signal <- c(phylo_signal_resp, phylo_signal_vars)

n_cases <- 200
n_species <- 75
missingness <- c(0.5, 0.1, 0.2, 0, 0.1)

mydesign <- design_size(predictor_types, dep_matrix)

mock_names <- sort(stringi::stri_rand_strings(n_species, 5))

# parameters for tree diversification
birth <- 0.8
death <- 0.4

# simulate the tree -----------------------------------------------
your_cases <- sample(mock_names, n_cases, replace = TRUE)
your_species <- sort(unique(your_cases))

mytree <- ape::rphylo(
    n = length(your_species),
    birth = birth, death = death, T0 = 100,
    fossils = FALSE
)
mytree$tip.label <- your_species
plot(mytree, cex = 0.5)
cor_phylo <- ape::vcv(mytree, corr = TRUE)



# sigmas for variance components ---------------------------------

sigma2_s <- c(1^2, 1^2, 1^2, 1^2, 1^2)
sigma2_e <- c(1.2^2, 1.2^2, 1.2^2, 1.2^2, 1.2^2)
sigma2_p <- (phylo_signal * (sigma2_s + sigma2_e)) / (1 - phylo_signal)

# ranom effects ---------------------------------------------------

# ranom effect coefs
u_s <- lapply(
    sigma2_s,
    function(x) {
        stats::rnorm(length(your_species),
            mean = 0, sd = sqrt(x)
        )
    }
)

u_p <- lapply(
    sigma2_p,
    function(x) {
        MASS::mvrnorm(
            n = 1,
            mu = rep(0, length(your_species)),
            Sigma = x * cor_phylo
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

# Z-matrix for species
Z <- stats::model.matrix(
    ~ 0 + factor(species),
    data.frame(species = your_cases)
)

# simulate the predictors ----------------------------------------
covars <- as.data.frame(matrix(0,
    nrow = n_cases,
    ncol = length(predictor_types),
    dimnames = list(NULL, var_names[-1])
))

# liab_ix <- unlist(mapply(
#     rep.int,
#     seq_along(mydesign$ns),
#     mydesign$ns
# ))

covars[, mydesign$sim_order[1]] <- intercept_p[1] +
    Z %*% u_s[[1]] +
    Z %*% u_p[[1]] +
    e[[1]]



for (vari in mydesign$sim_order[-1]) {
    if (predictor_types[vari] == "gaussian" ||
        predictor_types[vari] == "binary") {
        beta0 <- intercept_p[vari]

        which_x <- which(dep_matrix[vari, ] == 1)
        xi_formula <- paste(
            "~1",
            paste(names(covars)[which_x], collapse = " + "),
            sep = " + "
        )

        X <- model.matrix(
            as.formula(xi_formula),
            covars[, which_x, drop = FALSE]
        )

        beta <- c(beta0, beta_generator(
            beta_matrix = beta_matrix,
            ns = mydesign$ns,
            vari = vari,
            no_beta = FALSE
        ))

        liab_i <- X %*% beta +
            Z %*% u_s[[vari]] +
            Z %*% u_p[[vari]] +
            e[[vari]]

        if (predictor_types[vari] == "binary") {
            resp_i <- rbinom(n_cases, 1, plogis(liab_i))
        } else {
            resp_i <- liab_i
        }
    } else if (substring(predictor_types[vari], 1, 11) == "categorical") {
        beta0 <- seq(-1 * (intercept_p[vari] - 1),
            intercept_p[vari],
            length.out = mydesign$ns[vari]
        )

        which_x <- which(dep_matrix[vari, ] == 1)

        xi_formula <- paste(
            "~trait-1",
            paste(names(covars)[which_x], collapse = " + "),
            sep = " + "
        )
        covars_mv <- rbind(covars, covars)
        covars_mv$trait <- factor(rep(paste0("L", 1:mydesign$ns[vari]),
            each = n_cases
        ))

        X <- model.matrix(
            as.formula(xi_formula),
            covars_mv[, c(which_x, ncol(covars_mv))]
        )


        beta <- c(beta0, beta_generator(
            beta_matrix = beta_matrix,
            ns = mydesign$ns,
            vari = vari,
            no_beta = FALSE
        ))

        liab_i <- X %*% beta +
            rbind(Z, Z) %*% u_s[[vari]] +
            rbind(Z, Z) %*% u_p[[vari]] +
            sample(e[[vari]], n_cases, replace = TRUE)
        resp_i <- mnom_liab2cat(
            liab_i, n_cases,
            categories = LETTERS[1:(mydesign$ns[vari] + 1)]
        )
    } else {
        stop("Unknown predictor type: ", predictor_types[vari])
    }

    covars[, vari] <- resp_i
}

# simulate the response -------------------------------------------
if (response_type == "gaussian" || response_type == "binary") {
    beta0 <- intercept_r
    X <- model.matrix(
        ~ 1 + x1 + x2 + x3 + x4,
        covars
    )
    beta <- c(beta0, resp_betas)
    liab_r <- X %*% beta +
        Z %*% u_s[[5]] +
        Z %*% u_p[[5]] +
        e[[5]]


    if (response_type == "binary") {
        resp <- rbinom(n_cases, 1, plogis(liab_r))
    } else {
        resp <- liab_r
    }
} else if (substring(response_type, 1, 11) == "categorical") {
    beta0 <- seq(-1 * (intercept_r - 1),
        intercept_r,
        length.out = as.numeric(substring(response_type, 12, 12))
    )

    covars_mv <- rbind(covars, covars)
    covars_mv$trait <- factor(rep(paste0("L", 1:mydesign$ns[vari]),
        each = n_cases
    ))

    X <- model.matrix(
        ~ trait - 1 + x1 + x2 + x3 + x4,
        covars_mv
    )

    beta <- c(beta0, resp_betas)
    liab_r <- X %*% beta +
        rbind(Z, Z) %*% u_s[[5]] +
        rbind(Z, Z) %*% u_p[[5]] +
        sample(e[[5]], n_cases, replace = TRUE)
    resp <- mnom_liab2cat(
        liab_r, n_cases,
        categories = LETTERS[1:as.numeric(substring(response_type, 12, 12))]
    )
} else {
    stop("Unknown response type: ", response_type)
}

outdata <- cbind(spec = your_cases, resp = resp, covars)
ape::write.nexus(mytree, file = "tree.nex")
write.table(outdata, file = "data.csv", sep = ",", row.names = FALSE)

simGaussBernPois(
    covars = covars, vari = 3, beta_matrix = c(0.1,0.2,-0.1,0.1), beta0 = 0,
    mydesign = mydesign, u_s = u_s, u_p = u_p, e = e, Z = Z,
    family = "poisson", formula = NULL, Liab = F, event = T
)



beta_coefs <- list(
    x1 = list(0, 0.5, 0, c(-1, 2)),
    x2 = list(0, 0, 0, 0),
    x3 = list(0.1, 0, 0, 0),
    x4 = list(0, c(0.15, 0.33), 0, 0)
)
