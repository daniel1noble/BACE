mnom_liab2cat <- function (liability, n_cases, categories) {
    liab_wide <- matrix(liability, nrow = n_cases)
    IJ <- (1 / 3) * (diag(2) + matrix(1, 2, 2))
    c2 <- (16 * sqrt(3) / (15 * pi))^2

    delta <- diag(length(categories) - 1)
    delta <- rbind(rep(-1, length(categories) - 1), delta)
    Delta <- MASS::ginv(delta %*% t(delta)) %*% delta

    projection <- t(apply(liab_wide, 1, function(x) {
        Delta %*% (x / sqrt(1 + c2 * diag(IJ)))
    }))

    prob <- exp(projection) / rowSums(exp(projection))

    return(
        apply(prob, 1, function(x) {
            sample(categories, 1, prob = x)
        })
    )
}
