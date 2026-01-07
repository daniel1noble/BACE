data <- data.frame(
    yc <- sample(c("A", "B", "C"), size = 1000, replace = TRUE),
    y = rbinom(1000, size = 1, prob = 0.6)
)
head(data)

library(MCMCglmm)

model <- MCMCglmm(
    fixed = yc ~ 1,
    family = "categorical",
    data = data,
    verbose = FALSE,
    nitt = 100000,
    burnin = 40000,
    thin = 50,
    pl = TRUE,
    prior = list(
        R = list(V = diag(2), fix = 2, nu = 2.002)
    ),
    rcov = ~idh(trait):units
)

c <- 16 * sqrt(3) / (15 * pi)

summary(model)
mean(model$Liab)
mean(plogis(model$Liab))
mean(plogis(model$Liab/sqrt(1 + c^2 * model$VCV[, 1])))

mean(model$Sol[, 1] / sqrt(1 + c^2 * model$VCV[, 1]))
qlogis(0.6)

predict(model, type = "response", marginal = NULL)





## testing for liability correction

# simulate simple data from multinomial distribution - intercept only
# categories probabilities 0.2, 0.7, 0.1
set.seed(123)
n <- 800
p <- c(0.2, 0.7, 0.1)
y <- sample(1:3, size = n, replace = TRUE, prob = p)
data <- data.frame(y = factor(y))

# using J hadfield prior parametrisation
IJ <- 1 / 3 * (diag(2) + matrix(1, nrow = 2, ncol = 2))

library(MCMCglmm)
library(MASS)
model <- MCMCglmm(
    fixed = y ~ trait - 1,
    family = "categorical",
    data = data,
    verbose = FALSE,
    nitt = 140000,
    burnin = 40000,
    thin = 100,
    pl = TRUE,
    prior = list(
        R = list(V = IJ, fix = 1)
    ),
    rcov = ~us(trait):units
)

summary(model)
# recover original probabilities from the model
dim(model$Liab)


# calculating predicted probabilities manually
Delta <- cbind(c(-1, 1, 0), c(-1,0,1))
c2 <- ((16 * sqrt(3)) / (15 * pi))^2
corr <- c2
D <- ginv(Delta %*% t(Delta)) %*% Delta
Int <- t(apply(model$Sol, 1, function(x) D %*% (x / sqrt(1 + corr * diag(IJ)))))
summary(mcmc(exp(Int) / rowSums(exp(Int)))) # match with sim values 0.2 0.7 0.1


# using the Liab object instead
L <- model$Liab
n_traits <- 2
n_obs <- n
corr <- c2

for (i in 1:n_traits) {
    L[, ((i - 1) * (n_obs + 1) + 1):(i * (n_obs))] <- 
        L[, ((i - 1) * (n_obs + 1) + 1):(i * (n_obs))] /
        sqrt(1 + corr * IJ[i, i])    
}

L_mean <- NULL

for (i in 1:n_traits) {
    L_mean <- cbind(L_mean,
        rowMeans(L[, ((i - 1) * (n_obs + 1) + 1):(i * (n_obs))])
    )
}
P_mean <- exp(L_mean) / (rowSums(exp(L_mean))+1)
P_mean <- cbind(1/(rowSums(exp(L_mean))+1), P_mean)
summary(mcmc(P_mean))
