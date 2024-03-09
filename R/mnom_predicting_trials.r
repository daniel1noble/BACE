library(MCMCglmm)
library(MASS)


 <- rnorm(500, 0, 1.5)

l1 <- 0.2 * x - 1
l2 <- 1 - 2 * x

toy <- cbind(l1, l2, l3 = rowMeans(cbind(l1, l2)))
toy <- plogis(toy)
toy <- toy / rowSums(toy)
toy

toy2 <- data.frame(
    x = x,
    y = factor(apply(toy, 1, function(x) sample(c("A", "B", "C"), 1, prob = x)))
)

Delta <- cbind(c(-1, 1, 0), c(-1, 0, 1))
c2 <- (16 * sqrt(3) / (15 * pi))^2
D <- ginv(Delta %*% t(Delta)) %*% Delta

IJ <- (1/3) * (diag(2) + matrix(1, 2, 2))
prior = list(R = list(V = IJ, fix = 1))

model <- MCMCglmm(y ~ trait - 1,
    rcov = ~ us(trait):units,
    prior = prior, data = toy2, family = "categorical",
    verbose = FALSE, nitt = 300000, thin = 150, burnin = 30000,
    pr = T, pl = T, saveX = T, saveZ = T
)
summary(model)

liab <- model$Liab
X <- model$X
sol <- model$Sol

dim(liab)

Int <- t(apply(model$Sol, 1, function(x) {
D %*% (x/sqrt(1 + c2 * diag(IJ)))
}))

apply(as.matrix(model$Sol), 1, function(x) as.matrix(model$X) %*% x) -> pred
predm <- cbind(rowMeans(pred)[1:500], rowMeans(pred)[501:1000])


Int <- t(apply(model$Sol, 1, function(x) {
D %*% (x/sqrt(1 + c2 * diag(IJ)))
}))
summary(mcmc(exp(Int)/rowSums(exp(Int))))



t(apply(predm, 1, function(x) {
    D %*% (x / sqrt(1 + c2 * diag(IJ)))
})) -> Int2

exp(Int2) / rowSums(exp(Int2))

