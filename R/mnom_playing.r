Delta <- cbind(c(-1, 1, 0), c(-1, 0, 1))
c2 <- (16 * sqrt(3) / (15 * pi))^2
D <- ginv(Delta %*% t(Delta)) %*% Delta

IJ <- (1/3) * (diag(2) + matrix(1, 2, 2))

library(MASS)

y <- mvrnorm(1000, mu = c(0, 0), Sigma = IJ)

plot(y[, 1], y[, 2])

Int <- t(apply(y, 1, function(x) {
D %*% (x/sqrt(1 + c2 * diag(IJ)))
}))

pr <- exp(Int) / rowSums(exp(Int))

plot3d(pr)
TernaryPlot()
# TernaryPoints(pr)
TernaryDensityContour(pr, resolution = 50, nlevels = 8)
