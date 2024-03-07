

#Sims

y <- sample(c("red", "yel", "bla"), 10000, prob = c(0.5, 0.3, 0.2),replace = TRUE)
x <- rnorm(10000),
clus <- as.factor(sample(1:10, 10000, replace = TRUE))

data <- data.frame(y = as.factor(y), x = x, clus = clus)

#MCMCglmm

prior <-  prior <- list(R = list(V = diag(2), fix = 1),
                        G = list(G1 = list(V = diag(4), nu =4)))

mod <- MCMCglmm::MCMCglmm(y ~ -1+trait + trait:x, random = ~us(-1+trait + trait:x):clus, family = "categorical", nitt = 100000, pl = TRUE, pr = TRUE, thin = 5, prior = prior, burnin = 10000, verbose = FALSE, rcov = ~us(trait):units, data = data)

exp(0.5) - exp(0.2)
exp(0.3) - exp(0.2)

table(y) / sum(table(y))
exp(0.9482) / (1 + exp(0.9482))
p_red <- exp(0.9482) / (1+sum(c(exp(0.9482), exp(0.4296))))
p_yel <- exp(0.4296) / (1+sum(c(exp(0.9482), exp(0.4296))))
p_bla <- 1 - p_red - p_yel

X <- mod$X
pred <- mod$Liab
dim(pred)
means_pred <- apply(pred, 2, mean)
red  <- mean(means_pred[1:10000])
yel  <- mean(means_pred[10001:20000])
exp(1.006459) / (1+sum(c(exp(1.006459), exp(0.2766775))))

exp(0.5333103) / (1+sum(c(exp(0.5333103), exp(0.64866736))))
exp(0.64866736) / (1+sum(c(exp(0.5333103), exp(0.64866736))))


# begin imputation we select the most likely



# end imputation
imp <- sample(levels(y), 1, prob = c(p_red, p_yel, p_bla))


dat <- data.frame(num = 1:30, cat = sample(c("red", "yel", "bla"), 30, prob = c(0.5, 0.3, 0.2), replace = TRUE), x = rnorm(30), trait = sample(c("red", "yel", "bla"), 30, prob = c(0.5, 0.3, 0.2), replace = TRUE))

X <- model.matrix(~cat + x, data = dat)
X2 <- model.matrix(~num + x + trait, data = dat)

x <- rnorm(1000)
y <- x*0.5 + rnorm(1000)

data <- data.frame(y = y, x = x)

mod2 <- MCMCglmm::MCMCglmm(y ~ x,  family = "gaussian", nitt = 100000, pl = TRUE, pr = TRUE, thin = 5,  burnin = 10000, verbose = FALSE,  data = data)

     means <- apply(mod2$Liab, 2, mean)
means_pred <- predict(mod2, marginal = NULL)


length(means)
dim(means_pred)

x <- rnorm(1000)
y <- rpois(1000, exp(x*0.5))

data <- data.frame(y = y, x = x)

mod3 <- MCMCglmm::MCMCglmm(y ~ x,  family = "poisson", nitt = 100000, pl = TRUE, pr = TRUE, thin = 5,  burnin = 10000, verbose = FALSE,  data = data)

     means <- apply(mod3$Liab, 2, mean)
	 means <- apply(exp(mod3$Liab + 0.5*(mod3$VCV[,1])), 2, function(x) mean(x))
means_pred <- predict(mod3)
head(means)
head(means_pred)
