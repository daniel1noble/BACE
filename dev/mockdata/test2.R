# test for caterogical data

library(here)
library(tidyverse)
library(ape)
library(MCMCglmm)
library(MASS)


dat <- read.csv(here("R","mockdata","data.csv"))

tree <- read.nexus(here("R","mockdata","mock_data_tree.nex"))

#dat$x4 <- factor(dat$x4)
#dat$animal <- dat$spec

Ainv <- inverseA(tree)$Ainv

# number of y

n <- dim(dat)[1]

#number of categories
J <- length(unique(dat$x4))
K <- J-1 # the number of contrast
k <- 10000 # the number of iteration

number_fix_parameters <- ncol(dat) - 2 + (J-1)
number_random_effects <- 1
number_random_parameters <- 1 + (J - 1)

J_matrix <- array(1, dim = c(J, J) - 1) # matrix of ones
I_matrix <- diag(J - 1) #identiy matrix


# running a categorical model
IJ <- (I_matrix + J_matrix)/J # see Hadfield's Course notes p. 97
prior <- list(R = list(V = IJ, fix = 1),
              G = list(G1 = list(V = diag(number_random_parameters), nu = J+number_random_effects)))
prior2 <- list(R = list(V = IJ, fix = 1),
              G = list(G1 = list(V = diag(2), nu = 0.002)))

model <- MCMCglmm::MCMCglmm(x4 ~ -1 + trait, random = ~us(-1+trait):spec, family = "categorical",
                          nitt = 110000, pl = TRUE, pr = TRUE, thin = 10, prior = prior2, burnin = 10000,
                          ginverse=list(spec = Ainv),, verbose = FALSE, rcov = ~us(trait):units, data = dat)

summary(model)

# getting effects out
#> dim(model$Liab)
#[1] 10000   400
#test <- model$Liab[, c(1, 201)]

# look up list for levels in a category
cat_list <- levels(unique(factor(dat$x4)))

# n is the number of y
# K the number of category

# posterior

posterior.mean <- summary(model$Sol)$statistics[,1]
# get posterior mean of random effects
zpost.mean.extra <- posterior.mean[(ncol(model$X) + 1):length(posterior.mean)]
# remove the .Node columns
zpost.mean <- zpost.mean.extra[!grepl("\\.Node", names(zpost.mean.extra))]
#names <- names(zpost.mean)

new_names <- #substr(names, nchar(names)-4, nchar(names))
   paste(substr(names, nchar(names)-4, nchar(names)),
                   rep(c("B", "C"), each= length(names)/2),sep = "_")

names(zpost.mean) <- new_names

# getting y's name
matching <- paste(dat$spec, rep(c("B", "C"), each= 200),sep = "_")

random_eff <- zpost.mean[matching]



# Germany guy way
y_pred <- vector(mode = "character", n)
for(i in 1:n) {

  # get the column for i row's different realization)
  liab_i <- model$Liab[ , i+(0:(K-1))*n] # A-B,... A-C,....
  # adding randomm effect
  liab_i <- liab_i + random_eff[i+(0:(K-1))*n]

  # get prob for B and C (probably this is incorrect)?
  prob_i <- exp(liab_i)/ (1 + rowSums(exp(liab_i)))
  # get prob for A
  prob_all_i <- cbind(1 - rowSums(prob_i), prob_i)
  # average
  prob_mean_i <- colMeans(prob_all_i)
  # choose the biggest prob
  y_pred[i] <- cat_list[which(max(prob_mean_i) == prob_mean_i)]

}

# Jarrod's way
IJ <- (1/3) * (diag(2) + matrix(1, 2, 2))
Delta <- cbind(c(-1, 1, 0), c(-1, 0, 1))
c2 <- (16 * sqrt(3)/(15 * pi))^2
D <- ginv(Delta %*% t(Delta)) %*% Delta
Int <- t(apply(liab_i , 1, function(x) {
  D %*% (x/sqrt(1 + c2 * diag(IJ)))
}))
y_pred2 <- vector(mode = "character", n)
for(i in 1:n) {

  # get the column for i row's different realization)
  liab_i <- model$Liab[ , i+(0:(K-1))*n] # A-B,... A-C,....
  liab_i2 <- liab_i + random_eff[i+(0:(K-1))*n]
  # get prob for B and C (probably this is incorrect)?
  Int <- t(apply(liab_i , 1, function(x) {
    D %*% (x/sqrt(1 + c2 * diag(IJ)))
  }))

  Int2 <- t(apply(liab_i2 , 1, function(x) {
    D %*% (x/sqrt(1 + c2 * diag(IJ)))
  }))
  liab_test <- exp(Int+Int2)/rowSums(exp(Int+Int2))
  #prob_mean_i <- colMeans(liab_test)
  # average
  #prob_mean_i <- colMeans(prob_all_i)
  # choose the biggest prob
  y_pred2[i] <- cat_list[which(max(colMeans(liab_test)) == colMeans(liab_test))]

}

# they should match?
y_pred
y_pred2
dat$x4


