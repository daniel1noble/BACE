# test for caterogical data

library(here)
library(tidyverse)
library(ape)
library(MCMCglmm)


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

model <- MCMCglmm::MCMCglmm(x4 ~ -1 + trait + resp + x1 + x2 + x3, random = ~us(-1+trait + resp):spec, family = "categorical",
                          nitt = 110000, pl = TRUE, pr = TRUE, thin = 10, prior = prior, burnin = 10000,
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

y_pred <- vector(mode = "character", n)
for(i in 1:n) {

  # get the column for i row's different realization)
  liab_i <- model$Liab[ , i+(0:(K-1))*n] # A-B,... A-C,....
  # get prob for B and C (probably this is incorrect)
  prob_i <- exp(liab_i)/ (1 + rowSums(exp(liab_i)))
  # get one for A
  prob_all_i <- cbind(1 - rowSums(prob_i), prob_i)
  # average
  prob_mean_i <- colMeans(prob_all_i)
  # choose the biggest prob
  y_pred[i] <- cat_list[which(max(colMeans(prob_all_i)) == colMeans(prob_all_i))]

}

# they should match?
y_pred
dat$x4


