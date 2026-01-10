# test 1
## Packages used internally
library(MCMCglmm)
library(ape)
library(dplyr)
library(magrittr)
library(BACE)

set.seed(1)

## --- 1) Simulate data with missingness ---
sim <- sim_bace(
  response_type   = "gaussian",
  predictor_types = c("gaussian", "binary",  "threshold3"),
  n_cases         = 400,
  n_species       = 400,
  phylo_signal    = c(0.4, 0.2, 0.1, 0.1),   
  missingness     = c(0.25, 0.15, 0.15, 0.10) 
)

dat  <- sim$data
tree <- sim$tree

## --- 2) IMPORTANT: convert integer-coded categorical vars to factors ---
## simBACE returns:
## - multinomialK as character categories (OK; but I'd still factor() it)
## - binary / thresholdK as integers -> MUST convert or bace_imp will treat as gaussian/poisson

# binary predictor here is x2 (because predictor_types[2] == "binary")
dat$x2 <- factor(dat$x2, levels = sort(unique(dat$x2[!is.na(dat$x2)])))

# multinomial4 predictor here is x3 (characters like A/B/C/D)
#dat$x3 <- factor(dat$x3)

# threshold3 predictor here is x4 (integers 1/2/3) -> ordered factor
dat$x3 <- ordered(dat$x3, levels = sort(unique(dat$x3[!is.na(dat$x3)])))

## If you ever simulate binary/threshold RESPONSE, do the same to dat$y.

## --- 3) Run BACE imputation ---
## NOTE: simBACE uses column name "species" (lowercase),
## so use "~ 1 | species" (not "~ 1 | Species").
fit <- bace_imp(
  fixformula      = "y ~ x1 + x2 + x3",
  ran_phylo_form  = "~ 1 | species",
  phylo           = tree,
  data            = dat,
  runs            = 3,
  nitt            = 30000,
  burnin          = 10000,
  thin            = 100,
  verbose         = TRUE
)

## --- 4) Inspect outputs ---
names(fit)
names(fit$data)           # "Initial_Data", "Iter_1", "Iter_2", "Iter_3", ...

imp3 <- fit$data$Iter_3
colSums(is.na(imp3))      # should be ~0 for variables included in fixformula
str(imp3)