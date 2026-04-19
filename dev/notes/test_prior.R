library(MCMCglmm)
library(ape)

set.seed(123)
phylo <- rtree(10)
phylo <- compute.brlen(phylo, method = "Grafen")

n_total <- 30
data <- data.frame(
  y = rnorm(n_total),
  x = rnorm(n_total),
  Species = rep(phylo$tip.label, each = 3),
  Species2 = rep(phylo$tip.label, each = 3)
)

A <- inverseA(phylo, nodes = "ALL")$Ainv
I_species <- diag(10)
rownames(I_species) <- colnames(I_species) <- phylo$tip.label
I_species_sparse <- as(I_species, "dgCMatrix")

prior1 <- list(
  R = list(V = 1, nu = 0.002),
  G = list(
    list(V = 1, nu = 0.002),
    list(V = 1, nu = 0.002)
  )
)

cat("Testing prior with V=1 for both G elements...\n")
tryCatch({
  m <- MCMCglmm(
    y ~ x,
    random = ~ Species + Species2,
    data = data,
    ginverse = list(Species = A, Species2 = I_species_sparse),
    prior = prior1,
    nitt = 100, thin = 1, burnin = 0,
    verbose = FALSE
  )
  cat("SUCCESS!\n")
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})
