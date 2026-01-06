simdata <- simBACE_gaussian(beta_sparsity = 0.5, n_predictors = 5)
print_simBACE_summary(simdata)


simdata2 <- simBACE(
  response_type = "gaussian",
  predictor_types = c("gaussian", "binary", "multinomial3", "gaussian"),
  n_cases = 150,
  n_species = 50,
  phylo_signal = c(0.8, 0.5, 0.2, 0.2, 0.7),
  beta_sparsity = 0.3,
  missingness = c(0, 0.3, 0.4, 0.5, 0.2)
)
print_simBACE_summary(simdata2)


simdata3 <- simBACE(
  response_type = "ordinal4",  # 4-level ordinal response (1,2,3,4)
  predictor_types = c("gaussian", "binary"),
  n_cases = 200, n_species = 75, beta_sparsity = 0
)
print_simBACE_summary(simdata3)
