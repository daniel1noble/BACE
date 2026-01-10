simdata <- sim_bace(
    response_type = "gaussian",
    predictor_types = c("gaussian", "binary", "multinomial4", "poisson", "gaussian", "gaussian"),
    beta_resp = c(0.6, 0.4, -0.3, 0.3, 0.1, 0.8),
    intercepts = c(0, 0.2, 0.2, 0, -0.5, -0.3, 0),
    phylo_signal = c(0.8, 0.5, 0.3, 0.6, 0.9, 0.4, 0.2),
    n_cases = 750,
    n_species = 200,
    missingness = c(0.5, 0.3, 0, 0.1, 0.2, 0.2, 0)
)
results <- bace_imp(
    data = simdata$data,
    ran_phylo_form = "~ 1 | species",
    phylo = simdata$tree,
    runs = 50,
    fixformula = "y ~ x1 + x2 + x3 + x4 + x5"
)


conv <- assess_convergence(results, method = "summary", use_all_data = TRUE)
print(conv)
plot.bace_convergence(conv, type = "all")
plot_convergence_summary(conv)
plot_trace_convergence(conv)
plot_density_convergence(conv)
plot_bace_imputation_comparison(results, "x5")

conv <- assess_convergence(results, method = "energy")
print(conv)

conv <- assess_convergence(results, method = "wasserstein")
print(conv)
plot_wasserstein_convergence(conv)
