data <- sim_bace(
    n_cases = 500, n_species = 200,
    response_type = "gaussian",
    predictor_types = c("gaussian", "gaussian", "gaussian"),
    missingness = c(0.5, 0.2, 0, 0), phylo_signal = c(0.6, 0.3, 0.2, 0.1)
)

mod1 <- bace_imp(
    data = data$data,
    phylo = data$tree,
    fixformula = "y ~ x1 + x2 + x3",
    ran_phylo_form = "~ 1 |species",
    runs = 10, nitt = 25000, burnin = 7000, thin = 18
)


source("./dev/convergence/convergence_diagnostics.R")
source("./dev/convergence/convergence_plots.R")
conv1 <- assess_convergence(mod1, method = "summary")
plot.bace_convergence(conv1, type = "all")

conv1 <- assess_convergence(mod1, method = "summary.percentage")
plot.bace_convergence(conv1)

conv1 <- assess_convergence(mod1, method = "wasserstein")
conv1
plot_wasserstein_convergence(conv1)

conv1 <- assess_convergence(mod1, method = "energy")
conv1
plot_energy_convergence(conv1)

plot_bace_imputation_comparison(mod1, variable = "y")
