# Campaign figures in the style of ms Figures 3 and 4, for the 2026-11 design
# (n = 100/500; signal 0.3/0.7/0.9; regime floor/typical/good_covariate; mechanism mcar/mar/clade).
# Points = means over replicates, bars = +/- 1.96 Monte Carlo SE. Failed fits excluded (as in the ms).
# Usage: Rscript figures.R <results dir with per_dataset_*.csv and the .rds files> <out dir>
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
a <- commandArgs(TRUE); dir <- a[1]; out <- a[2]; dir.create(out, showWarnings = FALSE, recursive = TRUE)

rec <- read.csv(file.path(dir, "per_dataset_recovery.csv"))
slo <- read.csv(file.path(dir, "per_dataset_slopes.csv"))

## Per-value 95% interval coverage (numeric traits) is in the .rds files only, not in the CSVs.
fs <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
cov <- do.call(rbind, lapply(fs, function(f) {
  x <- readRDS(f); cc <- x$cell_coverage
  if (is.null(cc) || !nrow(cc)) return(NULL)
  data.frame(n_species = x$scenario$n_species, phylo_signal = x$scenario$phylo_signal,
             regime = x$scenario$regime, mechanism = x$scenario$mechanism,
             trait = cc$variable, value = cc$coverage)
}))

lab_trait <- c(x1 = "continuous 1", x2 = "continuous 2", cnt = "count",
               bin = "binary", ord = "ordinal (4)", cat3 = "categorical (3)")
lab_reg <- c(floor = "floor", typical = "typical", good_covariate = "good covariate")
summ <- function(d, by) {
  d <- d[!is.na(d$value), c(by, "value")]
  f <- function(fun) stats::aggregate(d["value"], d[by], fun)
  r <- f(mean); names(r)[names(r) == "value"] <- "mean"
  r$mcse <- f(function(v) sd(v) / sqrt(length(v)))$value
  r$n <- f(length)$value
  for (b in by) r[[b]] <- as.character(r[[b]])
  r
}
by <- c("n_species", "phylo_signal", "regime", "mechanism", "trait")
style <- function(p) p + theme_bw(base_size = 9) +
  scale_colour_viridis_d(end = 0.8, name = "phylogenetic signal") +
  scale_shape_manual(values = c(mcar = 16, mar = 17, clade = 15), name = "missingness") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
prep <- function(s) {
  s$n_species <- factor(s$n_species, c("100", "500"))
  s$trait <- factor(lab_trait[s$trait], lab_trait)
  s$regime <- factor(lab_reg[s$regime], lab_reg)
  s$mechanism <- factor(s$mechanism, c("mcar", "mar", "clade"))
  s
}
pd <- position_dodge(width = 0.75)
pt <- function(s, y) ggplot(s, aes(n_species, mean, colour = phylo_signal, shape = mechanism,
                                    group = interaction(phylo_signal, mechanism))) +
  geom_errorbar(aes(ymin = mean - 1.96 * mcse, ymax = mean + 1.96 * mcse), width = 0, position = pd, linewidth = 0.3) +
  geom_point(position = pd, size = 1.4) + facet_grid(regime ~ trait) + labs(x = "Number of species", y = y)

## Figure A (cf. ms Fig. 3): coverage, zRMSE, accuracy.
cA <- prep(summ(cov, by))
cB <- prep(summ(transform(rec[rec$kind == "numeric", ], value = zrmse), by))
cC <- prep(summ(transform(rec[rec$kind == "category", ], value = accuracy), by))
chance <- data.frame(trait = factor(lab_trait[c("bin", "ord", "cat3")], lab_trait), y = c(1/2, 1/4, 1/3))
pA <- style(pt(cA, "95% interval coverage") + geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40"))
pB <- style(pt(cB, "Standardised RMSE") + geom_hline(yintercept = 1, linetype = 3, colour = "grey40"))
pC <- style(pt(cC, "Accuracy") + geom_hline(data = chance, aes(yintercept = y), linetype = 3, colour = "grey40"))
figA <- (pA / (pB | pC)) + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(out, "fig_imputation.pdf"), figA, width = 11, height = 10)
ggsave(file.path(out, "fig_imputation.png"), figA, width = 11, height = 10, dpi = 200)

## Figure B (cf. ms Fig. 4): downstream slope, BACE vs complete cases vs oracle.
## The oracle is counted once per distinct dataset (clade/mar/mcar share the complete data).
s2 <- slo
s2$key <- paste(s2$arm, s2$n_species, s2$phylo_signal, s2$regime, s2$seed)
s2 <- s2[!(s2$arm == "oracle" & duplicated(s2$key)), ]
s2$mechanism[s2$arm == "oracle"] <- "mcar"            # oracle has no mechanism; plotted once
s2$trait <- s2$arm
lab_arm <- c(bace = "BACE", complete_case = "complete cases", oracle = "oracle (no missing data)")
mk <- function(v) { s <- summ(transform(s2, value = v), by); s$n_species <- factor(s$n_species, c("100", "500"))
  s$arm <- factor(lab_arm[s$trait], lab_arm); s$regime <- factor(lab_reg[s$regime], lab_reg)
  s$mechanism <- factor(s$mechanism, c("mcar", "mar", "clade")); s }
sb <- mk(s2$bias); sc <- mk(as.numeric(s2$covered))
pt2 <- function(s, y) ggplot(s, aes(n_species, mean, colour = phylo_signal, shape = mechanism,
                                     group = interaction(phylo_signal, mechanism))) +
  geom_errorbar(aes(ymin = mean - 1.96 * mcse, ymax = mean + 1.96 * mcse), width = 0, position = pd, linewidth = 0.3) +
  geom_point(position = pd, size = 1.4) + facet_grid(regime ~ arm) + labs(x = "Number of species", y = y)
pD <- style(pt2(sb, "Bias of slope") + geom_hline(yintercept = 0, linetype = 2, colour = "grey40"))
pE <- style(pt2(sc, "95% CI coverage") + geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40"))
figB <- (pD | pE) + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(out, "fig_downstream.pdf"), figB, width = 11, height = 6)
ggsave(file.path(out, "fig_downstream.png"), figB, width = 11, height = 6, dpi = 200)

write.csv(rbind(cbind(measure = "coverage", cA), cbind(measure = "zrmse", cB), cbind(measure = "accuracy", cC)),
          file.path(out, "fig_imputation_data.csv"), row.names = FALSE)
write.csv(rbind(cbind(measure = "bias", sb), cbind(measure = "coverage", sc)),
          file.path(out, "fig_downstream_data.csv"), row.names = FALSE)
cat("wrote figures to", out, "\n")
