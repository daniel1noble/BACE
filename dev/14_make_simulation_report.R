# =============================================================================
# 14_make_simulation_report.R
#   Read the saved simulation RDS outputs and produce (a) human-readable
#   per-replicate CSVs and (b) figures for the report.
#
#   Figures follow simulation-study convention (Morris, White & Crowther 2019,
#   Stat Med 38:2074): POINT ESTIMATES with Monte Carlo standard-error whiskers,
#   computed from the per-replicate values. Whiskers = +/- 1.96 * MCSE
#   (a 95% Monte Carlo interval for the estimate). Replicate count is stated
#   on every panel. Colours: Okabe-Ito colourblind-safe palette.
# =============================================================================

OUT <- file.path("dev", "simulation_results")
FIG <- file.path(OUT, "figures")
if (!dir.exists(FIG)) dir.create(FIG, recursive = TRUE)

A <- readRDS(file.path(OUT, "recovery_simulation.rds"))       # Study A (dev/12)
B <- readRDS(file.path(OUT, "response_type_simulation.rds"))  # Study B (dev/13)

# per-replicate CSVs (raw, inspectable)
utils::write.csv(A$per_rep, file.path(OUT, "recovery_per_rep.csv"), row.names = FALSE)
utils::write.csv(A$summary, file.path(OUT, "recovery_summary.csv"), row.names = FALSE)
utils::write.csv(B$per_rep, file.path(OUT, "response_type_per_rep.csv"), row.names = FALSE)
utils::write.csv(B$summary, file.path(OUT, "response_type_summary.csv"), row.names = FALSE)

A_REPS <- A$cfg$reps      # replicates per scenario, Study A
B_REPS <- B$cfg$reps      # replicates per scenario, Study B

# ---- palette + MCSE + dot-whisker helpers -----------------------------------
OI <- c(blue = "#0072B2", orange = "#E69F00", green = "#009E73",
        vermillion = "#D55E00", purple = "#CC79A7", sky = "#56B4E9",
        grey = "#8A8A8A", black = "#000000")
INK <- "#222222"; MUTE <- "#666666"; GRID <- "#E9E9E9"
MECH_COL <- c(MCAR = unname(OI["sky"]), MAR = unname(OI["vermillion"]))
Z <- 1.96

mcse_mean <- function(v) { v <- v[is.finite(v)]; stats::sd(v) / sqrt(length(v)) }
mcse_prop <- function(v) { v <- v[is.finite(v)]; p <- mean(v); sqrt(p * (1 - p) / length(v)) }

# est + mcse for a metric column over a subset of a per-rep frame.
# robust = TRUE uses the median and a MAD-based MCSE (SE of the median ~
# 1.2533 * MAD / sqrt(n)) — needed for poisson bias, which has a heavy
# extreme-count tail that corrupts the mean.
est_mcse <- function(df, col, prop = FALSE, robust = FALSE) {
  v <- df[[col]]; v <- v[is.finite(v)]
  if (robust) return(list(est = stats::median(v),
                          mcse = 1.2533 * stats::mad(v) / sqrt(length(v))))
  list(est = mean(v), mcse = if (prop) mcse_prop(v) else mcse_mean(v))
}

png_open <- function(f, w = 1500, h = 900) png(file.path(FIG, f), width = w,
  height = h, res = 190)

# one point-range at vertical position y
add_pr <- function(y, est, mcse, col) {
  segments(est - Z * mcse, y, est + Z * mcse, y, col = col, lwd = 2.6, lend = 1)
  points(est, y, pch = 21, bg = col, col = "white", cex = 1.6, lwd = 1.3)
}
vgrid <- function(at) abline(v = at, col = GRID, lwd = 1)

# =============================================================================
# Study A — parameter recovery (gaussian, known slope b1 = 0.8)
# =============================================================================
Ap <- A$per_rep
Ap$mean_bias[Ap$method == "oracle"] <- 0                    # oracle mean bias = 0 by def
meths <- c("bace", "complete_case", "oracle")               # bottom -> top
mlab  <- c(bace = "BACE", complete_case = "complete-case", oracle = "oracle")
ypos  <- stats::setNames(seq_along(meths), meths)
dodge <- c(MCAR = 0.16, MAR = -0.16)

panelA <- function(col, prop, xlab, main, ref, xlim, refj = FALSE) {
  plot(NULL, xlim = xlim, ylim = c(0.5, length(meths) + 0.5), yaxt = "n",
       xlab = xlab, ylab = "", main = main, col.lab = INK, col.main = INK,
       font.main = 1, col.axis = MUTE, cex.axis = 0.9)
  vgrid(pretty(xlim)); abline(v = ref, lty = 2, col = INK, lwd = 1.4)
  for (m in meths) for (mech in c("MCAR", "MAR")) {
    sub <- Ap[Ap$method == m & Ap$mechanism == mech, ]
    em <- est_mcse(sub, col, prop = prop)
    add_pr(ypos[m] + dodge[mech], em$est, em$mcse, MECH_COL[mech])
  }
  axis(2, at = ypos, labels = mlab[names(ypos)], las = 1, col.axis = INK,
       tick = FALSE)
}

png_open("figA_parameter_recovery.png", w = 1650, h = 820)
par(mfrow = c(1, 2), mar = c(4, 6.6, 2.6, 1.4), mgp = c(2.4, 0.6, 0),
    tcl = -0.3, oma = c(3, 0, 2.4, 0))
panelA("covered", TRUE, "95% CI coverage of true slope",
       "Slope coverage (nominal = 0.95)", ref = 0.95, xlim = c(0.70, 1.03))
panelA("slope_est", FALSE, "estimated slope (true b1 = 0.8)",
       "Slope estimate", ref = 0.8, xlim = c(0.72, 1.0))
mtext(sprintf("Study A: parameter recovery — %d replicates per scenario. Points = mean; whiskers = ±1.96 × Monte Carlo SE.",
              A_REPS), side = 1, outer = TRUE, line = 1.2, cex = 0.75, col = MUTE)
par(fig = c(0, 1, 0, 1), new = TRUE, mar = rep(0, 4)); plot.new()
legend("top", horiz = TRUE, pch = 21, pt.bg = MECH_COL, col = "white",
       pt.cex = 1.5, legend = names(MECH_COL), bty = "n", text.col = INK,
       cex = 0.9)
dev.off()

# Fig A2: marginal-mean bias (the MAR signature) — its own panel for emphasis
png_open("figA2_mean_bias.png", w = 1400, h = 760)
par(mar = c(4.2, 6.4, 3.4, 2), mgp = c(2.5, 0.6, 0), tcl = -0.3)
xlim <- c(-0.32, 0.06)
plot(NULL, xlim = xlim, ylim = c(0.5, length(meths) + 0.5), yaxt = "n",
     xlab = "bias in estimated mean of y", ylab = "",
     main = sprintf("Study A: marginal-mean bias — BACE corrects MAR selection (%d reps)", A_REPS),
     col.lab = INK, col.main = INK, font.main = 1, col.axis = MUTE)
vgrid(pretty(xlim)); abline(v = 0, lty = 2, col = INK, lwd = 1.4)
for (m in meths) for (mech in c("MCAR", "MAR")) {
  sub <- Ap[Ap$method == m & Ap$mechanism == mech, ]
  em <- est_mcse(sub, "mean_bias")
  add_pr(ypos[m] + dodge[mech], em$est, em$mcse, MECH_COL[mech])
}
axis(2, at = ypos, labels = mlab[names(ypos)], las = 1, col.axis = INK, tick = FALSE)
legend("bottomleft", pch = 21, pt.bg = MECH_COL, col = "white", pt.cex = 1.4,
       legend = names(MECH_COL), bty = "n", text.col = INK, cex = 0.82)
mtext("points = mean, whiskers = ±1.96×MCSE", side = 1, line = 2.6,
      cex = 0.72, col = MUTE)
dev.off()

# =============================================================================
# Study B — response-type recovery
# =============================================================================
Bp <- B$per_rep[B$per_rep$status == "ok", ]
ALL5 <- c("categorical", "ordinal", "binary", "poisson", "gaussian")  # bottom->top

# generic horizontal dot-whisker over (type x mechanism) for one metric column
panelB <- function(types, col, xlab, main, xlim, ref = NULL, chance = NULL,
                   robust = FALSE) {
  yp <- stats::setNames(seq_along(types), types)
  plot(NULL, xlim = xlim, ylim = c(0.5, length(types) + 0.5), yaxt = "n",
       xlab = xlab, ylab = "", main = main, col.lab = INK, col.main = INK,
       font.main = 1, col.axis = MUTE, cex.axis = 0.9)
  vgrid(pretty(xlim))
  if (!is.null(ref)) abline(v = ref, lty = 2, col = INK, lwd = 1.4)
  if (!is.null(chance)) for (i in seq_along(types)) {
    ch <- chance[types[i]]
    if (!is.na(ch)) segments(ch, yp[i] - 0.34, ch, yp[i] + 0.34, lty = 3,
                             col = MUTE, lwd = 1.3)
  }
  for (ty in types) for (mech in c("MCAR", "MAR")) {
    sub <- Bp[Bp$type == ty & Bp$mechanism == mech, ]
    em <- est_mcse(sub, col, robust = robust)
    add_pr(yp[ty] + dodge[mech], em$est, em$mcse, MECH_COL[mech])
  }
  axis(2, at = yp, labels = types, las = 1, col.axis = INK, tick = FALSE)
}
topleg <- function() { par(fig = c(0, 1, 0, 1), new = TRUE, mar = rep(0, 4))
  plot.new(); legend("top", horiz = TRUE, pch = 21, pt.bg = MECH_COL,
    col = "white", pt.cex = 1.5, legend = names(MECH_COL), bty = "n",
    text.col = INK, cex = 0.9) }

# Fig B1: recovery (continuous correlation | discrete balanced accuracy)
png_open("figB1_recovery_by_type.png", w = 1700, h = 850)
par(mfrow = c(1, 2), mar = c(4, 6.4, 2.6, 1), mgp = c(2.4, 0.6, 0),
    tcl = -0.3, oma = c(3, 0, 2.4, 0))
panelB(c("poisson", "gaussian"), "recovery",
       "correlation (imputed vs truth)", "Continuous / count: cell recovery",
       xlim = c(0, 1))
panelB(c("categorical", "binary", "ordinal"), "balanced_accuracy",
       "balanced accuracy", "Discrete: class recovery", xlim = c(0.25, 1),
       chance = c(binary = 0.5, ordinal = 1/3, categorical = 1/3))
mtext(sprintf("Study B: recovery by response distribution — %d replicates per scenario. Points = mean; whiskers = ±1.96 × MCSE; dotted line = chance.",
              B_REPS), side = 1, outer = TRUE, line = 1.2, cex = 0.75, col = MUTE)
topleg(); dev.off()

# Fig B2: BIAS and INTERVAL CALIBRATION for ALL response types
png_open("figB2_bias_calibration.png", w = 1750, h = 830)
par(mfrow = c(1, 2), mar = c(4, 6.4, 2.6, 1), mgp = c(2.4, 0.6, 0),
    tcl = -0.3, oma = c(3.4, 0, 2.4, 0))
# bias uses the MEDIAN (robust): poisson has a heavy extreme-count tail that
# makes the mean bias meaningless. x-range from the plotted whisker extremes.
bext <- unlist(lapply(ALL5, function(ty) lapply(c("MCAR", "MAR"), function(m) {
  e <- est_mcse(Bp[Bp$type == ty & Bp$mechanism == m, ], "bias", robust = TRUE)
  c(e$est - Z * e$mcse, e$est + Z * e$mcse) })))
brange <- range(c(0, bext), na.rm = TRUE) + c(-0.02, 0.02)
panelB(ALL5, "bias", "median standardised bias  (categorical: TV)",
       "Bias (nominal = 0)", xlim = brange, ref = 0, robust = TRUE)
panelB(ALL5, "cover95", "95% interval / set coverage",
       "Interval calibration (nominal = 0.95)", xlim = c(0.4, 1.02), ref = 0.95)
mtext(sprintf("Study B: bias & interval calibration — %d reps/scenario. Points = median (bias) / mean (coverage); whiskers = ±1.96 × MCSE. Categorical bias = total-variation distance.",
              B_REPS), side = 1, outer = TRUE, line = 1.1, cex = 0.70, col = MUTE)
topleg(); dev.off()

# Fig B3: timing by response type (averaged over mechanisms; MCSE whiskers)
png_open("figB3_timing.png", w = 1300, h = 720)
par(mar = c(4.2, 6.4, 3, 1.5), mgp = c(2.5, 0.6, 0), tcl = -0.3)
yp <- stats::setNames(seq_along(ALL5), ALL5)
xr <- c(0, max(Bp$runtime_sec, na.rm = TRUE) * 1.15)
plot(NULL, xlim = xr, ylim = c(0.5, length(ALL5) + 0.5), yaxt = "n",
     xlab = "BACE runtime (s per dataset)", ylab = "",
     main = sprintf("Study B: timing by response type (%d reps; mean ± 1.96×MCSE)", B_REPS),
     col.lab = INK, col.main = INK, font.main = 1, col.axis = MUTE)
vgrid(pretty(xr))
for (ty in ALL5) {
  em <- est_mcse(Bp[Bp$type == ty, ], "runtime_sec")
  col <- if (ty == "categorical") unname(OI["vermillion"]) else unname(OI["grey"])
  add_pr(yp[ty], em$est, em$mcse, col)
  text(em$est + Z * em$mcse, yp[ty], sprintf("  %.1fs", em$est), adj = 0,
       cex = 0.75, col = INK, xpd = NA)
}
axis(2, at = yp, labels = ALL5, las = 1, col.axis = INK, tick = FALSE)
dev.off()

cat("Figures written to", FIG, ":\n"); print(list.files(FIG, pattern = "\\.png$"))
cat("CSVs written to", OUT, "\n")
