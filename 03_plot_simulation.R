###############################################################################
## 03_plot_simulation.R
## ---------------------------------------------------------------------------
## Reconstruction MISE figure for Chapter 3, with the no-missingness baseline
## drawn as a reference line. Reads mc_raw_replications.rds (from
## 01_simulation_main.R). Edit BASE if your baseline values differ.
## Output: fig_mise_recon_final.png
###############################################################################

raw  <- readRDS("mc_raw_replications.rds")
BASE <- c(DS1 = 30.1, DS2 = 1.11)     # no-missingness baseline (02_baseline.R)
 
## x-axis order: baseline first, calibrated mechanisms in the middle, fixed last
xorder <- c("mcar", "thr_mar", "inc_mar", "thr_mnar", "inc_mnar")
xlabs  <- c("no\nmissing", "MCAR", "Threshold\nMAR", "Increasing\nMAR",
            "Threshold\nMNAR", "Increasing\nMNAR", "Fixed\nMAR")
 
val <- function(ds, scen) {
  x <- raw$ise_recon[raw$dataset == ds & raw$scenario == scen]
  x <- x[is.finite(x)]; if (length(x) == 0) NA else mean(x)
}
 
png("fig_mise_recon_final.png", width = 12, height = 5.4, units = "in", res = 300)
op <- par(mfrow = c(1, 2), mar = c(5, 5, 3, 1.5), oma = c(0, 0, 3.2, 0))
ds_ids   <- c("DS1", "DS2")
ds_title <- c("Dataset 1  (truth = mixed model)", "Dataset 2  (truth = FPCA)")
 
for (j in 1:2) {
  ds <- ds_ids[j]
  ## build the sequence of 7 points in x-order
  mid30 <- sapply(xorder, function(m) val(ds, sprintf("%s_30", m)))
  mid60 <- sapply(xorder, function(m) val(ds, sprintf("%s_60", m)))
  fm    <- val(ds, "fixed_mar_det")
  base  <- BASE[[ds]]
  y30 <- c(base, mid30, fm)     # baseline, 5 mechanisms at 30%, fixed
  y60 <- c(base, mid60, fm)     # baseline, 5 mechanisms at 60%, fixed
  x   <- seq_along(xlabs)       # 1..7
 
  ## y-limits: include everything, but check compression. If fixed/baseline are
  ## extreme, the middle still shows because the line + points span the range.
  ylim <- range(c(y30, y60), na.rm = TRUE)
  ylim[2] <- ylim[2] + 0.12 * diff(ylim)
 
  plot(NULL, xlim = c(0.7, 7.3), ylim = ylim, xaxt = "n",
       xlab = "", ylab = "MISE (reconstruction)", main = ds_title[j])
  axis(1, at = x, labels = xlabs, cex.axis = 0.72, padj = 0.5)
  grid(nx = NA, ny = NULL, col = "grey90")
 
  ## shade the MNAR region (positions 5-6)
  rect(4.5, ylim[1], 6.5, ylim[2], col = adjustcolor("#BC4749", 0.05), border = NA)
 
  ## connect all seven points with a line so the trend is visible
  lines(x, y30, lty = 1, lwd = 1.6, col = "grey45")
  lines(x, y60, lty = 2, lwd = 1.6, col = "grey45")
 
  ## colour points by class: baseline (green), MCAR (teal), MAR (blue),
  ## MNAR (red), fixed (grey)
  pcol <- c("#4a8f6d", "#4a8f6d", "#2A6F97", "#2A6F97",
            "#BC4749", "#BC4749", "#777777")
  points(x, y30, pch = 16, col = pcol, cex = 1.5)
  points(x, y60, pch = 15, col = pcol, cex = 1.5)
 
  legend("topleft", bty = "n", cex = 0.72,
         legend = c("dropout 30%", "dropout 60%"),
         pch = c(16, 15), lty = c(1, 2), col = "grey45")
}
mtext("Reconstruction MISE across mechanisms (100 replications, FVE 99%)",
      side = 3, line = 0.6, outer = TRUE, font = 2, cex = 1.15)
par(op); dev.off()
cat("saved fig_mise_recon_final.png (connected trend, baseline -> fixed)\n")
 
###############################################################################
## Second figure: MISE of the mean and first eigenfunction (Dataset 2),
## drawn as connected trends to match the reconstruction figure.
###############################################################################
summ <- readRDS("mc_summary.rds")
 
## baseline mu/phi from 02_baseline.R (no-missingness); edit if yours differ
BASE_MU  <- 0.578
BASE_PHI <- 0.067
 
xorder <- c("mcar", "thr_mar", "inc_mar", "thr_mnar", "inc_mnar")
xlabs  <- c("no\nmissing", "MCAR", "Threshold\nMAR", "Increasing\nMAR",
            "Threshold\nMNAR", "Increasing\nMNAR")
pcol   <- c("#4a8f6d", "#4a8f6d", "#2A6F97", "#2A6F97", "#BC4749", "#BC4749")
 
getv <- function(scen, col) {
  v <- summ[summ$dataset == "DS2" & summ$scenario == scen, col]
  if (length(v) == 0) NA else v[1]
}
 
png("fig_mise_components.png", width = 12, height = 5.2, units = "in", res = 300)
op <- par(mfrow = c(1, 2), mar = c(5, 5, 3, 1.2))
 
for (comp in c("MISE_mu", "MISE_phi")) {
  base <- if (comp == "MISE_mu") BASE_MU else BASE_PHI
  y30 <- c(base, sapply(xorder, function(m) getv(sprintf("%s_30", m), comp)))
  y60 <- c(base, sapply(xorder, function(m) getv(sprintf("%s_60", m), comp)))
  x   <- seq_along(xlabs)
  ylim <- range(c(y30, y60), na.rm = TRUE); ylim[2] <- ylim[2] + 0.12*diff(ylim)
  ylab <- if (comp == "MISE_mu") expression("MISE of " * hat(mu))
          else                    expression("MISE of " * hat(phi)[1])
  ttl  <- if (comp == "MISE_mu") "Mean function" else "First eigenfunction"
 
  plot(NULL, xlim = c(0.7, 6.3), ylim = ylim, xaxt = "n",
       xlab = "", ylab = ylab, main = ttl)
  axis(1, at = x, labels = xlabs, cex.axis = 0.72, padj = 0.5)
  grid(nx = NA, ny = NULL, col = "grey90")
  rect(4.5, ylim[1], 6.5, ylim[2], col = adjustcolor("#BC4749", 0.05), border = NA)
  pcol7 <- c("#4a8f6d", pcol)
  lines(x, y30, lty = 1, lwd = 1.6, col = "grey45")
  lines(x, y60, lty = 2, lwd = 1.6, col = "grey45")
  points(x, y30, pch = 16, col = pcol7, cex = 1.5)
  points(x, y60, pch = 15, col = pcol7, cex = 1.5)
  legend("topleft", bty = "n", cex = 0.72,
         legend = c("dropout 30%", "dropout 60%"),
         pch = c(16, 15), lty = c(1, 2), col = "grey45")
}
par(op); dev.off()
cat("saved fig_mise_components.png (connected trend)\n")
###############################################################################
