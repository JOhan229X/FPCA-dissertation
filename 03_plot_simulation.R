###############################################################################
## 03_plot_simulation.R
## ---------------------------------------------------------------------------
## Reconstruction MISE figure for Chapter 3, with the no-missingness baseline
## drawn as a reference line. Reads mc_raw_replications.rds (from
## 01_simulation_main.R). Edit BASE if your baseline values differ.
## Output: fig_mise_recon_final.png
###############################################################################

raw <- readRDS("mc_raw_replications.rds")

## baseline (no-missingness) reconstruction MISE from fpca_baseline_serial.R
BASE <- c(DS1 = 30.0, DS2 = 1.11)

mechs <- c("mcar", "thr_mar", "inc_mar", "thr_mnar", "inc_mnar")
labs  <- c("MCAR", "Threshold\nMAR", "Increasing\nMAR",
           "Threshold\nMNAR", "Increasing\nMNAR")
class <- c("MCAR", "MAR", "MAR", "MNAR", "MNAR")
cols  <- c(MCAR = "#4a8f6d", MAR = "#2A6F97", MNAR = "#BC4749")
pcol  <- cols[class]

val <- function(ds, scen) {
  x <- raw$ise_recon[raw$dataset == ds & raw$scenario == scen]
  x <- x[is.finite(x)]; if (length(x) == 0) NA else mean(x)
}

png("fig_mise_recon_final.png", width = 12, height = 5.4, units = "in", res = 300)
op <- par(mfrow = c(1, 2), mar = c(5, 5, 3, 1.5), oma = c(0, 0, 3.2, 0))
ds_ids   <- c("DS1", "DS2")
ds_title <- c("Dataset 1  (truth = mixed model)", "Dataset 2  (truth = FPCA)")

for (j in 1:2) {
  ds  <- ds_ids[j]
  v30 <- sapply(mechs, function(m) val(ds, sprintf("%s_30", m)))
  v60 <- sapply(mechs, function(m) val(ds, sprintf("%s_60", m)))
  fm  <- val(ds, "fixed_mar_det")
  base <- BASE[[ds]]

  ylim <- range(c(v30, v60, fm, base), na.rm = TRUE)
  ylim[1] <- min(ylim[1], base) - 0.05 * diff(ylim)
  ylim[2] <- ylim[2] + 0.15 * diff(ylim)

  plot(NULL, xlim = c(0.7, 5.3), ylim = ylim, xaxt = "n",
       xlab = "", ylab = "MISE (reconstruction)", main = ds_title[j])
  axis(1, at = 1:5, labels = labs, cex.axis = 0.8, padj = 0.5)
  grid(nx = NA, ny = NULL, col = "grey90")
  rect(3.5, ylim[1], 5.4, ylim[2], col = adjustcolor("#BC4749", 0.05), border = NA)
  text(4.5, ylim[2] - 0.08 * diff(ylim), "MNAR",
       col = adjustcolor("#BC4749", 0.7), cex = 0.85)

  ## baseline (no missing) reference line, bottom
  abline(h = base, lty = 4, lwd = 1.3, col = "#4a8f6d")
  text(0.75, base, sprintf("no missing = %.1f", base),
       adj = c(0, -0.3), cex = 0.6, col = "#4a8f6d")
  ## fixed MAR reference line, top
  abline(h = fm, lty = 3, col = "grey55")
  text(5.3, fm, sprintf("fixed MAR = %.1f", fm),
       adj = c(1, -0.3), cex = 0.6, col = "grey40")

  lines(1:5, v30, lty = 1, lwd = 1.5, col = "grey45")
  lines(1:5, v60, lty = 2, lwd = 1.5, col = "grey45")
  points(1:5, v30, pch = 16, col = pcol, cex = 1.5)
  points(1:5, v60, pch = 15, col = pcol, cex = 1.5)

  legend("topleft", bty = "n", cex = 0.72,
         legend = c("dropout 30%", "dropout 60%", "no missing", "fixed MAR"),
         pch = c(16, 15, NA, NA), lty = c(1, 2, 4, 3),
         col = c("grey45", "grey45", "#4a8f6d", "grey55"))
}
mtext("Reconstruction MISE across mechanisms (100 replications, FVE 99%)",
      side = 3, line = 0.6, outer = TRUE, font = 2, cex = 1.15)
par(op)
dev.off()
cat("saved fig_mise_recon_final.png (with baseline line)\n")
###############################################################################

###############################################################################
## Second figure: MISE of the estimated mean and first eigenfunction (DS2).
## Reads mc_summary.rds (from 01_simulation_main.R). Output: fig_mise_components.png
###############################################################################
summ <- readRDS("mc_summary.rds")
d2 <- summ[summ$dataset == "DS2" & summ$scenario != "fixed_mar_det", ]

ord <- c("mcar_30","mcar_60","thr_mar_30","thr_mar_60","inc_mar_30","inc_mar_60",
         "thr_mnar_30","thr_mnar_60","inc_mnar_30","inc_mnar_60")
d2 <- d2[match(ord, d2$scenario), ]
labs2 <- c("MCAR 30","MCAR 60","thrMAR 30","thrMAR 60","incMAR 30","incMAR 60",
           "thrMNAR 30","thrMNAR 60","incMNAR 30","incMNAR 60")

png("fig_mise_components.png", width = 12, height = 5, units = "in", res = 300)
op <- par(mfrow = c(1, 2), mar = c(7, 4.5, 3, 1))
barplot(d2$MISE_mu, names.arg = labs2, las = 2, col = "#4a8f6d",
        ylab = expression("MISE of " * hat(mu)), main = "Mean function",
        cex.names = 0.7)
barplot(d2$MISE_phi, names.arg = labs2, las = 2, col = "#2A6F97",
        ylab = expression("MISE of " * hat(phi)[1]), main = "First eigenfunction",
        cex.names = 0.7)
par(op); dev.off()
cat("saved fig_mise_components.png\n")
###############################################################################
