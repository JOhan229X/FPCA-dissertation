###############################################################################
## 06_adni_fpca_raw.R
## Same analysis as adni_log_figures.R but on the RAW ADAS11 score (no log).
## Produces:
##   adni_raw_mu_phi.png   -- mean + eigenfunctions
##   adni_raw_mode.png     -- mode of variation for FPC1
##   adni_raw_scores.png   -- distribution of FPC1 scores
###############################################################################

library(fdapace)

## ---- load cleaned data; rebuild from MSc_Data.csv if needed ----------------
if (file.exists("adni_clean.rds")) {
  dat <- readRDS("adni_clean.rds")
} else {
  raw <- read.csv("MSc_Data.csv", stringsAsFactors = FALSE)
  raw$date   <- as.Date(raw$EXAMDATE, format = "%d/%m/%Y")
  raw$ADAS11 <- suppressWarnings(as.numeric(raw$ADAS11))
  dat <- raw[!is.na(raw$ADAS11), c("RID", "date", "ADAS11")]
  dat <- dat[dat$RID %in% names(which(table(dat$RID) >= 2)), ]
  dat <- dat[order(dat$RID, dat$date), ]
  dat$years <- ave(as.numeric(dat$date), dat$RID,
                   FUN = function(x) (x - min(x)) / 365.25)
}
dat <- dat[order(dat$RID, dat$years), ]

## ---- RAW scale: use ADAS11 directly (no transform) -------------------------
ids <- unique(dat$RID)
Ly  <- lapply(ids, function(i) dat$ADAS11[dat$RID == i])   # raw, not log
Lt  <- lapply(ids, function(i) dat$years [dat$RID == i])
keep <- sapply(Ly, length) >= 2
Ly <- Ly[keep]; Lt <- Lt[keep]

fit <- FPCA(Ly, Lt, optns = list(dataType = "Sparse",
                                 FVEthreshold = 0.99, verbose = FALSE))
K <- fit$selectK; lam <- fit$lambda; wg <- fit$workGrid
cat(sprintf("raw-scale FPCA: K(99%%)=%d ; FPC1 = %.1f%% of variance\n",
            K, 100 * lam[1] / sum(lam)))

## FIG 1: mean + eigenfunctions ---------------------------------------------
png("adni_raw_mu_phi.png", width = 11, height = 4.6, units = "in", res = 300)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1.2))
plot(wg, fit$mu, type = "l", lwd = 3, col = "#BC4749",
     xlab = "Years since first visit", ylab = "ADAS11",
     main = "Mean decline (raw scale)"); grid(col = "grey90")
nshow <- min(K, 3)
matplot(wg, fit$phi[, seq_len(nshow), drop = FALSE], type = "l",
        lwd = 2, lty = 1, col = c("#2A6F97", "#4a8f6d", "#d9a441"),
        xlab = "Years since first visit", ylab = expression(phi(t)),
        main = "Principal modes (raw scale)")
abline(h = 0, col = "grey60", lty = 3); grid(col = "grey90")
legend("topright", bty = "n", cex = 0.8,
       legend = paste0("FPC", seq_len(nshow),
                       sprintf(" (%.0f%%)", 100 * lam[seq_len(nshow)] / sum(lam))),
       col = c("#2A6F97", "#4a8f6d", "#d9a441"), lwd = 2)
par(op); dev.off(); cat("saved adni_raw_mu_phi.png\n")

## FIG 2: mode of variation --------------------------------------------------
sd1 <- sqrt(lam[1])
hi <- fit$mu + 2 * sd1 * fit$phi[, 1]; lo <- fit$mu - 2 * sd1 * fit$phi[, 1]
png("adni_raw_mode.png", width = 7.5, height = 5, units = "in", res = 300)
op <- par(mar = c(4.5, 4.5, 3, 1))
plot(wg, fit$mu, type = "l", lwd = 3, col = "black", ylim = range(c(hi, lo)),
     xlab = "Years since first visit", ylab = "ADAS11",
     main = "Mode of variation (raw scale)")
lines(wg, hi, col = "#BC4749", lwd = 2, lty = 2)
lines(wg, lo, col = "#2A6F97", lwd = 2, lty = 2); grid(col = "grey90")
legend("topleft", bty = "n", cex = 0.85,
       legend = c("mean", "mean + 2SD x FPC1", "mean - 2SD x FPC1"),
       col = c("black", "#BC4749", "#2A6F97"), lwd = 2, lty = c(1, 2, 2))
par(op); dev.off(); cat("saved adni_raw_mode.png\n")

## FIG 3: score distribution -------------------------------------------------
xi1 <- fit$xiEst[, 1]
png("adni_raw_scores.png", width = 7, height = 4.6, units = "in", res = 300)
op <- par(mar = c(4.5, 4.5, 3, 1))
hist(xi1, breaks = 30, col = "#9dc3d4", border = "white",
     xlab = "FPC1 score", main = "FPC1 scores (raw scale)")
abline(v = 0, col = "#BC4749", lwd = 2)
par(op); dev.off(); cat("saved adni_raw_scores.png\n")

cat(sprintf("\nraw-scale FPC1 = %.1f%% (compare log-scale 96.0%%)\n",
            100 * lam[1] / sum(lam)))
###############################################################################
