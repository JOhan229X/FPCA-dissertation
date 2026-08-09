###############################################################################
## 05_adni_fpca_log.R
## Produces the three figures the Analysis chapter references on the LOG scale:
##   adni_fig_mu_phi.png   -- mean decline + eigenfunctions
##   adni_fig_mode.png     -- mode of variation for FPC1
##   adni_fig_scores.png   -- distribution of FPC1 scores
##
## Uses log(ADAS11 + 1). Run after adni_step1_explore.R (needs adni_clean.rds),
## or it will rebuild the cleaned data from MSc_Data.csv if that file is present.
###############################################################################

library(fdapace)

## ---- load cleaned data (RID, date, ADAS11, years); rebuild if needed -------
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

## ---- LOG transform: log(ADAS11 + 1) ----------------------------------------
dat$logadas <- log(dat$ADAS11 + 1)

## ---- reshape to fdapace list format ----------------------------------------
ids <- unique(dat$RID)
Ly  <- lapply(ids, function(i) dat$logadas[dat$RID == i])
Lt  <- lapply(ids, function(i) dat$years  [dat$RID == i])
keep <- sapply(Ly, length) >= 2
Ly <- Ly[keep]; Lt <- Lt[keep]

## ---- fit FPCA on the log scale ---------------------------------------------
fit <- FPCA(Ly, Lt, optns = list(dataType = "Sparse",
                                 FVEthreshold = 0.99, verbose = FALSE))
K   <- fit$selectK
lam <- fit$lambda
wg  <- fit$workGrid
cat(sprintf("log-scale FPCA: K(99%%)=%d ; FPC1 = %.1f%% of variance\n",
            K, 100 * lam[1] / sum(lam)))

## ===========================================================================
## FIG 1: mean decline + eigenfunctions  ->  adni_fig_mu_phi.png
##   wide device + short titles so nothing is cut
## ===========================================================================
png("adni_fig_mu_phi.png", width = 11, height = 4.6, units = "in", res = 300)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1.2))

## (a) mean
plot(wg, fit$mu, type = "l", lwd = 3, col = "#BC4749",
     xlab = "Years since first visit", ylab = "log(ADAS11 + 1)",
     main = "Mean decline")
grid(col = "grey90")

## (b) eigenfunctions (up to 3)
nshow <- min(K, 3)
matplot(wg, fit$phi[, seq_len(nshow), drop = FALSE], type = "l",
        lwd = 2, lty = 1, col = c("#2A6F97", "#4a8f6d", "#d9a441"),
        xlab = "Years since first visit", ylab = expression(phi(t)),
        main = "Principal modes")
abline(h = 0, col = "grey60", lty = 3); grid(col = "grey90")
legend("topright", bty = "n", cex = 0.8,
       legend = paste0("FPC", seq_len(nshow),
                       sprintf(" (%.0f%%)", 100 * lam[seq_len(nshow)] / sum(lam))),
       col = c("#2A6F97", "#4a8f6d", "#d9a441"), lwd = 2)
par(op)
dev.off()
cat("saved adni_fig_mu_phi.png\n")

## ===========================================================================
## FIG 2: mode of variation for FPC1  ->  adni_fig_mode.png
##   mean and mean +/- 2 sqrt(lambda_1) * phi_1
## ===========================================================================
sd1 <- sqrt(lam[1])
hi  <- fit$mu + 2 * sd1 * fit$phi[, 1]
lo  <- fit$mu - 2 * sd1 * fit$phi[, 1]

png("adni_fig_mode.png", width = 7.5, height = 5, units = "in", res = 300)
op <- par(mar = c(4.5, 4.5, 3, 1))
plot(wg, fit$mu, type = "l", lwd = 3, col = "black",
     ylim = range(c(hi, lo)),
     xlab = "Years since first visit", ylab = "log(ADAS11 + 1)",
     main = "Mode of variation: mean +/- 2 SD x FPC1")
lines(wg, hi, col = "#BC4749", lwd = 2, lty = 2)
lines(wg, lo, col = "#2A6F97", lwd = 2, lty = 2)
grid(col = "grey90")
legend("topleft", bty = "n", cex = 0.85,
       legend = c("mean", "mean + 2SD x FPC1 (high score)",
                  "mean - 2SD x FPC1 (low score)"),
       col = c("black", "#BC4749", "#2A6F97"), lwd = 2, lty = c(1, 2, 2))
par(op)
dev.off()
cat("saved adni_fig_mode.png\n")

## ===========================================================================
## FIG 3: distribution of FPC1 scores  ->  adni_fig_scores.png
## ===========================================================================
xi1 <- fit$xiEst[, 1]
png("adni_fig_scores.png", width = 7, height = 4.6, units = "in", res = 300)
op <- par(mar = c(4.5, 4.5, 3, 1))
hist(xi1, breaks = 30, col = "#9dc3d4", border = "white",
     xlab = "FPC1 score", main = "Distribution of FPC1 scores across patients")
abline(v = 0, col = "#BC4749", lwd = 2)
par(op)
dev.off()
cat("saved adni_fig_scores.png\n")

cat(sprintf("\nFPC1 score range %.2f to %.2f (SD %.2f)\n",
            min(xi1), max(xi1), sd(xi1)))
cat("All three log-scale figures written for Chapter 4.\n")
###############################################################################
