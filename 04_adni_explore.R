###############################################################################
## 04_adni_explore.R
## ---------------------------------------------------------------------------
## ADNI data preparation and missingness description for Chapter 4:
##   - reads MSc_Data.csv (RID, EXAMDATE, ADAS11)
##   - classifies each patient's missingness pattern (none / monotone /
##     intermittent)
##   - cleans the data (drop missing-ADAS11 rows and single-visit patients)
##   - saves the cleaned data as adni_clean.rds for the FPCA scripts
##   - draws the spaghetti plot of individual trajectories
##
## Output: adni_clean.rds, adni_fig_spaghetti.png ; prints the pattern table.
###############################################################################

raw <- read.csv("MSc_Data.csv", stringsAsFactors = FALSE)
raw$date   <- as.Date(raw$EXAMDATE, format = "%d/%m/%Y")
raw$ADAS11 <- suppressWarnings(as.numeric(raw$ADAS11))
raw <- raw[order(raw$RID, raw$date), ]

cat(sprintf("Raw: %d records, %d patients, %.1f%% ADAS11 missing\n",
            nrow(raw), length(unique(raw$RID)), 100*mean(is.na(raw$ADAS11))))

## ---- classify missingness pattern per patient ----
pattern <- sapply(split(raw, raw$RID), function(d) {
  m <- is.na(d$ADAS11[order(d$date)])
  if (!any(m)) return("none")
  ## monotone if all missing visits are at the end (once missing, stays missing)
  if (all(diff(which(m)) == 1) && max(which(m)) == length(m) &&
      min(which(m)) == length(m) - sum(m) + 1) return("monotone")
  "intermittent"
})
tab <- table(factor(pattern, levels = c("none","monotone","intermittent")))
cat("\nMissingness pattern across patients:\n")
print(tab)
cat(sprintf("  (%.0f%% none, %.0f%% monotone, %.0f%% intermittent)\n",
            100*tab[1]/sum(tab), 100*tab[2]/sum(tab), 100*tab[3]/sum(tab)))

## ---- informative check: mean ADAS11 just before a missing vs observed visit ----
by <- split(raw, raw$RID); prev_miss <- c(); prev_obs <- c()
for (d in by) {
  d <- d[order(d$date), ]
  for (j in 2:nrow(d)) {
    if (is.na(d$ADAS11[j-1])) next
    if (is.na(d$ADAS11[j])) prev_miss <- c(prev_miss, d$ADAS11[j-1])
    else                    prev_obs  <- c(prev_obs,  d$ADAS11[j-1])
  }
}
cat(sprintf("\nMean previous ADAS11 before a MISSING visit:  %.2f\n", mean(prev_miss)))
cat(sprintf("Mean previous ADAS11 before an OBSERVED visit: %.2f\n", mean(prev_obs)))

## ---- clean: drop missing-ADAS11 rows and single-visit patients ----
dat <- raw[!is.na(raw$ADAS11), c("RID", "date", "ADAS11")]
dat <- dat[dat$RID %in% names(which(table(dat$RID) >= 2)), ]
dat <- dat[order(dat$RID, dat$date), ]
dat$years <- ave(as.numeric(dat$date), dat$RID,
                 FUN = function(x) (x - min(x)) / 365.25)
saveRDS(dat, "adni_clean.rds")
cat(sprintf("\nCleaned: %d patients, %d records, median %d visits\n",
            length(unique(dat$RID)), nrow(dat),
            median(table(dat$RID))))

## ---- spaghetti plot ----
png("adni_fig_spaghetti.png", width = 8, height = 5, units = "in", res = 300)
op <- par(mar = c(4.5, 4.5, 2, 1))
ids <- unique(dat$RID)
set.seed(1); show <- sample(ids, min(120, length(ids)))   # sample for clarity
plot(NULL, xlim = range(dat$years), ylim = range(dat$ADAS11),
     xlab = "Years since first visit", ylab = "ADAS11",
     main = "ADNI ADAS11 trajectories")
for (i in show) {
  d <- dat[dat$RID == i, ]
  lines(d$years, d$ADAS11, col = adjustcolor("grey50", 0.35))
}
## smoothed population mean
lo <- loess(ADAS11 ~ years, data = dat, span = 0.5)
xg <- seq(min(dat$years), max(dat$years), length.out = 100)
lines(xg, predict(lo, xg), col = "#BC4749", lwd = 3)
legend("topleft", bty = "n", legend = c("individual patients", "population mean"),
       col = c("grey50", "#BC4749"), lwd = c(1, 3))
par(op); dev.off()
cat("saved adni_fig_spaghetti.png\n")
###############################################################################
