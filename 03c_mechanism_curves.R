###############################################################################
## 03c_mechanism_curves.R  
## ---------------------------------------------------------------------------
## Visualises the dropout mechanisms: probability of dropout vs its driver.
## Left: the time-based driver (MCAR) on the Dataset 1 domain.
## Right: the two value-based drivers on their actual value range, so the
## contrast between the smooth "increasing" driver and the "threshold" step is
## visible on a common scale.
## Illustrative beta0; beta1 as in the study (0.30 MCAR, 0.15 increasing, 2.0 threshold).
## Output: fig_mechanism_curves.png
###############################################################################

png("fig_mechanism_curves.png", width = 2400, height = 1000, res = 220)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

expit <- function(z) 1 / (1 + exp(-z))

## Left: time-based driver (MCAR), domain of Dataset 1
tt <- seq(0, 12, length.out = 400)
plot(tt, expit(-2.5 + 0.30 * tt), type = "l", lwd = 2.5, col = "seagreen",
     ylim = c(0, 1), xlab = "Visit time", ylab = "P(dropout at visit)",
     main = "Time-based driver (MCAR)")
grid(lty = 3, col = "grey80")

## Right: value-based drivers on their actual range
vv  <- seq(0, 40, length.out = 400)
med <- 20
plot(vv, expit(-4.0 + 0.15 * vv), type = "l", lwd = 2.5, col = "steelblue",
     ylim = c(0, 1), xlab = "Previous or current value",
     ylab = "P(dropout at visit)", main = "Value-based drivers")
lines(vv, expit(-1.5 + 2.0 * (vv > med)), lwd = 2.5, col = "firebrick")
abline(v = med, lty = 3, col = "grey50")
grid(lty = 3, col = "grey80")
legend("topleft", bty = "n", lwd = 2.5,
       col = c("steelblue", "firebrick"),
       legend = c("Increasing (driver: value)", "Threshold (step at median)"))

dev.off()
cat("saved fig_mechanism_curves.png\n")
###############################################################################
