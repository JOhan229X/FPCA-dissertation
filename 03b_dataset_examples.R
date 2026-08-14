###############################################################################
## 03b_dataset_examples.R  (self-contained)
## ---------------------------------------------------------------------------
## Spaghetti plot per dataset: 100 random individual trajectories, drawn as
## semi-transparent connected lines, showing the two data-generating models
## before any dropout.
##
## Self-contained: defines its own constants and gen_rep (identical to
## 01_simulation_main.R), so it runs on its own WITHOUT sourcing 01 and WITHOUT
## re-running the main simulation.
##
## Output: fig_dataset1_example.png, fig_dataset2_example.png
###############################################################################

## ---- constants (identical to 01_simulation_main.R) ----
N <- 700; TMAX <- 12; EVERY <- 2; JIT <- 1; SEED0 <- 20260000
SD_B <- sqrt(c(3,5,2,2)); SD_EPS1 <- 3
SD_XI <- sqrt(2);         SD_EPS2 <- 1
traj1 <- function(t, b) {
  y0 <- 10+b[1]; yinf <- 35+b[2]; tau <- 10+b[3]
  alpha <- 5+b[4]; theta <- 1+b[4]/3
  ts <- ifelse(t <= 0, 1e-8, t)
  yinf + (y0-yinf)/(1+exp(alpha*(log(ts)-log(tau))))^theta
}
mu_true  <- function(t) -0.5 + 1.5*sin(10*pi*(t-0.5)) + 4*(t-2)^3
phi_true <- function(t) sqrt(2)*cos(4*pi*t)

## ---- gen_rep: DS2 times generated INDEPENDENTLY on [0,1] (matches 01) ----
gen_rep <- function(seed) {
  set.seed(seed)
  theo <- seq(0, TMAX, by = EVERY)
  theo2 <- seq(0, 1, length.out = length(theo)); JIT2 <- JIT/TMAX
  tlist <- vector("list", N); blist <- vector("list", N); t2list <- vector("list", N)
  for (i in seq_len(N)) {
    blist[[i]] <- rnorm(4, 0, SD_B)
    tij <- theo + runif(length(theo), -JIT, JIT)
    tlist[[i]] <- sort(pmin(pmax(tij, 0), TMAX))
    t2 <- theo2 + runif(length(theo2), -JIT2, JIT2)
    t2list[[i]] <- sort(pmin(pmax(t2, 0), 1))
  }
  ds1 <- do.call(rbind, lapply(seq_len(N), function(i){
    t <- tlist[[i]]; ys <- traj1(t, blist[[i]])
    data.frame(id=i, time=t, truth=ys, y=ys+rnorm(length(t),0,SD_EPS1))}))
  xis <- rnorm(N, 0, SD_XI)
  ds2 <- do.call(rbind, lapply(seq_len(N), function(i){
    t01 <- t2list[[i]]; ft <- mu_true(t01)+xis[i]*phi_true(t01)
    data.frame(id=i, time=t01, truth=ft, y=ft+rnorm(length(t01),0,SD_EPS2))}))
  list(ds1=ds1, ds2=ds2, b=blist, xi=xis)
}

ro <- gen_rep(SEED0 + 1)

## ==========================================
## Dataset 1: 100 random individual trajectories (spaghetti)
## ==========================================
png("fig_dataset1_example.png", width = 7.5, height = 5, units = "in", res = 300)
op <- par(mar = c(4, 4.2, 3, 1))
d1 <- ro$ds1
plot(NULL, xlim = c(0, TMAX), ylim = range(d1$y),
     xlab = "t", ylab = "Simulated marker",
     main = "Dataset 1: 100 random individuals")
set.seed(1)
for (i in sample(unique(d1$id), 100)) {
  d <- d1[d1$id == i, ]
  lines(d$time, d$y, col = adjustcolor("black", 0.25), lwd = 0.7)
}
par(op); dev.off(); cat("saved fig_dataset1_example.png\n")

## ==========================================
## Dataset 2: 100 random individual trajectories (spaghetti)
## ==========================================
png("fig_dataset2_example.png", width = 7.5, height = 5, units = "in", res = 300)
op <- par(mar = c(4, 4.2, 3, 1))
d2 <- ro$ds2
plot(NULL, xlim = c(0, 1), ylim = range(d2$y),
     xlab = "t (rescaled to [0,1])", ylab = "Simulated marker",
     main = "Dataset 2: 100 random individuals")
set.seed(1)
for (i in sample(unique(d2$id), 100)) {
  d <- d2[d2$id == i, ]
  lines(d$time, d$y, col = adjustcolor("black", 0.25), lwd = 0.7)
}
## overlay the population mean so the single-component structure is visible
tg2 <- seq(0, 1, length.out = 200)
lines(tg2, mu_true(tg2), col = "#BC4749", lwd = 2.5)
legend("topleft", bty = "n", cex = 0.8,
       legend = c("individuals", "population mean"),
       col = c(adjustcolor("black",0.5), "#BC4749"), lwd = c(0.7, 2.5))
par(op); dev.off(); cat("saved fig_dataset2_example.png\n")

cat("\nBoth dataset spaghetti figures written.\n")
###############################################################################
