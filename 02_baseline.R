###############################################################################
## 02_baseline.R
## ---------------------------------------------------------------------------
## No-missingness baseline: FPCA on the COMPLETE data (before any dropout),
## over 100 replications. Gives the reference point against which every
## dropout scenario is judged in Chapter 3.
##
## Requires the helper functions and constants from 01_simulation_main.R,
## so source that script first:
##     source("01_simulation_main.R")
## (that defines gen_rep, trapz, GRID1, GRID2, mu_true, phi_true, traj1,
##  FVE, SEED0 and the constants they use).
##
## Output: mc_baseline_nomissing.rds ; prints the baseline MISE per dataset.
###############################################################################

library(fdapace)
 
## ---- constants (identical to 01_simulation_main.R) ----
N <- 700; TMAX <- 12; EVERY <- 2; JIT <- 1; FVE <- 0.99; SEED0 <- 20260000
GRID1 <- seq(0, TMAX, length.out = 101)
GRID2 <- seq(0, 1,    length.out = 101)
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
trapz <- function(x, y){o<-order(x);x<-x[o];y<-y[o]
  if(length(x)<2)return(NA);sum(diff(x)*(head(y,-1)+tail(y,-1))/2)}
 
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
 
analyse_complete <- function(df, dataset, rep_obj) {
  ids <- sort(unique(df$id))
  Ly  <- lapply(ids, function(i) df$y[df$id == i])
  Lt  <- lapply(ids, function(i) df$time[df$id == i])
  keep <- sapply(Ly, length) >= 2
  Ly <- Ly[keep]; Lt <- Lt[keep]; ids_fit <- ids[keep]
  if (length(Ly) < 10) return(NULL)
  fit <- tryCatch(FPCA(Ly, Lt, optns=list(dataType="Sparse",
                       FVEthreshold=FVE, verbose=FALSE)), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  wg <- fit$workGrid; yhat <- fitted(fit)
  grid <- if (dataset=="DS1") GRID1 else GRID2
  ise_i <- numeric(length(ids_fit))
  for (k in seq_along(ids_fit)) {
    i <- ids_fit[k]
    fhat <- approx(wg, yhat[k,], xout=grid, rule=2)$y
    ftru <- if (dataset=="DS1") traj1(grid, rep_obj$b[[i]])
            else mu_true(grid) + rep_obj$xi[i]*phi_true(grid)
    ise_i[k] <- trapz(grid, (fhat-ftru)^2)
  }
  res <- data.frame(dataset=dataset, ise_recon=mean(ise_i, na.rm=TRUE),
                    ise_mu=NA_real_, ise_phi=NA_real_)
  if (dataset=="DS2") {
    mh <- approx(wg, fit$mu, xout=GRID2, rule=2)$y
    res$ise_mu <- trapz(GRID2, (mh-mu_true(GRID2))^2)
    ph <- approx(wg, fit$phi[,1], xout=GRID2, rule=2)$y
    pt <- phi_true(GRID2); if (sum(ph*pt)<0) ph <- -ph
    nrm <- sqrt(trapz(GRID2, ph^2)); if (is.finite(nrm)&&nrm>0) ph <- ph/nrm
    res$ise_phi <- trapz(GRID2, (ph-pt)^2)
  }
  res
}
 
## run replications until BOTH datasets have 100 clean baseline values.
## DS2 never explodes, so its first 100 reps are kept. DS1 can explode on
## extreme random-effect draws, so extra reps (new seeds) are generated until
## 100 clean DS1 values are collected. This matches the 100-replication design
## used everywhere else.
ds1_rows <- list(); ds2_rows <- list()
n_ds1 <- 0; n_ds2 <- 0; r <- 0; tries <- 0
while ((n_ds1 < 100 || n_ds2 < 100) && tries < 400) {
  r <- r + 1; tries <- tries + 1
  ro <- tryCatch(gen_rep(SEED0 + r), error=function(e) NULL)
  if (is.null(ro)) next
  ## DS1
  if (n_ds1 < 100) {
    res1 <- tryCatch(analyse_complete(ro$ds1, "DS1", ro), error=function(e) NULL)
    if (!is.null(res1) && is.finite(res1$ise_recon) && res1$ise_recon <= 500) {
      res1$rep <- r; ds1_rows[[length(ds1_rows)+1]] <- res1; n_ds1 <- n_ds1 + 1
    }
  }
  ## DS2
  if (n_ds2 < 100) {
    res2 <- tryCatch(analyse_complete(ro$ds2, "DS2", ro), error=function(e) NULL)
    if (!is.null(res2)) { res2$rep <- r; ds2_rows[[length(ds2_rows)+1]] <- res2; n_ds2 <- n_ds2 + 1 }
  }
  if (tries %% 20 == 0) cat(sprintf("  tries %d: DS1 clean=%d, DS2=%d\n", tries, n_ds1, n_ds2))
}
raw_base <- rbind(do.call(rbind, ds1_rows), do.call(rbind, ds2_rows))
cat(sprintf("Collected DS1 n=%d, DS2 n=%d (after %d replications).\n",
            n_ds1, n_ds2, tries))
 
## aggregate
agg <- function(x){x<-x[is.finite(x)];if(!length(x))return(c(NA,NA));c(mean(x),sd(x)/sqrt(length(x)))}
for (ds in c("DS1","DS2")) {
  d <- raw_base[raw_base$dataset==ds,]; a <- agg(d$ise_recon)
  cat(sprintf("\n[%s] no-missing baseline (n=%d):\n", ds, nrow(d)))
  cat(sprintf("  reconstruction MISE = %.3f (SE %.3f)\n", a[1], a[2]))
  if (ds=="DS2") {
    cat(sprintf("  MISE_mu  = %.4f\n", agg(d$ise_mu)[1]))
    cat(sprintf("  MISE_phi = %.4f\n", agg(d$ise_phi)[1]))
  }
}
saveRDS(raw_base, "mc_baseline_nomissing.rds")
cat("\nSaved: mc_baseline_nomissing.rds\n")
###############################################################################
