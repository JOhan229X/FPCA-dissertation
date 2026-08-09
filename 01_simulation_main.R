###############################################################################
## 01_simulation_main.R
## ---------------------------------------------------------------------------
## Monte Carlo simulation study for the dissertation
##   "Robustness of FPCA to Informative Missingness" (MATH5872M).
##
## Generates two longitudinal datasets, imposes six dropout mechanisms
## (MCAR / MAR / MNAR) at two dropout rates, fits FPCA to each, and records
## the integrated squared error (ISE) against the known truth. Averaging the
## ISE over R replications gives the mean integrated squared error (MISE).
##
## Reported quantities (averaged over replications):
##   both datasets : MISE of the reconstruction (overall / observed / missing)
##   Dataset 2     : MISE of the estimated mean mu-hat and eigenfunction phi-hat_1
##
## Output files: mc_raw_replications.rds  (one row per rep x scenario)
##               mc_summary.csv / .rds    (aggregated MISE per scenario)
##
## Reproducible: replication r uses seed SEED0 + r; calibration uses seed 2024.
###############################################################################

library(fdapace)
library(parallel)

## ===========================================================================
## 0. Settings
## ===========================================================================
PILOT   <- FALSE                 # TRUE = 3 reps (timing); FALSE = 100 reps
R_REPS  <- if (PILOT) 3 else 100
N_CORES <- max(1, detectCores() - 1)

N     <- 700                     # number of individuals
TMAX  <- 12                      # Dataset 1 time domain [0, 12]
EVERY <- 2                       # nominal visit spacing (years)
JIT   <- 1                       # uniform jitter on visit times
FVE   <- 0.99                    # fraction-of-variance-explained threshold
SEED0 <- 20260000                # replication r uses seed SEED0 + r

GRID1 <- seq(0, TMAX, length.out = 101)   # DS1 evaluation grid [0, 12]
GRID2 <- seq(0, 1,    length.out = 101)   # DS2 evaluation grid [0, 1]

## ---- truth: Dataset 1 (five-parameter logistic mixed model) ----
SD_B    <- sqrt(c(3, 5, 2, 2))   # random-effect SDs (B = diag(3,5,2,2) is a VARIANCE)
SD_EPS1 <- 3                     # measurement-error SD
traj1 <- function(t, b) {
  y0 <- 10 + b[1]; yinf <- 35 + b[2]; tau <- 10 + b[3]
  alpha <- 5 + b[4]; theta <- 1 + b[4] / 3
  ts <- ifelse(t <= 0, 1e-8, t)
  yinf + (y0 - yinf) / (1 + exp(alpha * (log(ts) - log(tau))))^theta
}

## ---- truth: Dataset 2 (single-component Karhunen-Loeve) ----
SD_XI   <- sqrt(2)               # score SD, xi ~ N(0, 2)
SD_EPS2 <- 1                     # measurement-error SD
mu_true  <- function(t) -0.5 + 1.5 * sin(10 * pi * (t - 0.5)) + 4 * (t - 2)^3
phi_true <- function(t) sqrt(2) * cos(4 * pi * t)

## six mechanisms: slope beta1 for each (mcar uses time; others use value drivers)
MECHS   <- c(mcar = 0.30, thr_mar = 2.0, inc_mar = 0.15,
             thr_mnar = 2.0, inc_mnar = 0.15)
TARGETS <- c(0.30, 0.60)         # target overall dropout rates

inv_logit <- function(z) 1 / (1 + exp(-z))
trapz <- function(x, y) { o <- order(x); x <- x[o]; y <- y[o]
  if (length(x) < 2) return(NA); sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) }

## ===========================================================================
## 1. Generate one replication (both datasets share the visit times)
## ===========================================================================
gen_rep <- function(seed) {
  set.seed(seed)
  theo <- seq(0, TMAX, by = EVERY)
  tlist <- vector("list", N); blist <- vector("list", N)
  for (i in seq_len(N)) {
    blist[[i]] <- rnorm(4, 0, SD_B)
    tij <- theo + runif(length(theo), -JIT, JIT)
    tlist[[i]] <- sort(pmin(pmax(tij, 0), TMAX))
  }
  ds1 <- do.call(rbind, lapply(seq_len(N), function(i) {
    t <- tlist[[i]]; ys <- traj1(t, blist[[i]])
    data.frame(id = i, time = t, truth = ys,
               y = ys + rnorm(length(t), 0, SD_EPS1)) }))
  xis <- rnorm(N, 0, SD_XI)
  ds2 <- do.call(rbind, lapply(seq_len(N), function(i) {
    t01 <- tlist[[i]] / TMAX
    ft  <- mu_true(t01) + xis[i] * phi_true(t01)
    data.frame(id = i, time = t01, truth = ft,
               y = ft + rnorm(length(t01), 0, SD_EPS2)) }))
  list(ds1 = ds1, ds2 = ds2, b = blist, xi = xis)
}

## ===========================================================================
## 2. Impose monotone dropout (beta0 supplied; not recalibrated per replication)
## ===========================================================================
apply_dropout <- function(df, mech, b0, b1, nu) {
  df <- df[order(df$id, df$time), ]
  obs <- logical(nrow(df))
  for (id in unique(df$id)) {
    idx <- which(df$id == id); t <- df$time[idx]; y <- df$y[idx]; n <- length(idx)
    keep <- rep(TRUE, n)
    for (j in seq_len(n)) {
      pr <- switch(mech,
                   mcar      = inv_logit(b0 + b1 * t[j]),
                   fixed_mar = if (j >= 2 && y[j-1] > nu) 1 else 0,
                   thr_mar   = inv_logit(b0 + b1 * (if (j >= 2) as.numeric(y[j-1] > nu) else 0)),
                   inc_mar   = inv_logit(b0 + b1 * (if (j >= 2) y[j-1] else y[1])),
                   thr_mnar  = inv_logit(b0 + b1 * as.numeric(y[j] > nu)),
                   inc_mnar  = inv_logit(b0 + b1 * y[j]))
      if (runif(1) < pr) { keep[j:n] <- FALSE; break }   # monotone: drop rest
    }
    obs[idx] <- keep
  }
  df$observed <- obs; df
}

## ===========================================================================
## 3. One-off calibration of beta0 on a reference replication (seed 2024)
##    beta0 is chosen so the realised dropout rate hits the target, then held
##    fixed across all replications.
## ===========================================================================
cat("Calibrating beta0 once on a reference replication ...\n")
ref <- gen_rep(2024)
calibrate <- function(df, mech, b1, nu, target) {
  f <- function(b0) { set.seed(123)
    1 - mean(apply_dropout(df, mech, b0, b1, nu)$observed) - target }
  out <- tryCatch(uniroot(f, lower = -50, upper = 50, tol = 1e-3),
                  error = function(e) NULL)
  if (is.null(out)) NA else out$root
}
build_scen <- function(df) {
  nu <- median(df$y); rows <- list()
  for (m in names(MECHS)) for (tg in TARGETS) {
    b0 <- calibrate(df, m, MECHS[[m]], nu, tg)
    rows[[length(rows)+1]] <- data.frame(
      scenario = sprintf("%s_%d", m, round(tg*100)),
      mech = m, b1 = MECHS[[m]], b0 = b0, nu = nu, stringsAsFactors = FALSE) }
  rows[[length(rows)+1]] <- data.frame(scenario = "fixed_mar_det",
    mech = "fixed_mar", b1 = 0, b0 = 0, nu = nu, stringsAsFactors = FALSE)
  do.call(rbind, rows)
}
SCEN1 <- build_scen(ref$ds1)
SCEN2 <- build_scen(ref$ds2)
cat(sprintf("  done. nu(DS1)=%.2f  nu(DS2)=%.2f\n\n", SCEN1$nu[1], SCEN2$nu[1]))

## ===========================================================================
## 4. Analyse one dataset x one scenario -> one row of ISE metrics
## ===========================================================================
analyse_scenario <- function(df, sc, dataset, rep_obj) {
  scen <- apply_dropout(df, sc$mech, sc$b0, sc$b1, sc$nu)
  rate <- 1 - mean(scen$observed)

  o <- scen[scen$observed, ]; ids <- sort(unique(o$id))
  Ly <- lapply(ids, function(i) o$y[o$id == i])
  Lt <- lapply(ids, function(i) o$time[o$id == i])
  keep <- sapply(Ly, length) >= 2
  Ly <- Ly[keep]; Lt <- Lt[keep]; ids_fit <- ids[keep]
  if (length(Ly) < 10) return(NULL)

  fit <- tryCatch(FPCA(Ly, Lt, optns = list(dataType = "Sparse",
                                            FVEthreshold = FVE, verbose = FALSE)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  wg <- fit$workGrid; yhat <- fitted(fit)
  grid <- if (dataset == "DS1") GRID1 else GRID2

  ## reconstruction ISE per person on the common grid, then averaged
  ise_i <- numeric(length(ids_fit))
  for (k in seq_along(ids_fit)) {
    i <- ids_fit[k]
    fhat <- approx(wg, yhat[k, ], xout = grid, rule = 2)$y
    ftru <- if (dataset == "DS1") traj1(grid, rep_obj$b[[i]])
            else mu_true(grid) + rep_obj$xi[i] * phi_true(grid)
    ise_i[k] <- trapz(grid, (fhat - ftru)^2)
  }

  ## squared error split into observed / missing portions
  se_o <- c(); se_m <- c()
  for (k in seq_along(ids_fit)) {
    i <- ids_fit[k]; di <- scen[scen$id == i, ]
    pr <- approx(wg, yhat[k, ], xout = di$time, rule = 2)$y
    io <- di$observed
    if (sum(io)  >= 2) se_o <- c(se_o, trapz(di$time[io],  (pr[io]  - di$truth[io])^2))
    if (sum(!io) >= 2) se_m <- c(se_m, trapz(di$time[!io], (pr[!io] - di$truth[!io])^2))
  }

  res <- data.frame(dataset = dataset, scenario = sc$scenario, rate = rate,
                    K = fit$selectK, n_fit = length(ids_fit),
                    ise_recon = mean(ise_i, na.rm = TRUE),
                    ise_obs = mean(se_o, na.rm = TRUE),
                    ise_miss = if (length(se_m)) mean(se_m, na.rm = TRUE) else NA,
                    ise_mu = NA, ise_phi = NA)

  ## Dataset 2: mean and eigenfunction errors (truth known)
  if (dataset == "DS2") {
    mh <- approx(wg, fit$mu, xout = GRID2, rule = 2)$y
    res$ise_mu <- trapz(GRID2, (mh - mu_true(GRID2))^2)
    ph <- approx(wg, fit$phi[, 1], xout = GRID2, rule = 2)$y
    pt <- phi_true(GRID2)
    if (sum(ph * pt) < 0) ph <- -ph                 # sign-align to truth
    nrm <- sqrt(trapz(GRID2, ph^2)); if (is.finite(nrm) && nrm > 0) ph <- ph / nrm
    res$ise_phi <- trapz(GRID2, (ph - pt)^2)
  }
  res
}

## ===========================================================================
## 5. One full replication = both datasets x all scenarios
## ===========================================================================
run_rep <- function(r) {
  ro <- gen_rep(SEED0 + r); out <- list()
  for (s in seq_len(nrow(SCEN1)))
    out[[length(out)+1]] <- analyse_scenario(ro$ds1, SCEN1[s, ], "DS1", ro)
  for (s in seq_len(nrow(SCEN2)))
    out[[length(out)+1]] <- analyse_scenario(ro$ds2, SCEN2[s, ], "DS2", ro)
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(NULL)
  cbind(rep = r, do.call(rbind, out))
}

## ===========================================================================
## 6. Run all replications in parallel
## ===========================================================================
cat(sprintf("Running %d replications on %d cores ...\n", R_REPS, N_CORES))
t0 <- Sys.time()
cl <- makeCluster(N_CORES)
clusterEvalQ(cl, library(fdapace))
clusterExport(cl, c("N","TMAX","EVERY","JIT","FVE","SEED0","GRID1","GRID2",
                    "SD_B","SD_EPS1","traj1","SD_XI","SD_EPS2","mu_true","phi_true",
                    "MECHS","TARGETS","inv_logit","trapz","gen_rep","apply_dropout",
                    "analyse_scenario","SCEN1","SCEN2","run_rep"))
reps <- parLapply(cl, seq_len(R_REPS),
                  function(r) tryCatch(run_rep(r), error = function(e) NULL))
stopCluster(cl)
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("Finished in %.1f min (%.2f min/rep).\n", el, el / R_REPS))

raw <- do.call(rbind, Filter(Negate(is.null), reps))
saveRDS(raw, "mc_raw_replications.rds")
cat(sprintf("Collected %d rows from %d replications.\n",
            nrow(raw), length(unique(raw$rep))))

## ===========================================================================
## 7. Aggregate ISE -> MISE (mean over reps) + Monte Carlo SE
## ===========================================================================
agg <- function(x) c(mean(x, na.rm = TRUE),
                     sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
summ <- do.call(rbind, lapply(
  split(raw, list(raw$dataset, raw$scenario), drop = TRUE), function(d) {
    a_rec <- agg(d$ise_recon); a_o <- agg(d$ise_obs); a_m <- agg(d$ise_miss)
    out <- data.frame(dataset = d$dataset[1], scenario = d$scenario[1],
                      n_rep = nrow(d), rate_mean = mean(d$rate), rate_sd = sd(d$rate),
                      K_median = median(d$K),
                      MISE_recon = a_rec[1], MISE_recon_se = a_rec[2],
                      MISE_obs = a_o[1], MISE_miss = a_m[1], MISE_miss_se = a_m[2],
                      MISE_mu = NA, MISE_mu_se = NA, MISE_phi = NA, MISE_phi_se = NA)
    if (d$dataset[1] == "DS2") {
      a_mu <- agg(d$ise_mu); a_ph <- agg(d$ise_phi)
      out$MISE_mu <- a_mu[1]; out$MISE_mu_se <- a_mu[2]
      out$MISE_phi <- a_ph[1]; out$MISE_phi_se <- a_ph[2]
    }
    out
  }))
rownames(summ) <- NULL

cat("\n\n=========== RECONSTRUCTION MISE (both datasets) ===========\n\n")
print(summ[, c("dataset","scenario","n_rep","rate_mean","rate_sd","K_median",
               "MISE_recon","MISE_recon_se","MISE_obs","MISE_miss")],
      row.names = FALSE, digits = 3)
cat("\n\n=========== COMPONENTS MISE (Dataset 2) ===========\n\n")
d2 <- summ[summ$dataset == "DS2", ]
print(d2[, c("scenario","MISE_mu","MISE_mu_se","MISE_phi","MISE_phi_se")],
      row.names = FALSE, digits = 3)

write.csv(summ, "mc_summary.csv", row.names = FALSE)
saveRDS(summ, "mc_summary.rds")
cat("\nSaved: mc_raw_replications.rds, mc_summary.csv / .rds\n")
###############################################################################
