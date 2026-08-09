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

## check the required objects are in memory
need <- c("gen_rep", "trapz", "GRID1", "GRID2", "mu_true", "phi_true", "traj1",
          "FVE", "SEED0")
missing <- need[!sapply(need, exists)]
if (length(missing))
  stop("Missing objects: ", paste(missing, collapse = ", "),
       "\n-> source('01_simulation_main.R') first.")

analyse_complete <- function(df, dataset, rep_obj) {
  ids <- sort(unique(df$id))
  Ly  <- lapply(ids, function(i) df$y[df$id == i])
  Lt  <- lapply(ids, function(i) df$time[df$id == i])
  keep <- sapply(Ly, length) >= 2
  Ly <- Ly[keep]; Lt <- Lt[keep]; ids_fit <- ids[keep]
  if (length(Ly) < 10) return(NULL)

  fit <- tryCatch(FPCA(Ly, Lt, optns = list(dataType = "Sparse",
                       FVEthreshold = FVE, verbose = FALSE)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  wg <- fit$workGrid; yhat <- fitted(fit)
  grid <- if (dataset == "DS1") GRID1 else GRID2

  ise_i <- numeric(length(ids_fit))
  for (k in seq_along(ids_fit)) {
    i <- ids_fit[k]
    fhat <- approx(wg, yhat[k, ], xout = grid, rule = 2)$y
    ftru <- if (dataset == "DS1") traj1(grid, rep_obj$b[[i]])
            else mu_true(grid) + rep_obj$xi[i] * phi_true(grid)
    ise_i[k] <- trapz(grid, (fhat - ftru)^2)
  }
  res <- data.frame(dataset = dataset,
                    ise_recon = mean(ise_i, na.rm = TRUE),
                    ise_mu = NA_real_, ise_phi = NA_real_)
  if (dataset == "DS2") {
    mh <- approx(wg, fit$mu, xout = GRID2, rule = 2)$y
    res$ise_mu <- trapz(GRID2, (mh - mu_true(GRID2))^2)
    ph <- approx(wg, fit$phi[, 1], xout = GRID2, rule = 2)$y
    pt <- phi_true(GRID2)
    if (sum(ph * pt) < 0) ph <- -ph
    nrm <- sqrt(trapz(GRID2, ph^2)); if (is.finite(nrm) && nrm > 0) ph <- ph / nrm
    res$ise_phi <- trapz(GRID2, (ph - pt)^2)
  }
  res
}

## run 100 replications, serially
rows <- list()
for (r in 1:100) {
  ro <- tryCatch(gen_rep(SEED0 + r), error = function(e) NULL)
  if (is.null(ro)) next
  for (spec in list(list(ro$ds1, "DS1"), list(ro$ds2, "DS2"))) {
    res <- tryCatch(analyse_complete(spec[[1]], spec[[2]], ro),
                    error = function(e) NULL)
    if (!is.null(res)) { res$rep <- r; rows[[length(rows)+1]] <- res }
  }
  if (r %% 20 == 0) cat("  done", r, "replications\n")
}
raw_base <- do.call(rbind, rows)

## drop DS1 explosions (extreme random-effect draws)
bad <- raw_base$dataset == "DS1" &
       (!is.finite(raw_base$ise_recon) | raw_base$ise_recon > 500)
bad_reps <- unique(raw_base$rep[bad])
raw_base <- raw_base[!(raw_base$rep %in% bad_reps), ]
cat(sprintf("Dropped %d exploded replications; %d remain.\n",
            length(bad_reps), length(unique(raw_base$rep))))

## aggregate
agg <- function(x) { x <- x[is.finite(x)]
  if (length(x) == 0) return(c(NA, NA)); c(mean(x), sd(x)/sqrt(length(x))) }
for (ds in c("DS1", "DS2")) {
  d <- raw_base[raw_base$dataset == ds, ]; a <- agg(d$ise_recon)
  cat(sprintf("\n[%s] no-missing baseline (n=%d):\n", ds, nrow(d)))
  cat(sprintf("  reconstruction MISE = %.3f (SE %.3f)\n", a[1], a[2]))
  if (ds == "DS2") {
    cat(sprintf("  MISE_mu  = %.4f\n", agg(d$ise_mu)[1]))
    cat(sprintf("  MISE_phi = %.4f\n", agg(d$ise_phi)[1]))
  }
}
saveRDS(raw_base, "mc_baseline_nomissing.rds")
cat("\nSaved: mc_baseline_nomissing.rds\n")
###############################################################################
