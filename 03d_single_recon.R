###############################################################################
## 03d_single_recon.R
## ---------------------------------------------------------------------------
## One individual from Dataset 1, true curve vs FPCA reconstruction, under
## increasing MAR and increasing MNAR at 60%. Same individual in both panels,
## so the only difference is the mechanism.
##
## Requires the objects from 01_simulation_main.R in memory:
##     source("01_simulation_main.R")
## (needs gen_rep, apply_dropout, traj1, SCEN1, GRID1, FVE, SEED0).
##
## IMPORTANT: run it, LOOK at the figure, and keep it only if the MNAR panel
## shows a visible departure. If no individual is dramatic, drop the figure.
##
## Output: fig_single_recon.png
###############################################################################

library(fdapace)
need <- c("gen_rep","apply_dropout","traj1","SCEN1","GRID1","FVE","SEED0")
if (any(!sapply(need, exists)))
  stop("source('01_simulation_main.R') first.")

## calibrated intercepts for the two 60% mechanisms
b0_mar  <- SCEN1$b0[SCEN1$scenario == "inc_mar_60"]
b1_mar  <- SCEN1$b1[SCEN1$scenario == "inc_mar_60"]
nu_mar  <- SCEN1$nu[SCEN1$scenario == "inc_mar_60"]
b0_mnar <- SCEN1$b0[SCEN1$scenario == "inc_mnar_60"]
b1_mnar <- SCEN1$b1[SCEN1$scenario == "inc_mnar_60"]
nu_mnar <- SCEN1$nu[SCEN1$scenario == "inc_mnar_60"]

## one replication of the full (complete) Dataset 1
## If no clear individual is found, bump SEED_OFFSET to 2, 3, ... and rerun.
SEED_OFFSET <- 1
ro  <- gen_rep(SEED0 + SEED_OFFSET)
ds1 <- ro$ds1

## apply the two mechanisms
d_mar  <- apply_dropout(ds1, "inc_mar",  b0_mar,  b1_mar,  nu_mar)
d_mnar <- apply_dropout(ds1, "inc_mnar", b0_mnar, b1_mnar, nu_mnar)

## fit FPCA to each retained dataset
fitFPCA <- function(dd) {
  o <- dd[dd$observed, ]; ids <- sort(unique(o$id))
  Ly <- lapply(ids, function(i) o$y[o$id==i])
  Lt <- lapply(ids, function(i) o$time[o$id==i])
  keep <- sapply(Ly, length) >= 2
  list(fit = FPCA(Ly[keep], Lt[keep],
                  optns=list(dataType="Sparse", FVEthreshold=FVE, verbose=FALSE)),
       ids = ids[keep])
}
FM <- fitFPCA(d_mar); FN <- fitFPCA(d_mnar)

## Selection: find an individual where MAR RETAINS the rise but MNAR does NOT.
## increasing MAR depends on the PREVIOUS value, so the individual leaves only
## after a high value has been recorded -> the rise enters the data.
## increasing MNAR depends on the CURRENT value, so the drop happens at the
## high-value visit -> the rise is never recorded. The right individual is one
## whose largest RETAINED value is much higher under MAR than under MNAR, AND
## which keeps at least two points under BOTH mechanisms so both FPCA fits
## include it.
maxobs <- function(dd) sapply(split(seq_len(nrow(dd)), dd$id), function(ix){
  yy <- dd$y[ix][dd$observed[ix]]; if(length(yy)) max(yy) else -Inf })
ncnt <- function(dd) sapply(split(dd$observed, dd$id), sum)

ids_all  <- sort(unique(ds1$id))
max_mar  <- maxobs(d_mar)[as.character(ids_all)]
max_mnar <- maxobs(d_mnar)[as.character(ids_all)]
n_mar    <- ncnt(d_mar)[as.character(ids_all)]
n_mnar   <- ncnt(d_mnar)[as.character(ids_all)]
n_full   <- as.numeric(table(ds1$id))

gap <- max_mar - max_mnar
## both fits must contain the individual: keep >= 2 points under each mechanism
in_both <- ids_all %in% FM$ids & ids_all %in% FN$ids
ok  <- gap > 6 & n_mar >= 3 & n_mnar >= 2 & n_full >= 6 & in_both
cand <- ids_all[ok]
cand <- cand[order(gap[ok], decreasing = TRUE)]
cat("candidates (MAR records a higher value than MNAR, both fits keep them):",
    length(cand), "\n")
if (!length(cand))
  stop("no clear MAR-records-rise / MNAR-misses-rise individual with both fits; ",
       "try SEED0+2, +3, ... or drop the figure.")

i <- cand[1]
cat("chosen individual:", i,
    " (max retained: MAR=", round(max_mar[as.character(i)],1),
    " MNAR=", round(max_mnar[as.character(i)],1),
    "; points MAR=", n_mar[as.character(i)], " MNAR=", n_mnar[as.character(i)], ")\n")

## data for the chosen individual
b  <- ro$b[[i]]; xt <- traj1(GRID1, b)
di_full <- ds1[ds1$id==i, ]
di_mar  <- d_mar[d_mar$id==i, ]
di_mnar <- d_mnar[d_mnar$id==i, ]
km <- which(FM$ids==i); kn <- which(FN$ids==i)
rm_ <- approx(FM$fit$workGrid, fitted(FM$fit)[km,], xout=GRID1, rule=2)$y
rn_ <- approx(FN$fit$workGrid, fitted(FN$fit)[kn,], xout=GRID1, rule=2)$y

## plot
png("fig_single_recon.png", width=2400, height=1000, res=220)
par(mfrow=c(1,2), mar=c(4.5,4.5,3,1))
yl <- range(xt, rm_, rn_, di_full$y)

## panel 1: increasing MAR
plot(GRID1, xt, type="l", col="black", lwd=1.5, ylim=yl,
     xlab="t", ylab="Simulated marker", main="Increasing MAR, 60%")
points(di_mar$time[di_mar$observed],  di_mar$y[di_mar$observed],  pch=16)
points(di_mar$time[!di_mar$observed], di_mar$y[!di_mar$observed], pch=1)
lines(GRID1, rm_, col="steelblue", lwd=2.5)
grid(lty=3, col="grey80")
legend("topleft", bty="n", cex=0.8,
       legend=c("true trajectory","retained obs","removed obs","FPCA reconstruction"),
       lty=c(1,NA,NA,1), pch=c(NA,16,1,NA), lwd=c(1.5,NA,NA,2.5),
       col=c("black","black","black","steelblue"))

## panel 2: increasing MNAR
plot(GRID1, xt, type="l", col="black", lwd=1.5, ylim=yl,
     xlab="t", ylab="Simulated marker", main="Increasing MNAR, 60%")
points(di_mnar$time[di_mnar$observed],  di_mnar$y[di_mnar$observed],  pch=16)
points(di_mnar$time[!di_mnar$observed], di_mnar$y[!di_mnar$observed], pch=1)
lines(GRID1, rn_, col="firebrick", lwd=2.5)
grid(lty=3, col="grey80")

dev.off()
cat("saved fig_single_recon.png  --  LOOK at it before keeping the figure.\n")
###############################################################################
