###############################################################################
## 07_adni_missingness_test.R
##
## Two approaches, following the missing-data literature:
##   1. Little's MCAR test (Little 1988) via naniar::mcar_test
##   2. Logistic regression of the missingness indicator on observed variables
##      (previous ADAS11, time, visit number) -- if missingness depends on the
##      PREVIOUS observed value, that is evidence for MAR rather than MCAR; the
##      current (would-be) value is unobservable, so MNAR cannot be tested.
##
## Uses the RAW data MSc_Data.csv (with the missing ADAS11 rows kept).
###############################################################################

## ---- packages ----
## install.packages("naniar")   # if not installed
have_naniar <- requireNamespace("naniar", quietly = TRUE)

raw <- read.csv("MSc_Data.csv", stringsAsFactors = FALSE)
raw$date   <- as.Date(raw$EXAMDATE, format = "%d/%m/%Y")
raw$ADAS11 <- suppressWarnings(as.numeric(raw$ADAS11))
raw <- raw[order(raw$RID, raw$date), ]

cat("Total records:", nrow(raw),
    "| missing ADAS11:", sum(is.na(raw$ADAS11)),
    sprintf("(%.1f%%)\n", 100*mean(is.na(raw$ADAS11))))

## ===========================================================================
## 1. Little's MCAR test
##    naniar::mcar_test expects a data frame of the variables whose joint
##    missingness is being tested. Here ADAS11 is the variable with missingness;
##    to give the test something to compare against, we add per-visit covariates
##    that are always observed (time since first visit, visit number).
## ===========================================================================
raw$years <- ave(as.numeric(raw$date), raw$RID,
                 FUN = function(x) (x - min(x)) / 365.25)
raw$visit <- ave(raw$RID, raw$RID, FUN = seq_along)

if (have_naniar) {
  cat("\n=== Little's MCAR test (naniar::mcar_test) ===\n")
  mc <- naniar::mcar_test(raw[, c("ADAS11", "years", "visit")])
  print(mc)
  cat("\nInterpretation: a small p-value rejects MCAR (missingness is NOT\n")
  cat("completely at random); it does NOT distinguish MAR from MNAR.\n")
} else {
  cat("\n[naniar not installed: run install.packages('naniar') to get\n")
  cat(" Little's MCAR test. Proceeding with the logistic diagnostic below.]\n")
}

## ===========================================================================
## 2. Logistic regression of the missingness indicator
##    Build, for each visit that has a PREVIOUS visit, an indicator of whether
##    THIS visit's ADAS11 is missing, and regress it on the previous observed
##    ADAS11, the time, and the visit number. Dependence on the previous
##    (observed) value is the signature of MAR.
## ===========================================================================
cat("\n=== Logistic regression: P(missing at visit j) ~ observed covariates ===\n")

by <- split(raw, raw$RID)
rows <- list()
for (d in by) {
  d <- d[order(d$date), ]
  if (nrow(d) < 2) next
  for (j in 2:nrow(d)) {
    prev <- d$ADAS11[j-1]
    if (is.na(prev)) next                 # need previous observed value
    rows[[length(rows)+1]] <- data.frame(
      missing_now = as.integer(is.na(d$ADAS11[j])),
      prev_adas   = prev,
      years       = d$years[j],
      visit       = j)
  }
}
reg <- do.call(rbind, rows)
cat(sprintf("Rows available for the model: %d (missing rate %.1f%%)\n",
            nrow(reg), 100*mean(reg$missing_now)))

fit <- glm(missing_now ~ prev_adas + years + visit,
           data = reg, family = binomial)
cat("\nCoefficients (log-odds of being missing):\n")
print(summary(fit)$coefficients)

## odds ratios for interpretability
cat("\nOdds ratios (exp(coef)):\n")
print(round(exp(coef(fit)), 4))

cat("\nInterpretation:\n")
cat("- prev_adas: if its coefficient is significant, missingness depends on the\n")
cat("  PREVIOUS observed score -> evidence of MAR (dependence on observed data).\n")
cat("  A negative coefficient means higher previous score -> LESS likely missing\n")
cat("  (i.e. better-scoring, lower-ADAS patients miss slightly more), matching\n")
cat("  the earlier finding that pre-missing values are a little lower.\n")
cat("- years / visit: dependence on time or visit number is MCAR-like with\n")
cat("  respect to the outcome.\n")
cat("- MNAR (dependence on the CURRENT, unobserved value) cannot be tested here,\n")
cat("  because that value is exactly what is missing.\n")

## quick summary numbers used in the text
cat(sprintf("\nMean previous ADAS11 when THIS visit is missing:    %.2f\n",
            mean(reg$prev_adas[reg$missing_now==1])))
cat(sprintf("Mean previous ADAS11 when THIS visit is observed:   %.2f\n",
            mean(reg$prev_adas[reg$missing_now==0])))
###############################################################################
