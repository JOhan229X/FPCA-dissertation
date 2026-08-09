# Robustness of FPCA to Informative Missingness

R code for my MSc dissertation (MATH5872M, University of Leeds):
a simulation study and real-data application on how robust functional
principal component analysis (FPCA) is to missing data.

## Requirements
- R (>= 4.0)
- packages: `fdapace`, `parallel` (base), `naniar` (for the MCAR test)

```r
install.packages(c("fdapace", "naniar"))
```

## Scripts (run in order)

| Script | What it does |
|--------|--------------|
| `01_simulation_main.R`       | Main Monte Carlo: two datasets, six dropout mechanisms, 100 replications; saves reconstruction / mean / eigenfunction MISE |
| `02_baseline.R`              | No-missingness baseline (source `01` first) |
| `03_plot_simulation.R`       | Figures: reconstruction MISE (with baseline) and component MISE |
| `04_adni_explore.R`          | ADNI cleaning, missingness pattern, spaghetti plot |
| `05_adni_fpca_log.R`         | ADNI FPCA on log(ADAS11+1) with figures (primary analysis) |
| `06_adni_fpca_raw.R`         | ADNI FPCA on the raw scale (sensitivity, appendix) |
| `07_adni_missingness_test.R` | Little's MCAR test + logistic regression of the missingness indicator |

`02` and `03` reuse objects from `01`, so run `source("01_simulation_main.R")` first.

## Data
The ADNI data (`MSc_Data.csv`) are **not included** here, as they are subject to
a data use agreement. The scripts expect that file in the working directory.

## Reproducibility
All random seeds are fixed inside the scripts (replication r uses seed
`SEED0 + r`; calibration uses seed 2024), so the results are reproducible.
