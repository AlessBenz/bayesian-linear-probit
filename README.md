# Bayesian Ordinal Probit — Gibbs + Metropolis-Hastings (R)

This repo is a cleaned, GitHub-friendly structure extracted from the original single Rmd:
`analysis/02_original_report.Rmd`.

## What’s inside

- `analysis/01_report.Rmd`: a reproducible report template
- `analysis/00_run_all.R`: runs the pipeline end-to-end (no knitting)
- `R/ordinal_probit_sampler.R`: Gibbs sampler for β and z + MH updates for cutpoints γ
- `R/posterior_predictive.R`: posterior predictive simulation helpers
- `R/data_prep.R`: data loading + preprocessing

## How to run

1. Put the dataset CSV in `data/raw/vino_qualita.csv` (see `data/README.md`).
2. In R (from the repo root):

```r
install.packages(c("here","truncnorm","MASS","bayesplot","coda","ggplot2","rmarkdown"))
source("analysis/00_run_all.R")
rmarkdown::render("analysis/01_report.Rmd")
```

## Notes / small fixes compared to the original Rmd

- Uses relative paths via `here::here()` instead of a local Windows path.
- (Optional) Adds an intercept column and standardizes only non-intercept predictors.
- Posterior predictive simulation can add probit noise `N(0,1)` before discretizing (recommended).
