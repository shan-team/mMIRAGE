# mMIRAGE

`mmirage` implements a two-trait extension of MIRAGE for joint rare-variant
gene discovery from variant-level case-control summary counts.

## Install locally

From the repository root:

```r
install.packages(".", repos = NULL, type = "source")
```

or:

```bash
R CMD INSTALL .
```

## Minimal example

```r
library(mmirage)

sim <- simulate_case_control_two_trait(seed = 1)
fit <- mirage_two_trait(
  sim$trait1, sim$trait2,
  n1_trait1 = 3000, n0_trait1 = 3000,
  n1_trait2 = 3000, n0_trait2 = 3000,
  log_bf_col_trait1 = "log_bf",
  log_bf_col_trait2 = "log_bf",
  pi_init = sim$parameters$pi,
  eta_init_trait1 = sim$parameters$eta_trait1,
  eta_init_trait2 = sim$parameters$eta_trait2
)

head(fit$posterior)
fit$parameters$pi
```

## Reproducibility script

The manuscript-style one-replicate reproducibility script is:

```bash
Rscript inst/scripts/reproduce_one_case.R
```

It generates one 5% high-background simulation, fits joint mMIRAGE and a
single-trait MIRAGE comparator, runs burden tests, and writes summary files
under `reproducibility_output/`.

## High-background variant demo

A smaller standalone demo is included to explain why high-background variants
can favor mMIRAGE/MIRAGE over an mTADA-style case-count aggregation:

```bash
Rscript inst/scripts/background_vs_nonbackground_demo.R
```

The script reads the included simulation files under `inst/extdata/`, compares
causal, high-background noncausal, and ordinary null genes, and writes a compact
table plus plot under `background_demo_output/`. The mTADA-like score is an
illustrative case-count aggregation for this diagnostic demo, not a full
reimplementation of the mTADA package.
