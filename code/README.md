# Code

Save command-line scripts and shared R code here.

## Core two-trait MIRAGE code

- `two_trait_mirage.R`: fits the two-trait gene-level four-state MIRAGE model.
  The four gene states are null in both traits (`00`), trait 1 only (`10`),
  trait 2 only (`01`), and shared risk (`11`).
- `simulate_two_trait_mirage.R`: creates small simulated data sets for examples
  and smoke tests.
- `smoke_test_two_trait_mirage.R`: runs a minimal end-to-end check from
  simulation to model fit.

To run the smoke test from the repository root:

```sh
Rscript code/smoke_test_two_trait_mirage.R
```
