# Reproducibility Scripts

This folder contains lightweight scripts for reproducing the manuscript-style
workflow without running the full 100-replicate simulation.

## One-case simulation

From the repository root, first install the package:

```bash
R CMD INSTALL .
```

Then run:

```bash
Rscript inst/scripts/reproduce_one_case.R
```

The script:

1. Simulates one two-trait case-control rare-variant dataset.
2. Uses the manuscript's 5% high-background variant setting with a 5-fold
   background-rate multiplier.
3. Fits joint mMIRAGE.
4. Fits a single-trait MIRAGE comparator through an independent fixed prior.
5. Runs gene-level burden tests for both traits.
6. Writes summary files to `reproducibility_output/`.

This script is intended as a transparent reproducibility example. It does not
repeat the full 100-replicate simulation grid.
