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

## High-background variant demo

To show why high-background variants can change the comparison between
mMIRAGE/MIRAGE and an mTADA-style aggregation, run:

```bash
Rscript inst/scripts/background_vs_nonbackground_demo.R
```

The script uses the fixed example data in:

```text
inst/extdata/background_demo_variants.csv
inst/extdata/background_demo_truth.csv
```

It compares three gene classes:

1. `causal_risk`: true disease-risk genes with case enrichment.
2. `high_background_noncausal`: non-risk genes whose variants have elevated
   background counts in both cases and controls.
3. `ordinary_null`: non-risk genes without elevated background counts.

The mTADA-like score is an illustrative aggregation, not a full reimplementation
of the mTADA package. It intentionally uses case counts and the reported
mutation rate only, so high-background noncausal genes can look artificially
strong after aggregation. The MIRAGE/mMIRAGE fit uses both case and control
counts, so genes with similarly elevated counts in cases and controls are
down-weighted. Output files are written to `background_demo_output/`.
