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
inst/extdata/background_demo_meta.csv
```

These files are not newly generated toy data. They are the first replicate from
the original control-informative simulation used in the manuscript exploration:
500 genes, 15 variants per gene per trait, 3,000 cases and 3,000 controls per
trait, 5% high-background variants, and a 5-fold background-rate multiplier.

The script compares three gene classes:

1. `trait1_risk`: true Trait 1 risk genes with case enrichment.
2. `nonrisk_with_high_background`: non-risk genes containing high-background
   variants, whose counts are elevated in both cases and controls.
3. `ordinary_nonrisk`: non-risk genes without high-background variants.

The mTADA-like score is an illustrative aggregation, not a full reimplementation
of the mTADA package. It intentionally uses case counts and the reported
mutation rate only, so high-background noncausal genes can look artificially
strong after aggregation. The MIRAGE/mMIRAGE fit uses both case and control
counts, so genes with similarly elevated counts in cases and controls are
down-weighted. Output files are written to `background_demo_output/`, including
`background_demo_auc.csv`.
