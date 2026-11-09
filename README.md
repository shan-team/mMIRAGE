# mMIRAGE

**An R package for joint rare-variant gene discovery in two traits.**

mMIRAGE extends single-trait MIRAGE to analyze variant-level case-control
allele counts from two studies. It estimates the probability that each gene
is associated with trait 1, trait 2, or both traits. The method combines
functional annotations with case-control evidence and allows information
sharing between traits without requiring the same variants to be observed
in both studies.

## Basic idea

The model has two levels:

1. **Gene level.** Each gene belongs to one of four risk states:
   `00` (neither trait), `10` (trait 1 only), `01` (trait 2 only), or
   `11` (both traits). Their proportions are denoted by `pi`.
2. **Variant level.** Within a risk gene, a variant contributes to disease risk
   with an annotation-specific probability `eta`. Nonrisk variants can occur
   in both risk and nonrisk genes. Variant relative risks follow
   annotation-specific Gamma distributions.

For each variant, the likelihood models its case allele count conditional on
its total case-control allele count. Relative risks are integrated over the
Gamma distribution to obtain variant Bayes factors, which are combined within
genes and across traits. The four-state model lets evidence from one trait
inform inference for the other when the fitted gene-risk architecture supports
sharing. `pi` and `eta` are estimated by expectation-maximization (EM).

## Installation

R version 4.1.0 or later is required. Install from GitHub:

```r
# Run once if remotes is not installed:
# install.packages("remotes")
remotes::install_github("shan-team/mMIRAGE")
library(mmirage)
```

Alternatively, clone or download this repository and run `R CMD INSTALL .`
from its root directory.

## Input data

Provide **two data frames**, one for each trait. Each row represents a variant
within a gene, with the following columns:

| Column | Description |
|---|---|
| `Gene` | Gene identifier, using the same naming convention in both studies. |
| `ID` | Variant identifier. Optional; identifiers are generated if absent. |
| `No.case` | Alternate allele count in cases, as a nonnegative integer. |
| `No.contr` | Alternate allele count in controls, as a nonnegative integer. |
| `category` | Functional annotation category, such as `LoF` or `Missense`. |

These are **allele counts, not numbers of variant carriers**. Column names can
be changed through `gene_col`, `id_col`, `case_col`, `control_col`, and
`category_col`.

The studies may contain different variants or annotation categories. By default,
the analysis uses the union of genes observed in either study. A gene with no
usable variants in one study receives neutral evidence (gene BF = 1) for that
study; missing variants are not filled with artificial zero counts.

Sample sizes must be supplied separately. The software does not infer them
from the allele-count tables. Study data should be quality controlled and gene
identifiers harmonized before fitting. The model assumes conditional independence
of the two studies given gene-risk state.

## Parameters supplied by the user

| Argument | Meaning | Default |
|---|---|---|
| `n1_trait1`, `n0_trait1` | Numbers of cases and controls for trait 1. | Required |
| `n1_trait2`, `n0_trait2` | Numbers of cases and controls for trait 2. | Required |
| `gamma_trait1`, `gamma_trait2` | Mean relative risk of a risk variant in each trait, optionally by annotation. | `3` |
| `sigma_trait1`, `sigma_trait2` | Rate of the Gamma relative-risk distribution, optionally by annotation. | `2` |
| `pi_init` | Starting proportions in state order `00`, `10`, `01`, `11`. | `(0.90, 0.04, 0.04, 0.02)` |
| `eta_init_trait1`, `eta_init_trait2` | Starting risk-variant probabilities, optionally by annotation. | `0.1` |
| `max_iter` | Maximum number of EM iterations. | `500` |
| `tolerance_parameter` | Maximum absolute parameter-change threshold. | `1e-6` |
| `tolerance_loglik` | Relative log-likelihood-change threshold. | `1e-8` |

The Gamma distribution uses **shape = `gamma * sigma` and rate = `sigma`**, so
its mean is `gamma` and variance is `gamma / sigma`. These hyperparameters are
specified by the user, not estimated by this EM routine. Choose them based on
scientific knowledge or a justified external analysis.

A scalar applies to every annotation category. A named vector specifies
category-specific values, for example:

```r
gamma_trait1 <- c(LoF = 6, Missense = 3)
sigma_trait1 <- c(LoF = 2, Missense = 2)
```

`pi_init` and the `eta_init_*` arguments are starting values: the corresponding
parameters are estimated by default, rather than fixed at these values.
Simulation truth, baseline allele frequencies, and mutation rates are **not
required inputs**. EM runs from one starting point per function call; to assess
multiple starts, fit them explicitly and compare final log likelihoods.

## Minimal working example

The following small, illustrative count tables demonstrate the interface; they
are not intended to establish statistical power or calibration.

```r
library(mmirage)

trait1 <- data.frame(
  Gene = c("G1", "G1", "G2", "G2", "G3", "G3"),
  ID = paste0("t1_v", 1:6),
  No.case = c(12, 5, 1, 2, 3, 1),
  No.contr = c(1, 1, 2, 2, 3, 2),
  category = rep(c("LoF", "Missense"), 3)
)
trait2 <- data.frame(
  Gene = c("G1", "G1", "G2", "G2", "G3", "G3"),
  ID = paste0("t2_v", 1:6),
  No.case = c(8, 4, 2, 1, 2, 1),
  No.contr = c(1, 1, 2, 2, 2, 2),
  category = rep(c("LoF", "Missense"), 3)
)

fit <- mirage_two_trait(
  data_trait1 = trait1,
  data_trait2 = trait2,
  n1_trait1 = 3000, n0_trait1 = 3000,
  n1_trait2 = 3000, n0_trait2 = 3000,
  trait_names = c("Trait1", "Trait2"),
  gamma_trait1 = c(LoF = 6, Missense = 3),
  sigma_trait1 = c(LoF = 2, Missense = 2),
  gamma_trait2 = c(LoF = 6, Missense = 3),
  sigma_trait2 = c(LoF = 2, Missense = 2),
  max_iter = 10000
)

print(fit)
fit$parameters
fit$posterior

# Select trait-1 genes at nominal 5% Bayesian FDR:
subset(fit$bfdr$trait1, bayesian_fdr <= 0.05)

# Use the shared-state posterior for shared-risk discoveries:
subset(fit$bfdr$pleiotropy, bayesian_fdr <= 0.05)

# Check numerical convergence:
fit$diagnostics$converged
fit$diagnostics$n_iter
fit$diagnostics$loglik_non_decreasing
```

## Outputs

`mirage_two_trait()` returns an object of class `mirage_two_trait`.

| Component | Contents |
|---|---|
| `fit$parameters` | Estimated four-state proportions `pi`, marginal risk-gene proportions `delta_trait1` and `delta_trait2`, and annotation-specific `eta_trait1` and `eta_trait2`. |
| `fit$posterior` | One row per gene: four-state probabilities `tau_00`, `tau_10`, `tau_01`, `tau_11`, marginal probabilities `PP_trait1` and `PP_trait2`, and shared-risk probability `PP_pleiotropy`. |
| `fit$gene_bf` | Trait-specific gene BFs and log BFs, variant counts, and study-coverage indicators. |
| `fit$bfdr` | Separate ranked tables for `trait1`, `trait2`, and `pleiotropy`, containing posterior probability, local FDR, cumulative Bayesian FDR, and rank. |
| `fit$diagnostics` | EM convergence flag, iteration count, likelihood monotonicity, boundary-pi flags, annotation-weight diagnostics, and missing-study information. |
| `fit$traces` | Log likelihood, maximum parameter change, pi and eta trajectories over EM iterations. |
| `fit$variant_info` | Standardized variant tables with log BFs and conditional variant-risk probabilities. |

A Bayesian FDR table ranks genes by posterior probability and reports the
cumulative mean of `1 - posterior_probability`. Select the leading genes whose
cumulative Bayesian FDR is below the desired threshold. A shared-risk discovery
uses `PP_pleiotropy`; it is not the same criterion as intersecting the two
marginal discovery lists.

## Convergence and interpretation

EM is marked converged only when **both** parameter change and relative
log-likelihood change pass the requested thresholds. Reaching `max_iter` alone
does not establish convergence. Inspect `fit$diagnostics` and the parameter
trajectories; increase the iteration limit or compare alternative starts when
needed.

Numerical convergence does not establish global optimality, parameter
identifiability, or posterior/FDR calibration. Boundary estimates and sparse
case-control information require particular care. Posterior probabilities are
computed using fitted parameters; this implementation does not sample a full
Bayesian posterior over pi and eta.

## Additional options

- **Precomputed BFs:** supply the names of finite natural-log variant BF columns
  through `log_bf_col_trait1` and `log_bf_col_trait2`. They are used unchanged;
  verify their calculation and provenance. Without these columns, the function
  computes Gamma-integrated BFs from counts.
- **Fixed parameters:** `fixed_pi` fixes the complete four-state vector. The
  current `estimate_eta` flag applies jointly to both traits: supplying either
  `fixed_eta_trait1` or `fixed_eta_trait2` disables eta estimation for both.
  To fix eta deliberately, provide both named vectors.
- **Gene alignment:** `gene_universe` can be `"union"`, `"intersection"`,
  `"trait1"`, or `"trait2"`.

Use `formals(mirage_two_trait)` to inspect the complete function interface.
The package also provides `gene_burden_test()` for a gene-level Fisher burden
baseline and `bfdr_table()` for constructing a ranked Bayesian FDR table.

## License

MIT; see [LICENSE](LICENSE).
