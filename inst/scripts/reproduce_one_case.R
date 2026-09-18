#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mmirage)
})

outdir <- "reproducibility_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

seed <- 20260813
n_case <- 3000
n_control <- 3000

sim <- simulate_case_control_two_trait(
  n_genes = 500,
  variants_per_gene = 15,
  pi = c("00" = 0.94, "10" = 0.04, "01" = 0.01, "11" = 0.01),
  eta_trait1 = c("C1" = 0.05, "C2" = 0.20, "C3" = 0.50),
  eta_trait2 = c("C1" = 0.03, "C2" = 0.15, "C3" = 0.40),
  n_case_trait1 = n_case,
  n_control_trait1 = n_control,
  n_case_trait2 = n_case,
  n_control_trait2 = n_control,
  background_prob = 0.05,
  background_multiplier = 5,
  seed = seed
)

joint_fit <- mirage_two_trait(
  sim$trait1, sim$trait2,
  n1_trait1 = n_case, n0_trait1 = n_control,
  n1_trait2 = n_case, n0_trait2 = n_control,
  log_bf_col_trait1 = "log_bf",
  log_bf_col_trait2 = "log_bf",
  pi_init = sim$parameters$pi,
  eta_init_trait1 = sim$parameters$eta_trait1,
  eta_init_trait2 = sim$parameters$eta_trait2,
  max_iter = 100
)

delta1 <- unname(sim$parameters$pi[["10"]] + sim$parameters$pi[["11"]])
delta2 <- unname(sim$parameters$pi[["01"]] + sim$parameters$pi[["11"]])
independent_pi <- c(
  "00" = (1 - delta1) * (1 - delta2),
  "10" = delta1 * (1 - delta2),
  "01" = (1 - delta1) * delta2,
  "11" = delta1 * delta2
)

single_fit <- mirage_two_trait(
  sim$trait1, sim$trait2,
  n1_trait1 = n_case, n0_trait1 = n_control,
  n1_trait2 = n_case, n0_trait2 = n_control,
  log_bf_col_trait1 = "log_bf",
  log_bf_col_trait2 = "log_bf",
  fixed_pi = independent_pi,
  eta_init_trait1 = sim$parameters$eta_trait1,
  eta_init_trait2 = sim$parameters$eta_trait2,
  max_iter = 100
)

burden1 <- gene_burden_test(sim$trait1, n_case, n_control)
burden2 <- gene_burden_test(sim$trait2, n_case, n_control)

truth <- sim$truth[match(joint_fit$posterior$Gene, sim$truth$Gene), ]
burden1 <- burden1[match(joint_fit$posterior$Gene, burden1$Gene), ]
burden2 <- burden2[match(joint_fit$posterior$Gene, burden2$Gene), ]

auc_table <- data.frame(
  method = c(
    "burden_test", "single_trait_MIRAGE", "joint_mMIRAGE",
    "burden_test", "single_trait_MIRAGE", "joint_mMIRAGE",
    "single_trait_MIRAGE", "joint_mMIRAGE"
  ),
  target = c(
    "trait1", "trait1", "trait1",
    "trait2", "trait2", "trait2",
    "shared", "shared"
  ),
  auc = c(
    binary_auc(burden1$score, truth$trait1_risk),
    binary_auc(single_fit$posterior$PP_trait1, truth$trait1_risk),
    binary_auc(joint_fit$posterior$PP_trait1, truth$trait1_risk),
    binary_auc(burden2$score, truth$trait2_risk),
    binary_auc(single_fit$posterior$PP_trait2, truth$trait2_risk),
    binary_auc(joint_fit$posterior$PP_trait2, truth$trait2_risk),
    binary_auc(single_fit$posterior$PP_pleiotropy, truth$pleiotropic),
    binary_auc(joint_fit$posterior$PP_pleiotropy, truth$pleiotropic)
  )
)

parameter_table <- data.frame(
  parameter = names(sim$parameters$pi),
  truth = as.numeric(sim$parameters$pi),
  joint_mMIRAGE = as.numeric(joint_fit$parameters$pi[names(sim$parameters$pi)]),
  single_trait_fixed_prior = as.numeric(independent_pi[names(sim$parameters$pi)])
)

fdr_counts <- data.frame(
  method = c("burden_test", "single_trait_MIRAGE", "joint_mMIRAGE"),
  trait1_fdr05 = c(
    sum(burden1$bh_fdr <= 0.05, na.rm = TRUE),
    sum(single_fit$bfdr$trait1$bayesian_fdr <= 0.05, na.rm = TRUE),
    sum(joint_fit$bfdr$trait1$bayesian_fdr <= 0.05, na.rm = TRUE)
  ),
  trait2_fdr05 = c(
    sum(burden2$bh_fdr <= 0.05, na.rm = TRUE),
    sum(single_fit$bfdr$trait2$bayesian_fdr <= 0.05, na.rm = TRUE),
    sum(joint_fit$bfdr$trait2$bayesian_fdr <= 0.05, na.rm = TRUE)
  ),
  shared_fdr05 = c(
    NA_integer_,
    sum(single_fit$bfdr$pleiotropy$bayesian_fdr <= 0.05, na.rm = TRUE),
    sum(joint_fit$bfdr$pleiotropy$bayesian_fdr <= 0.05, na.rm = TRUE)
  )
)

write.csv(auc_table, file.path(outdir, "one_case_auc.csv"), row.names = FALSE)
write.csv(parameter_table, file.path(outdir, "one_case_parameter_estimates.csv"), row.names = FALSE)
write.csv(fdr_counts, file.path(outdir, "one_case_fdr05_counts.csv"), row.names = FALSE)
write.csv(joint_fit$posterior, file.path(outdir, "one_case_joint_posterior.csv"), row.names = FALSE)

cat("\nOne-case reproducibility run complete.\n")
cat("Output directory:", normalizePath(outdir), "\n\n")
print(auc_table)
cat("\nParameter estimates:\n")
print(parameter_table)
cat("\nFDR < 5% discovery counts:\n")
print(fdr_counts)
