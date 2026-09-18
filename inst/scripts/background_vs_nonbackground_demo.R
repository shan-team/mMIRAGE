#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mmirage)
})

`%||%` <- function(x, y) if (length(x) == 0 || is.na(x)) y else x

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", file_arg[1] %||% NA_character_)
script_dir <- if (is.na(script_path)) getwd() else dirname(normalizePath(script_path))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
extdata_dir <- file.path(repo_root, "inst", "extdata")

variant_file <- file.path(extdata_dir, "background_demo_variants.csv")
truth_file <- file.path(extdata_dir, "background_demo_truth.csv")
meta_file <- file.path(extdata_dir, "background_demo_meta.csv")
outdir <- "background_demo_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(variant_file) || !file.exists(truth_file)) {
  stop(
    "Missing demo data. Expected files under inst/extdata/: ",
    "background_demo_variants.csv and background_demo_truth.csv.",
    call. = FALSE
  )
}

variants <- read.csv(variant_file, stringsAsFactors = FALSE)
truth <- read.csv(truth_file, stringsAsFactors = FALSE)
meta <- if (file.exists(meta_file)) read.csv(meta_file, stringsAsFactors = FALSE) else NULL

required <- c(
  "trait", "ID", "Gene", "No.case", "No.contr", "category", "log_bf",
  "mutation_rate", "count_rate", "high_background", "true_z"
)
missing <- setdiff(required, names(variants))
if (length(missing) > 0) {
  stop("background_demo_variants.csv is missing columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

trait1 <- variants[variants$trait == "trait1", ]
trait2 <- variants[variants$trait == "trait2", ]

n_case1 <- if (!is.null(meta)) meta$n_case1[1] else 3000
n_control1 <- if (!is.null(meta)) meta$n_control1[1] else 3000
n_case2 <- if (!is.null(meta)) meta$n_case2[1] else 3000
n_control2 <- if (!is.null(meta)) meta$n_control2[1] else 3000

pi_init <- c("00" = 0.94, "10" = 0.04, "01" = 0.01, "11" = 0.01)
eta1 <- c("C1" = 0.05, "C2" = 0.20, "C3" = 0.50)
eta2 <- c("C1" = 0.03, "C2" = 0.15, "C3" = 0.40)
rr_mean <- c("C1" = 3, "C2" = 3, "C3" = 5)
rr_sigma <- c("C1" = 1, "C2" = 1, "C3" = 1)

fit <- mirage_two_trait(
  trait1,
  trait2,
  n1_trait1 = n_case1,
  n0_trait1 = n_control1,
  n1_trait2 = n_case2,
  n0_trait2 = n_control2,
  gamma_trait1 = rr_mean,
  sigma_trait1 = rr_sigma,
  gamma_trait2 = rr_mean,
  sigma_trait2 = rr_sigma,
  log_bf_col_trait1 = "log_bf",
  log_bf_col_trait2 = "log_bf",
  pi_init = pi_init,
  eta_init_trait1 = eta1,
  eta_init_trait2 = eta2,
  max_iter = 100
)

mtada_like_gene_score <- function(d, n_case) {
  agg <- aggregate(cbind(No.case, mutation_rate) ~ Gene + category, d, sum)
  gene_rows <- lapply(split(agg, agg$Gene), function(z) {
    lambda <- 2 * n_case * z$mutation_rate
    pvals <- ppois(
      z$No.case - 1,
      lambda = pmax(lambda, .Machine$double.xmin),
      lower.tail = FALSE
    )
    data.frame(
      Gene = z$Gene[1],
      mtada_like_score = sum(-log10(pmax(pvals, .Machine$double.xmin))),
      mtada_like_min_p = min(pvals),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, gene_rows)
}

case_control_gene_p <- function(d, n_case, n_control) {
  agg <- aggregate(cbind(No.case, No.contr) ~ Gene, d, sum)
  pvals <- vapply(seq_len(nrow(agg)), function(i) {
    test <- prop.test(
      x = c(agg$No.case[i], agg$No.contr[i]),
      n = c(2 * n_case, 2 * n_control),
      alternative = "greater",
      correct = FALSE
    )
    test$p.value
  }, numeric(1))
  data.frame(
    Gene = agg$Gene,
    burden_case_control_p = pvals,
    stringsAsFactors = FALSE
  )
}

count_summary <- aggregate(
  cbind(No.case, No.contr, high_background, true_z) ~ Gene,
  trait1,
  sum
)
rate_summary <- aggregate(cbind(mutation_rate, count_rate) ~ Gene, trait1, sum)
gene_summary <- Reduce(function(x, y) merge(x, y, by = "Gene"), list(
  truth,
  count_summary,
  rate_summary,
  mtada_like_gene_score(trait1, n_case1),
  case_control_gene_p(trait1, n_case1, n_control1),
  fit$posterior[, c("Gene", "PP_trait1", "PP_trait2", "PP_pleiotropy")]
))

gene_summary$background_gene <- gene_summary$high_background > 0
gene_summary$gene_group <- ifelse(
  gene_summary$trait1_risk,
  "trait1_risk",
  ifelse(gene_summary$background_gene, "nonrisk_with_high_background", "ordinary_nonrisk")
)
gene_summary$case_control_ratio <- (gene_summary$No.case + 0.5) /
  (gene_summary$No.contr + 0.5)
gene_summary <- gene_summary[order(-gene_summary$mtada_like_score), ]

group_summary <- aggregate(
  cbind(
    No.case, No.contr, high_background, mutation_rate, count_rate,
    mtada_like_score, burden_case_control_p, PP_trait1, case_control_ratio
  ) ~ gene_group,
  gene_summary,
  mean
)

write.csv(gene_summary, file.path(outdir, "background_demo_gene_scores.csv"), row.names = FALSE)
write.csv(group_summary, file.path(outdir, "background_demo_group_summary.csv"), row.names = FALSE)

png(file.path(outdir, "background_demo_score_comparison.png"),
    width = 1500, height = 1100, res = 180)
cols <- c(
  trait1_risk = "#1f77b4",
  nonrisk_with_high_background = "#d95f02",
  ordinary_nonrisk = "#7f7f7f"
)
plot(
  gene_summary$mtada_like_score,
  gene_summary$PP_trait1,
  pch = 19,
  col = cols[gene_summary$gene_group],
  xlab = "mTADA-like aggregated score (case counts + mutation rate)",
  ylab = "MIRAGE posterior probability for Trait 1",
  main = "Original 5% high-background simulation, replicate 1"
)
legend("topleft", legend = names(cols), col = cols, pch = 19, bty = "n")
dev.off()

cat("\nBackground vs non-background comparison using the original simulation data.\n")
cat("Input data:\n")
cat("  ", normalizePath(variant_file), "\n", sep = "")
cat("  ", normalizePath(truth_file), "\n", sep = "")
if (!is.null(meta)) {
  cat("\nSimulation setting:\n")
  print(meta)
}
cat("\nOutput directory:\n")
cat("  ", normalizePath(outdir), "\n\n", sep = "")

cat("Group-level averages for Trait 1 genes:\n")
print(group_summary)

cat("\nTop genes by mTADA-like aggregation:\n")
print(head(gene_summary[, c(
  "Gene", "gene_group", "state", "No.case", "No.contr", "high_background",
  "true_z", "mtada_like_score", "burden_case_control_p", "PP_trait1"
)], 15))

cat("\nTop genes by MIRAGE posterior:\n")
mirage_rank <- gene_summary[order(-gene_summary$PP_trait1), ]
print(head(mirage_rank[, c(
  "Gene", "gene_group", "state", "No.case", "No.contr", "high_background",
  "true_z", "mtada_like_score", "burden_case_control_p", "PP_trait1"
)], 15))
