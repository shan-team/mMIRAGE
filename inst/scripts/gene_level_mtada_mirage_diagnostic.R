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
outdir <- "gene_diagnostic_output"
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

mtada_like_gene_score <- function(d, n_case, suffix) {
  agg <- aggregate(cbind(No.case, mutation_rate) ~ Gene + category, d, sum)
  rows <- lapply(split(agg, agg$Gene), function(z) {
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
  out <- do.call(rbind, rows)
  names(out)[names(out) == "mtada_like_score"] <- paste0("mtada_like_score_", suffix)
  names(out)[names(out) == "mtada_like_min_p"] <- paste0("mtada_like_min_p_", suffix)
  out
}

variant_gene_summary <- function(d, suffix) {
  count_summary <- aggregate(
    cbind(No.case, No.contr, high_background, true_z) ~ Gene,
    d,
    sum
  )
  rate_summary <- aggregate(cbind(mutation_rate, count_rate) ~ Gene, d, sum)
  out <- merge(count_summary, rate_summary, by = "Gene")
  names(out)[names(out) != "Gene"] <- paste0(names(out)[names(out) != "Gene"], "_", suffix)
  out
}

gene_scores <- Reduce(function(x, y) merge(x, y, by = "Gene"), list(
  truth,
  variant_gene_summary(trait1, "trait1"),
  variant_gene_summary(trait2, "trait2"),
  mtada_like_gene_score(trait1, n_case1, "trait1"),
  mtada_like_gene_score(trait2, n_case2, "trait2"),
  fit$posterior[, c("Gene", "PP_trait1", "PP_trait2", "PP_pleiotropy")],
  fit$gene_bf[, c("Gene", "log_B_trait1", "log_B_trait2")]
))

gene_scores$trait1_case_control_ratio <- (gene_scores$No.case_trait1 + 0.5) /
  (gene_scores$No.contr_trait1 + 0.5)
gene_scores$trait2_case_control_ratio <- (gene_scores$No.case_trait2 + 0.5) /
  (gene_scores$No.contr_trait2 + 0.5)

gene_scores$mtada_rank_trait1 <- rank(-gene_scores$mtada_like_score_trait1, ties.method = "first")
gene_scores$mirage_rank_trait1 <- rank(-gene_scores$PP_trait1, ties.method = "first")
gene_scores$mtada_rank_trait2 <- rank(-gene_scores$mtada_like_score_trait2, ties.method = "first")
gene_scores$mirage_rank_trait2 <- rank(-gene_scores$PP_trait2, ties.method = "first")

risk_gene_trait1 <- gene_scores$Gene[which.max(ifelse(
  gene_scores$trait1_risk,
  gene_scores$mtada_like_score_trait1 - scale(gene_scores$PP_trait1)[, 1],
  -Inf
))]
nonrisk_high_background_trait1 <- gene_scores$Gene[which.max(ifelse(
  !gene_scores$trait1_risk & gene_scores$high_background_trait1 > 0,
  gene_scores$mtada_like_score_trait1 - scale(gene_scores$PP_trait1)[, 1],
  -Inf
))]
risk_gene_trait2 <- gene_scores$Gene[which.max(ifelse(
  gene_scores$trait2_risk,
  gene_scores$mtada_like_score_trait2 - scale(gene_scores$PP_trait2)[, 1],
  -Inf
))]
nonrisk_high_background_trait2 <- gene_scores$Gene[which.max(ifelse(
  !gene_scores$trait2_risk & gene_scores$high_background_trait2 > 0,
  gene_scores$mtada_like_score_trait2 - scale(gene_scores$PP_trait2)[, 1],
  -Inf
))]

selected <- data.frame(
  comparison = c(
    "trait1_risk_gene",
    "trait1_nonrisk_high_background_gene",
    "trait2_risk_gene",
    "trait2_nonrisk_high_background_gene"
  ),
  Gene = c(
    risk_gene_trait1,
    nonrisk_high_background_trait1,
    risk_gene_trait2,
    nonrisk_high_background_trait2
  ),
  stringsAsFactors = FALSE
)

selected_gene_summary <- merge(selected, gene_scores, by = "Gene", all.x = TRUE)
selected_gene_summary <- selected_gene_summary[order(selected_gene_summary$comparison), ]

selected_variant_evidence <- variants[variants$Gene %in% selected$Gene, ]
selected_variant_evidence$case_control_ratio <- (selected_variant_evidence$No.case + 0.5) /
  (selected_variant_evidence$No.contr + 0.5)
selected_variant_evidence <- merge(
  selected[, c("comparison", "Gene")],
  selected_variant_evidence,
  by = "Gene",
  all.x = TRUE
)
selected_variant_evidence <- selected_variant_evidence[
  order(selected_variant_evidence$comparison, selected_variant_evidence$trait,
        -selected_variant_evidence$log_bf),
]

category_evidence <- aggregate(
  cbind(No.case, No.contr, mutation_rate, count_rate, high_background, true_z) ~
    Gene + trait + category,
  variants[variants$Gene %in% selected$Gene, ],
  sum
)
category_evidence$case_control_ratio <- (category_evidence$No.case + 0.5) /
  (category_evidence$No.contr + 0.5)
category_evidence <- merge(selected[, c("comparison", "Gene")], category_evidence, by = "Gene")

write.csv(gene_scores, file.path(outdir, "all_gene_scores.csv"), row.names = FALSE)
write.csv(selected_gene_summary, file.path(outdir, "selected_gene_summary.csv"), row.names = FALSE)
write.csv(category_evidence, file.path(outdir, "selected_category_evidence.csv"), row.names = FALSE)
write.csv(selected_variant_evidence, file.path(outdir, "selected_variant_evidence.csv"), row.names = FALSE)

cat("\nGene-level mTADA-like vs MIRAGE diagnostic complete.\n")
cat("Input data:\n")
cat("  ", normalizePath(variant_file), "\n", sep = "")
cat("  ", normalizePath(truth_file), "\n", sep = "")
cat("\nOutput directory:\n")
cat("  ", normalizePath(outdir), "\n\n", sep = "")

cat("Selected gene-level examples:\n")
print(selected_gene_summary[, c(
  "comparison", "Gene", "state", "trait1_risk", "trait2_risk",
  "No.case_trait1", "No.contr_trait1", "high_background_trait1",
  "mtada_like_score_trait1", "PP_trait1", "mtada_rank_trait1", "mirage_rank_trait1",
  "No.case_trait2", "No.contr_trait2", "high_background_trait2",
  "mtada_like_score_trait2", "PP_trait2", "mtada_rank_trait2", "mirage_rank_trait2"
)])

cat("\nInterpretation guide:\n")
cat("* mTADA-like aggregation uses case counts plus mutation rate.\n")
cat("* MIRAGE/mMIRAGE uses case and control counts at the variant level.\n")
cat("* A non-risk high-background gene can have a high mTADA-like score if its case counts are high relative to the original mutation rate.\n")
cat("* MIRAGE/mMIRAGE downweights that pattern when controls are also elevated.\n")
