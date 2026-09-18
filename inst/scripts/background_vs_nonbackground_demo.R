#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mmirage)
})

args <- commandArgs(trailingOnly = TRUE)
write_data <- "--write-data" %in% args

`%||%` <- function(x, y) if (length(x) == 0 || is.na(x)) y else x

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", file_arg[1] %||% NA_character_)
script_dir <- if (is.na(script_path)) getwd() else dirname(normalizePath(script_path))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
extdata_dir <- file.path(repo_root, "inst", "extdata")
dir.create(extdata_dir, recursive = TRUE, showWarnings = FALSE)

variant_file <- file.path(extdata_dir, "background_demo_variants.csv")
truth_file <- file.path(extdata_dir, "background_demo_truth.csv")
outdir <- "background_demo_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

simulate_background_demo <- function(seed = 20260917) {
  set.seed(seed)
  n_case <- 5000
  n_control <- 5000
  variants_per_gene <- 15
  genes <- sprintf("GENE%02d", 1:30)
  gene_group <- c(
    rep("causal_risk", 6),
    rep("high_background_noncausal", 8),
    rep("ordinary_null", 16)
  )
  truth <- data.frame(
    Gene = genes,
    gene_group = gene_group,
    trait1_risk = gene_group == "causal_risk",
    high_background_gene = gene_group == "high_background_noncausal",
    stringsAsFactors = FALSE
  )

  rows <- vector("list", length(genes))
  for (i in seq_along(genes)) {
    categories <- sample(c("LoF", "Mis"), variants_per_gene, replace = TRUE,
      prob = c(0.35, 0.65))
    reported_rate <- runif(variants_per_gene, 0.00008, 0.00025)
    high_background <- gene_group[i] == "high_background_noncausal" &
      runif(variants_per_gene) < 0.75
    count_rate <- reported_rate * ifelse(high_background, 8, 1)
    is_causal_variant <- gene_group[i] == "causal_risk" &
      rbinom(variants_per_gene, 1, ifelse(categories == "LoF", 0.55, 0.35)) == 1
    rr <- ifelse(is_causal_variant, ifelse(categories == "LoF", 9.0, 6.0), 1.0)
    case_count <- rpois(variants_per_gene, 2 * n_case * count_rate * rr)
    control_count <- rpois(variants_per_gene, 2 * n_control * count_rate)
    total_count <- case_count + control_count
    p_null <- n_case / (n_case + n_control)
    p_alt <- ifelse(categories == "LoF", 5, 3) * n_case /
      (ifelse(categories == "LoF", 5, 3) * n_case + n_control)
    log_bf <- dbinom(case_count, total_count, p_alt, log = TRUE) -
      dbinom(case_count, total_count, p_null, log = TRUE)
    rows[[i]] <- data.frame(
      ID = paste0("v", seq_len(variants_per_gene), "_", genes[i]),
      Gene = genes[i],
      category = categories,
      No.case = case_count,
      No.contr = control_count,
      mutation_rate = reported_rate,
      count_rate = count_rate,
      high_background_variant = high_background,
      true_causal_variant = is_causal_variant,
      log_bf = log_bf,
      stringsAsFactors = FALSE
    )
  }
  variants_trait1 <- do.call(rbind, rows)

  variants_trait2 <- variants_trait1
  variants_trait2$ID <- paste0("trait2_", variants_trait2$ID)
  variants_trait2$No.case <- rpois(nrow(variants_trait2), 2 * n_case * variants_trait2$mutation_rate)
  variants_trait2$No.contr <- rpois(nrow(variants_trait2), 2 * n_control * variants_trait2$mutation_rate)
  variants_trait2$log_bf <- 0
  variants_trait2$high_background_variant <- FALSE
  variants_trait2$true_causal_variant <- FALSE

  list(trait1 = variants_trait1, trait2 = variants_trait2, truth = truth)
}

if (write_data || !file.exists(variant_file) || !file.exists(truth_file)) {
  demo <- simulate_background_demo()
  variants_out <- rbind(
    cbind(trait = "trait1", demo$trait1),
    cbind(trait = "trait2", demo$trait2)
  )
  write.csv(variants_out, variant_file, row.names = FALSE)
  write.csv(demo$truth, truth_file, row.names = FALSE)
}

variants <- read.csv(variant_file, stringsAsFactors = FALSE)
truth <- read.csv(truth_file, stringsAsFactors = FALSE)
trait1 <- variants[variants$trait == "trait1", ]
trait2 <- variants[variants$trait == "trait2", ]

n_case <- 5000
n_control <- 5000

fit <- mirage_two_trait(
  trait1,
  trait2,
  n1_trait1 = n_case,
  n0_trait1 = n_control,
  n1_trait2 = n_case,
  n0_trait2 = n_control,
  log_bf_col_trait1 = "log_bf",
  log_bf_col_trait2 = "log_bf",
  fixed_pi = c("00" = 0.90, "10" = 0.08, "01" = 0.01, "11" = 0.01),
  fixed_eta_trait1 = c("LoF" = 0.35, "Mis" = 0.15),
  fixed_eta_trait2 = c("LoF" = 0.35, "Mis" = 0.15),
  max_iter = 5
)

mtada_like_gene_score <- function(d) {
  agg <- aggregate(cbind(No.case, mutation_rate) ~ Gene + category, d, sum)
  gene_rows <- lapply(split(agg, agg$Gene), function(z) {
    lambda <- 2 * n_case * z$mutation_rate
    pvals <- ppois(z$No.case - 1, lambda = pmax(lambda, .Machine$double.xmin),
      lower.tail = FALSE)
    data.frame(
      Gene = z$Gene[1],
      mtada_like_score = sum(-log10(pmax(pvals, .Machine$double.xmin))),
      mtada_like_min_p = min(pvals),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, gene_rows)
}

case_control_gene_p <- function(d) {
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
  cbind(No.case, No.contr, high_background_variant, true_causal_variant) ~ Gene,
  trait1,
  sum
)
rate_summary <- aggregate(cbind(mutation_rate, count_rate) ~ Gene, trait1, sum)
gene_summary <- Reduce(function(x, y) merge(x, y, by = "Gene"), list(
  truth,
  count_summary,
  rate_summary,
  mtada_like_gene_score(trait1),
  case_control_gene_p(trait1),
  fit$posterior[, c("Gene", "PP_trait1", "PP_pleiotropy")]
))
gene_summary$case_control_ratio <- (gene_summary$No.case + 0.5) /
  (gene_summary$No.contr + 0.5)
gene_summary <- gene_summary[order(-gene_summary$mtada_like_score), ]

group_summary <- aggregate(
  cbind(No.case, No.contr, mutation_rate, count_rate, mtada_like_score,
        PP_trait1, case_control_ratio) ~ gene_group,
  gene_summary,
  mean
)

write.csv(gene_summary, file.path(outdir, "background_demo_gene_scores.csv"), row.names = FALSE)
write.csv(group_summary, file.path(outdir, "background_demo_group_summary.csv"), row.names = FALSE)

png(file.path(outdir, "background_demo_score_comparison.png"), width = 1500, height = 1100, res = 180)
cols <- c(
  causal_risk = "#1f77b4",
  high_background_noncausal = "#d95f02",
  ordinary_null = "#7f7f7f"
)
plot(
  gene_summary$mtada_like_score,
  gene_summary$PP_trait1,
  pch = 19,
  col = cols[gene_summary$gene_group],
  xlab = "mTADA-like aggregated score (case counts only)",
  ylab = "MIRAGE posterior probability (case-control)",
  main = "High-background genes look strong after aggregation but not in MIRAGE"
)
legend("topleft", legend = names(cols), col = cols, pch = 19, bty = "n")
text(
  gene_summary$mtada_like_score,
  gene_summary$PP_trait1,
  labels = gene_summary$Gene,
  pos = 4,
  cex = 0.65
)
dev.off()

cat("\nBackground vs non-background demo complete.\n")
cat("Simulation data:\n")
cat("  ", normalizePath(variant_file), "\n", sep = "")
cat("  ", normalizePath(truth_file), "\n", sep = "")
cat("Output directory:\n")
cat("  ", normalizePath(outdir), "\n\n", sep = "")

cat("Group-level averages:\n")
print(group_summary)

cat("\nTop genes by mTADA-like aggregation (case counts only):\n")
print(head(gene_summary[, c(
  "Gene", "gene_group", "No.case", "No.contr", "high_background_variant",
  "true_causal_variant", "mtada_like_score", "burden_case_control_p", "PP_trait1"
)], 10))

cat("\nTop genes by MIRAGE posterior (case-control counts):\n")
mirage_rank <- gene_summary[order(-gene_summary$PP_trait1), ]
print(head(mirage_rank[, c(
  "Gene", "gene_group", "No.case", "No.contr", "high_background_variant",
  "true_causal_variant", "mtada_like_score", "burden_case_control_p", "PP_trait1"
)], 10))
