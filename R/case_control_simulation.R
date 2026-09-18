simulate_case_control_two_trait <- function(
    n_genes = 500,
    variants_per_gene = 15,
    pi = c("00" = 0.94, "10" = 0.04, "01" = 0.01, "11" = 0.01),
    eta_trait1 = c("C1" = 0.05, "C2" = 0.20, "C3" = 0.50),
    eta_trait2 = c("C1" = 0.03, "C2" = 0.15, "C3" = 0.40),
    category_prob = c("C1" = 0.60, "C2" = 0.30, "C3" = 0.10),
    rr_mean = c("C1" = 3, "C2" = 3, "C3" = 5),
    rr_sigma = c("C1" = 1, "C2" = 1, "C3" = 1),
    n_case_trait1 = 3000,
    n_control_trait1 = 3000,
    n_case_trait2 = 3000,
    n_control_trait2 = 3000,
    background_prob = 0.05,
    background_multiplier = 5,
    max_count_rate = 0.05,
    fast_log_bf = TRUE,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pi <- normalize_pi(pi)
  categories <- names(category_prob)
  eta1 <- expand_by_category(eta_trait1, categories, "eta_trait1")
  eta2 <- expand_by_category(eta_trait2, categories, "eta_trait2")
  rr_mean <- expand_by_category(rr_mean, categories, "rr_mean")
  rr_sigma <- expand_by_category(rr_sigma, categories, "rr_sigma")
  category_prob <- category_prob / sum(category_prob)

  genes <- sprintf("GENE%05d", seq_len(n_genes))
  states <- sample(names(pi), n_genes, replace = TRUE, prob = pi)
  active1 <- states %in% c("10", "11")
  active2 <- states %in% c("01", "11")

  make_trait <- function(trait, active, eta, n_case, n_control) {
    rows <- vector("list", n_genes)
    for (i in seq_len(n_genes)) {
      cats <- sample(categories, variants_per_gene, replace = TRUE, prob = category_prob)
      q <- pmin(rbeta(variants_per_gene, shape1 = 0.35, shape2 = 80) * 0.04, 0.01)
      high_background <- runif(variants_per_gene) < background_prob
      count_rate <- pmin(q * ifelse(high_background, background_multiplier, 1), max_count_rate)
      z <- if (active[i]) rbinom(variants_per_gene, 1, eta[cats]) else rep(0L, variants_per_gene)
      rr <- rep(1, variants_per_gene)
      idx <- z == 1
      if (any(idx)) {
        rr[idx] <- rgamma(
          sum(idx),
          shape = rr_mean[cats[idx]] * rr_sigma[cats[idx]],
          rate = rr_sigma[cats[idx]]
        )
      }
      case_count <- rpois(variants_per_gene, 2 * n_case * count_rate * rr)
      control_count <- rpois(variants_per_gene, 2 * n_control * count_rate)
      total_count <- case_count + control_count
      log_bf <- rep(NA_real_, variants_per_gene)
      if (fast_log_bf) {
        p_null <- n_case / (n_case + n_control)
        p_alt <- rr_mean[cats] * n_case / (rr_mean[cats] * n_case + n_control)
        log_bf <- dbinom(case_count, total_count, p_alt, log = TRUE) -
          dbinom(case_count, total_count, p_null, log = TRUE)
      }
      rows[[i]] <- data.frame(
        ID = paste0("t", trait, "_", genes[i], "_v", seq_len(variants_per_gene)),
        Gene = genes[i],
        No.case = case_count,
        No.contr = control_count,
        category = cats,
        log_bf = log_bf,
        mutation_rate = q,
        count_rate = count_rate,
        high_background = high_background,
        true_z = z,
        true_rr = rr,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }

  list(
    trait1 = make_trait(1, active1, eta1, n_case_trait1, n_control_trait1),
    trait2 = make_trait(2, active2, eta2, n_case_trait2, n_control_trait2),
    truth = data.frame(
      Gene = genes,
      state = states,
      trait1_risk = active1,
      trait2_risk = active2,
      pleiotropic = states == "11",
      stringsAsFactors = FALSE
    ),
    parameters = list(
      pi = pi,
      eta_trait1 = eta1,
      eta_trait2 = eta2,
      category_prob = category_prob,
      rr_mean = rr_mean,
      rr_sigma = rr_sigma,
      n_case_trait1 = n_case_trait1,
      n_control_trait1 = n_control_trait1,
      n_case_trait2 = n_case_trait2,
      n_control_trait2 = n_control_trait2,
      background_prob = background_prob,
      background_multiplier = background_multiplier
    )
  )
}

gene_burden_test <- function(data, n_case, n_control, gene_col = "Gene",
                             case_col = "No.case", control_col = "No.contr",
                             variants_per_gene = NULL) {
  required <- c(gene_col, case_col, control_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("data are missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  d <- data.frame(
    Gene = as.character(data[[gene_col]]),
    case = as.numeric(data[[case_col]]),
    control = as.numeric(data[[control_col]]),
    stringsAsFactors = FALSE
  )
  agg <- aggregate(cbind(case, control) ~ Gene, d, sum)
  if (is.null(variants_per_gene)) {
    variant_counts <- aggregate(case ~ Gene, d, length)
    names(variant_counts)[2] <- "variant_count"
    agg <- merge(agg, variant_counts, by = "Gene", all.x = TRUE)
  } else {
    agg$variant_count <- variants_per_gene
  }
  p_value <- odds_ratio <- numeric(nrow(agg))
  for (i in seq_len(nrow(agg))) {
    case_an <- 2 * n_case * agg$variant_count[i]
    control_an <- 2 * n_control * agg$variant_count[i]
    tab <- matrix(
      c(
        agg$case[i], case_an - agg$case[i],
        agg$control[i], control_an - agg$control[i]
      ),
      nrow = 2,
      byrow = TRUE
    )
    ft <- fisher.test(tab, alternative = "two.sided")
    p_value[i] <- ft$p.value
    odds_ratio[i] <- unname(ft$estimate)
  }
  agg$odds_ratio <- odds_ratio
  agg$p_value <- p_value
  agg$bh_fdr <- p.adjust(p_value, method = "BH")
  agg$score <- -log10(pmax(p_value, .Machine$double.xmin))
  agg
}
