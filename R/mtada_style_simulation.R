# mTADA-style simulation and benchmarking for the two-trait MIRAGE model.
#
# The generator follows the four-state, gene-level de novo model described by
# Nguyen et al. (2020). It uses synthetic mutation rates only.

mtada_state_probabilities <- function(pi_trait1 = 0.05, pi_trait2 = 0.03,
                                      pi_pleiotropic = 0.01) {
  if (pi_trait1 < 0 || pi_trait2 < 0 ||
      pi_pleiotropic < 0 ||
      pi_pleiotropic > min(pi_trait1, pi_trait2) ||
      pi_trait1 + pi_trait2 - pi_pleiotropic > 1) {
    stop("Invalid single-trait or pleiotropic risk proportions.", call. = FALSE)
  }
  c(
    "00" = 1 - pi_trait1 - pi_trait2 + pi_pleiotropic,
    "10" = pi_trait1 - pi_pleiotropic,
    "01" = pi_trait2 - pi_pleiotropic,
    "11" = pi_pleiotropic
  )
}

simulate_synthetic_mutation_rates <- function(n_genes, categories,
                                              median_rates = NULL,
                                              log_sd = 0.9) {
  if (is.null(median_rates)) {
    median_rates <- rep(1e-5, length(categories))
    names(median_rates) <- categories
  }
  median_rates <- expand_by_category(
    median_rates, categories, "median_rates"
  )
  rates <- vapply(
    categories,
    function(category) {
      rlnorm(n_genes, log(median_rates[[category]]), log_sd)
    },
    numeric(n_genes)
  )
  colnames(rates) <- categories
  rates
}

log_bf_denovo_tada <- function(count, n_trios, mutation_rate,
                               gamma_mean, beta) {
  if (any(count < 0) || any(mutation_rate <= 0) ||
      n_trios <= 0 || gamma_mean <= 0 || beta <= 0) {
    stop("Counts and de novo model parameters are outside their valid ranges.",
         call. = FALSE)
  }
  lambda <- 2 * n_trios * mutation_rate
  dnbinom(
    count,
    size = gamma_mean * beta,
    prob = beta / (beta + lambda),
    log = TRUE
  ) - dpois(count, lambda = lambda, log = TRUE)
}

simulate_mtada_style <- function(
    n_genes = 19538,
    categories = c("MiD", "LoF"),
    mutation_rates = NULL,
    median_mutation_rates = c("MiD" = 1.5e-5, "LoF" = 0.5e-5),
    mutation_rate_log_sd = 0.9,
    pi_trait1 = 0.05,
    pi_trait2 = 0.03,
    pi_pleiotropic = 0.01,
    n_trios_trait1 = 5122,
    n_trios_trait2 = 1077,
    gamma_mean_trait1 = c("MiD" = 20, "LoF" = 50),
    gamma_mean_trait2 = c("MiD" = 12, "LoF" = 2),
    beta_trait1 = c("MiD" = 0.8, "LoF" = 0.8),
    beta_trait2 = c("MiD" = 0.8, "LoF" = 0.8),
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pi <- mtada_state_probabilities(
    pi_trait1, pi_trait2, pi_pleiotropic
  )
  gamma1 <- expand_by_category(
    gamma_mean_trait1, categories, "gamma_mean_trait1"
  )
  gamma2 <- expand_by_category(
    gamma_mean_trait2, categories, "gamma_mean_trait2"
  )
  beta1 <- expand_by_category(beta_trait1, categories, "beta_trait1")
  beta2 <- expand_by_category(beta_trait2, categories, "beta_trait2")

  if (is.null(mutation_rates)) {
    mutation_rates <- simulate_synthetic_mutation_rates(
      n_genes,
      categories,
      median_mutation_rates,
      mutation_rate_log_sd
    )
  } else {
    mutation_rates <- as.matrix(mutation_rates)
    if (nrow(mutation_rates) != n_genes ||
        !all(categories %in% colnames(mutation_rates))) {
      stop("mutation_rates must have one row per gene and named category columns.",
           call. = FALSE)
    }
    mutation_rates <- mutation_rates[, categories, drop = FALSE]
  }
  if (any(!is.finite(mutation_rates)) || any(mutation_rates <= 0)) {
    stop("mutation_rates must be finite and positive.", call. = FALSE)
  }

  genes <- sprintf("GENE%05d", seq_len(n_genes))
  states <- sample(names(pi), n_genes, replace = TRUE, prob = pi)
  active1 <- states %in% c("10", "11")
  active2 <- states %in% c("01", "11")

  simulate_trait <- function(trait, active, n_trios, gamma_mean, beta) {
    rows <- vector("list", length(categories))
    for (j in seq_along(categories)) {
      category <- categories[j]
      rr <- rep(1, n_genes)
      rr[active] <- rgamma(
        sum(active),
        shape = gamma_mean[[category]] * beta[[category]],
        rate = beta[[category]]
      )
      mu <- mutation_rates[, category]
      count <- rpois(n_genes, 2 * n_trios * mu * rr)
      rows[[j]] <- data.frame(
        ID = paste0("t", trait, "_", genes, "_", category),
        Gene = genes,
        No.case = count,
        No.contr = 0,
        category = category,
        log_bf = log_bf_denovo_tada(
          count, n_trios, mu,
          gamma_mean[[category]], beta[[category]]
        ),
        mutation_rate = mu,
        true_rr = rr,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }

  list(
    trait1 = simulate_trait(1, active1, n_trios_trait1, gamma1, beta1),
    trait2 = simulate_trait(2, active2, n_trios_trait2, gamma2, beta2),
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
      pi_trait1 = pi_trait1,
      pi_trait2 = pi_trait2,
      pi_pleiotropic = pi_pleiotropic,
      n_trios_trait1 = n_trios_trait1,
      n_trios_trait2 = n_trios_trait2,
      gamma_mean_trait1 = gamma1,
      gamma_mean_trait2 = gamma2,
      beta_trait1 = beta1,
      beta_trait2 = beta2
    )
  )
}

binary_auc <- function(score, truth) {
  truth <- as.logical(truth)
  if (!any(truth) || all(truth)) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  n1 <- sum(truth)
  n0 <- sum(!truth)
  (sum(ranks[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

selection_metrics <- function(score, truth, threshold = 0.8) {
  selected <- score >= threshold
  positives <- sum(truth)
  discoveries <- sum(selected)
  data.frame(
    threshold = threshold,
    discoveries = discoveries,
    true_positives = sum(selected & truth),
    power = if (positives > 0) sum(selected & truth) / positives else NA_real_,
    fdp = if (discoveries > 0) sum(selected & !truth) / discoveries else 0,
    stringsAsFactors = FALSE
  )
}

bfdr_metrics <- function(score, truth, target = 0.05) {
  ord <- order(score, decreasing = TRUE)
  estimated <- cumsum(1 - score[ord]) / seq_along(score)
  keep <- which(estimated <= target)
  selected <- if (length(keep)) ord[seq_len(max(keep))] else integer()
  data.frame(
    bfdr_target = target,
    discoveries = length(selected),
    true_positives = sum(truth[selected]),
    power = if (sum(truth) > 0) sum(truth[selected]) / sum(truth) else NA_real_,
    observed_fdp = if (length(selected)) mean(!truth[selected]) else 0,
    estimated_bfdr = if (length(selected)) estimated[length(selected)] else 0,
    stringsAsFactors = FALSE
  )
}

fit_mtada_style_mirage <- function(sim, estimate_pi = TRUE,
                                   pi_init = NULL, max_iter = 200) {
  if (is.null(pi_init)) {
    pi_init <- sim$parameters$pi
  }
  fit <- mirage_two_trait(
    sim$trait1,
    sim$trait2,
    n1_trait1 = sim$parameters$n_trios_trait1,
    n0_trait1 = sim$parameters$n_trios_trait1,
    n1_trait2 = sim$parameters$n_trios_trait2,
    n0_trait2 = sim$parameters$n_trios_trait2,
    log_bf_col_trait1 = "log_bf",
    log_bf_col_trait2 = "log_bf",
    fixed_eta_trait1 = setNames(
      rep(1, length(sim$parameters$gamma_mean_trait1)),
      names(sim$parameters$gamma_mean_trait1)
    ),
    fixed_eta_trait2 = setNames(
      rep(1, length(sim$parameters$gamma_mean_trait2)),
      names(sim$parameters$gamma_mean_trait2)
    ),
    pi_init = pi_init,
    fixed_pi = if (estimate_pi) NULL else sim$parameters$pi,
    max_iter = max_iter
  )

  truth <- sim$truth[match(fit$posterior$Gene, sim$truth$Gene), ]
  separate1 <- single_trait_posterior_from_log_bf(
    fit$gene_bf$log_B_trait1, sim$parameters$pi_trait1
  )
  separate2 <- single_trait_posterior_from_log_bf(
    fit$gene_bf$log_B_trait2, sim$parameters$pi_trait2
  )

  auc <- data.frame(
    target = c("trait1", "trait1", "trait2", "trait2", "pleiotropy"),
    method = c("joint", "separate", "joint", "separate", "joint"),
    auc = c(
      binary_auc(fit$posterior$PP_trait1, truth$trait1_risk),
      binary_auc(separate1, truth$trait1_risk),
      binary_auc(fit$posterior$PP_trait2, truth$trait2_risk),
      binary_auc(separate2, truth$trait2_risk),
      binary_auc(fit$posterior$PP_pleiotropy, truth$pleiotropic)
    ),
    stringsAsFactors = FALSE
  )

  score_sets <- list(
    joint_trait1 = list(fit$posterior$PP_trait1, truth$trait1_risk),
    separate_trait1 = list(separate1, truth$trait1_risk),
    joint_trait2 = list(fit$posterior$PP_trait2, truth$trait2_risk),
    separate_trait2 = list(separate2, truth$trait2_risk),
    joint_pleiotropy = list(
      fit$posterior$PP_pleiotropy, truth$pleiotropic
    )
  )
  threshold <- do.call(
    rbind,
    lapply(names(score_sets), function(name) {
      out <- selection_metrics(
        score_sets[[name]][[1]], score_sets[[name]][[2]], 0.8
      )
      out$analysis <- name
      out
    })
  )
  rownames(threshold) <- NULL

  bfdr <- do.call(
    rbind,
    lapply(names(score_sets), function(name) {
      out <- bfdr_metrics(
        score_sets[[name]][[1]], score_sets[[name]][[2]], 0.05
      )
      out$analysis <- name
      out
    })
  )
  rownames(bfdr) <- NULL

  list(
    simulation = sim,
    fit = fit,
    separate_posterior = data.frame(
      Gene = fit$posterior$Gene,
      PP_trait1 = separate1,
      PP_trait2 = separate2,
      stringsAsFactors = FALSE
    ),
    metrics = list(
      auc = auc,
      pp_threshold = threshold,
      bfdr = bfdr,
      pi = data.frame(
        state = names(sim$parameters$pi),
        truth = as.numeric(sim$parameters$pi),
        estimate = as.numeric(
          fit$parameters$pi[names(sim$parameters$pi)]
        ),
        stringsAsFactors = FALSE
      )
    )
  )
}

run_mtada_style_grid <- function(
    pi3_values = c(0, 0.01, 0.02, 0.03),
    n_replicates = 10,
    n_genes = 2000,
    seed = 20260619,
    ...) {
  rows <- vector("list", length(pi3_values) * n_replicates)
  k <- 0
  for (pi3 in pi3_values) {
    for (replicate_id in seq_len(n_replicates)) {
      k <- k + 1
      sim <- simulate_mtada_style(
        n_genes = n_genes,
        pi_pleiotropic = pi3,
        seed = seed + k,
        ...
      )
      result <- fit_mtada_style_mirage(sim)
      auc <- result$metrics$auc
      bfdr <- result$metrics$bfdr
      rows[[k]] <- rbind(
        data.frame(
          pi3 = pi3,
          replicate = replicate_id,
          analysis = paste(auc$method, auc$target, sep = "_"),
          metric = "auc",
          value = auc$auc,
          stringsAsFactors = FALSE
        ),
        data.frame(
          pi3 = pi3,
          replicate = replicate_id,
          analysis = bfdr$analysis,
          metric = "bfdr_observed_fdp",
          value = bfdr$observed_fdp,
          stringsAsFactors = FALSE
        ),
        data.frame(
          pi3 = pi3,
          replicate = replicate_id,
          analysis = bfdr$analysis,
          metric = "bfdr_power",
          value = bfdr$power,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  do.call(rbind, rows)
}
