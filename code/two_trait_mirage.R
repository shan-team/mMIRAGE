# Two-trait MIRAGE gene-level four-state model.
#
# This first implementation keeps pleiotropy at the gene level. Variant risk
# probabilities and relative-risk priors are trait specific.

log_sum_exp <- function(x) {
  if (length(x) == 0) return(-Inf)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

row_log_sum_exp <- function(mat) {
  apply(mat, 1, log_sum_exp)
}

safe_log_eta_mix <- function(log_bf, eta) {
  if (any(!is.finite(log_bf))) {
    stop("log_bf values must be finite.", call. = FALSE)
  }
  if (!is.numeric(eta) || length(eta) != 1 || is.na(eta) || eta < 0 || eta > 1) {
    stop("eta values must be numeric probabilities in [0, 1].", call. = FALSE)
  }
  if (eta == 0) return(rep(0, length(log_bf)))
  if (eta == 1) return(log_bf)
  vapply(log_bf, function(lb) log_sum_exp(c(log1p(-eta), log(eta) + lb)), numeric(1))
}

standardize_trait_data <- function(data, trait_name, gene_col, id_col, case_col,
                                   control_col, category_col, log_bf_col = NULL) {
  if (!is.data.frame(data)) {
    stop(sprintf("%s data must be a data frame.", trait_name), call. = FALSE)
  }
  required <- c(gene_col, case_col, control_col, category_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(sprintf("%s data are missing required columns: %s",
                 trait_name, paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (!is.null(log_bf_col) && !log_bf_col %in% names(data)) {
    stop(sprintf("%s log_bf_col '%s' is not present.", trait_name, log_bf_col),
         call. = FALSE)
  }
  out <- data.frame(
    ID = if (id_col %in% names(data)) as.character(data[[id_col]]) else paste0(trait_name, "_var_", seq_len(nrow(data))),
    Gene = as.character(data[[gene_col]]),
    No.case = as.numeric(data[[case_col]]),
    No.contr = as.numeric(data[[control_col]]),
    category = as.character(data[[category_col]]),
    stringsAsFactors = FALSE
  )
  if (any(!nzchar(out$Gene)) || any(is.na(out$Gene))) {
    stop(sprintf("%s data contain missing gene identifiers.", trait_name), call. = FALSE)
  }
  if (any(is.na(out$No.case)) || any(is.na(out$No.contr)) ||
      any(out$No.case < 0) || any(out$No.contr < 0)) {
    stop(sprintf("%s counts must be nonnegative numeric values.", trait_name), call. = FALSE)
  }
  if (any(out$No.case != floor(out$No.case)) || any(out$No.contr != floor(out$No.contr))) {
    stop(sprintf("%s counts must be integer-valued.", trait_name), call. = FALSE)
  }
  if (any(!nzchar(out$category)) || any(is.na(out$category))) {
    stop(sprintf("%s data contain missing annotation categories.", trait_name), call. = FALSE)
  }
  out$log_bf <- if (!is.null(log_bf_col)) as.numeric(data[[log_bf_col]]) else NA_real_
  if (!is.null(log_bf_col) && any(!is.finite(out$log_bf))) {
    stop(sprintf("%s log BF values must be finite.", trait_name), call. = FALSE)
  }
  out
}

expand_by_category <- function(value, categories, name) {
  if (is.null(value)) return(NULL)
  if (length(value) == 1 && is.null(names(value))) {
    out <- rep(as.numeric(value), length(categories))
    names(out) <- categories
    return(out)
  }
  if (!is.null(names(value))) {
    missing <- setdiff(categories, names(value))
    if (length(missing) > 0) {
      stop(sprintf("%s is missing categories: %s", name, paste(missing, collapse = ", ")),
           call. = FALSE)
    }
    out <- as.numeric(value[categories])
    names(out) <- categories
    return(out)
  }
  if (length(value) != length(categories)) {
    stop(sprintf("%s must be scalar, named by category, or length %d.",
                 name, length(categories)), call. = FALSE)
  }
  out <- as.numeric(value)
  names(out) <- categories
  out
}

log_variant_bf_integrated <- function(var_case, var_control, gamma_mean, sigma,
                                      n_case, n_control, upper = 100) {
  if (gamma_mean <= 0 || sigma <= 0) {
    stop("gamma and sigma must be positive.", call. = FALSE)
  }
  total <- var_case + var_control
  log_null <- dbinom(var_case, total, n_case / (n_case + n_control), log = TRUE)
  integrand <- function(rr) {
    exp(dbinom(var_case, total, rr * n_case / (rr * n_case + n_control), log = TRUE) +
          dgamma(rr, shape = gamma_mean * sigma, rate = sigma, log = TRUE))
  }
  risk <- integrate(integrand, lower = 0, upper = upper, stop.on.error = FALSE)$value
  if (!is.finite(risk) || risk <= 0) {
    return(-Inf)
  }
  log(risk) - log_null
}

add_variant_log_bf <- function(trait_data, n_case, n_control, gamma, sigma,
                               cache = TRUE) {
  categories <- unique(trait_data$category)
  gamma_by_cat <- expand_by_category(gamma, categories, "gamma")
  sigma_by_cat <- expand_by_category(sigma, categories, "sigma")
  if (all(is.finite(trait_data$log_bf))) return(trait_data)

  cache_env <- new.env(parent = emptyenv())
  log_bf <- numeric(nrow(trait_data))
  for (i in seq_len(nrow(trait_data))) {
    cat_i <- trait_data$category[i]
    key <- paste(trait_data$No.case[i], trait_data$No.contr[i], n_case, n_control,
                 gamma_by_cat[[cat_i]], sigma_by_cat[[cat_i]], sep = "|")
    if (cache && exists(key, envir = cache_env, inherits = FALSE)) {
      log_bf[i] <- get(key, envir = cache_env)
    } else {
      log_bf[i] <- log_variant_bf_integrated(
        trait_data$No.case[i], trait_data$No.contr[i],
        gamma_by_cat[[cat_i]], sigma_by_cat[[cat_i]],
        n_case, n_control
      )
      if (cache) assign(key, log_bf[i], envir = cache_env)
    }
  }
  trait_data$log_bf <- log_bf
  trait_data
}

align_trait_data <- function(trait_data, genes) {
  trait_data$gene_id <- match(trait_data$Gene, genes)
  trait_data$category_id <- match(trait_data$category, unique(trait_data$category))
  trait_data
}

compute_log_gene_bf <- function(trait_data, genes, eta) {
  log_B <- rep(0, length(genes))
  names(log_B) <- genes
  if (nrow(trait_data) == 0) return(log_B)
  eta_by_variant <- eta[trait_data$category]
  contrib <- numeric(nrow(trait_data))
  for (i in seq_len(nrow(trait_data))) {
    contrib[i] <- safe_log_eta_mix(trait_data$log_bf[i], eta_by_variant[i])
  }
  sums <- rowsum(contrib, trait_data$gene_id, reorder = FALSE)
  log_B[as.integer(rownames(sums))] <- as.numeric(sums[, 1])
  log_B
}

compute_rho <- function(trait_data, eta) {
  if (nrow(trait_data) == 0) return(numeric())
  eta_by_variant <- eta[trait_data$category]
  out <- numeric(nrow(trait_data))
  for (i in seq_len(nrow(trait_data))) {
    e <- eta_by_variant[i]
    if (e == 0) {
      out[i] <- 0
    } else if (e == 1) {
      out[i] <- 1
    } else {
      log_den <- log_sum_exp(c(log1p(-e), log(e) + trait_data$log_bf[i]))
      out[i] <- exp(log(e) + trait_data$log_bf[i] - log_den)
    }
  }
  out
}

compute_joint_quantities <- function(log_B1, log_B2, pi) {
  log_pi <- log(pi)
  log_weights <- cbind(
    `00` = log_pi[["00"]],
    `10` = log_pi[["10"]] + log_B1,
    `01` = log_pi[["01"]] + log_B2,
    `11` = log_pi[["11"]] + log_B1 + log_B2
  )
  log_den <- row_log_sum_exp(log_weights)
  tau <- exp(log_weights - log_den)
  colnames(tau) <- c("00", "10", "01", "11")
  list(tau = tau, loglik = sum(log_den), log_den = log_den)
}

update_eta_trait <- function(trait_data, genes, w, rho, eta_old, min_denom = 1e-12) {
  eta_new <- eta_old
  diagnostics <- data.frame(category = names(eta_old), denominator = NA_real_,
                            weak = FALSE, stringsAsFactors = FALSE)
  for (cat in names(eta_old)) {
    idx <- which(trait_data$category == cat)
    if (length(idx) == 0) {
      diagnostics$denominator[diagnostics$category == cat] <- 0
      diagnostics$weak[diagnostics$category == cat] <- TRUE
      next
    }
    denom <- sum(w[trait_data$gene_id[idx]])
    numer <- sum(w[trait_data$gene_id[idx]] * rho[idx])
    diagnostics$denominator[diagnostics$category == cat] <- denom
    diagnostics$weak[diagnostics$category == cat] <- denom <= min_denom
    if (denom > min_denom) {
      eta_new[cat] <- min(max(numer / denom, 0), 1)
    }
  }
  list(eta = eta_new, diagnostics = diagnostics)
}

bfdr_table <- function(genes, pp, target = "target") {
  ord <- order(pp, decreasing = TRUE)
  pp_sorted <- pp[ord]
  data.frame(
    Gene = genes[ord],
    target = target,
    posterior_probability = pp_sorted,
    local_fdr = 1 - pp_sorted,
    bayesian_fdr = cumsum(1 - pp_sorted) / seq_along(pp_sorted),
    rank = seq_along(pp_sorted),
    stringsAsFactors = FALSE
  )
}

single_trait_posterior_from_log_bf <- function(log_B, delta) {
  log_num <- log(delta) + log_B
  log_den <- vapply(log_B, function(lb) log_sum_exp(c(log1p(-delta), log(delta) + lb)), numeric(1))
  exp(log_num - log_den)
}

normalize_pi <- function(pi) {
  if (is.null(names(pi))) names(pi) <- c("00", "10", "01", "11")
  pi <- pi[c("00", "10", "01", "11")]
  if (any(is.na(pi)) || any(pi < 0) || sum(pi) <= 0) {
    stop("pi values must be nonnegative and have positive sum.", call. = FALSE)
  }
  pi / sum(pi)
}

validate_eta <- function(eta, categories, name, default = 0.1) {
  if (is.null(eta)) eta <- default
  out <- expand_by_category(eta, categories, name)
  if (any(is.na(out)) || any(out < 0) || any(out > 1)) {
    stop(sprintf("%s values must be probabilities in [0, 1].", name), call. = FALSE)
  }
  out
}

mirage_two_trait <- function(
    data_trait1,
    data_trait2,
    n1_trait1,
    n0_trait1,
    n1_trait2,
    n0_trait2,
    gamma_trait1 = 3,
    sigma_trait1 = 2,
    gamma_trait2 = 3,
    sigma_trait2 = 2,
    trait_names = c("ASD", "SCZ"),
    gene_col = "Gene",
    id_col = "ID",
    case_col = "No.case",
    control_col = "No.contr",
    category_col = "category",
    log_bf_col_trait1 = NULL,
    log_bf_col_trait2 = NULL,
    gene_universe = c("union", "intersection", "trait1", "trait2"),
    estimate_pi = TRUE,
    estimate_eta = TRUE,
    pi_init = NULL,
    eta_init_trait1 = NULL,
    eta_init_trait2 = NULL,
    fixed_pi = NULL,
    fixed_eta_trait1 = NULL,
    fixed_eta_trait2 = NULL,
    dirichlet_alpha = NULL,
    cache_variant_bf = TRUE,
    tolerance_loglik = 1e-8,
    tolerance_parameter = 1e-6,
    max_iter = 500,
    verbose = FALSE) {
  gene_universe <- match.arg(gene_universe)
  if (length(trait_names) != 2) stop("trait_names must have length 2.", call. = FALSE)
  if (max_iter < 1) stop("max_iter must be positive.", call. = FALSE)

  d1 <- standardize_trait_data(data_trait1, trait_names[1], gene_col, id_col,
                               case_col, control_col, category_col, log_bf_col_trait1)
  d2 <- standardize_trait_data(data_trait2, trait_names[2], gene_col, id_col,
                               case_col, control_col, category_col, log_bf_col_trait2)
  d1 <- add_variant_log_bf(d1, n1_trait1, n0_trait1, gamma_trait1, sigma_trait1, cache_variant_bf)
  d2 <- add_variant_log_bf(d2, n1_trait2, n0_trait2, gamma_trait2, sigma_trait2, cache_variant_bf)

  genes <- switch(
    gene_universe,
    union = sort(unique(c(d1$Gene, d2$Gene))),
    intersection = sort(intersect(unique(d1$Gene), unique(d2$Gene))),
    trait1 = sort(unique(d1$Gene)),
    trait2 = sort(unique(d2$Gene))
  )
  if (length(genes) == 0) stop("No genes remain after gene alignment.", call. = FALSE)
  d1 <- d1[d1$Gene %in% genes, , drop = FALSE]
  d2 <- d2[d2$Gene %in% genes, , drop = FALSE]
  d1 <- align_trait_data(d1, genes)
  d2 <- align_trait_data(d2, genes)

  cat1 <- unique(d1$category)
  cat2 <- unique(d2$category)
  eta1 <- if (!is.null(fixed_eta_trait1)) {
    validate_eta(fixed_eta_trait1, cat1, "fixed_eta_trait1")
  } else {
    validate_eta(eta_init_trait1, cat1, "eta_init_trait1")
  }
  eta2 <- if (!is.null(fixed_eta_trait2)) {
    validate_eta(fixed_eta_trait2, cat2, "fixed_eta_trait2")
  } else {
    validate_eta(eta_init_trait2, cat2, "eta_init_trait2")
  }
  if (!is.null(fixed_eta_trait1) || !is.null(fixed_eta_trait2)) estimate_eta <- FALSE

  pi <- if (!is.null(fixed_pi)) {
    normalize_pi(fixed_pi)
  } else if (!is.null(pi_init)) {
    normalize_pi(pi_init)
  } else {
    normalize_pi(c("00" = 0.90, "10" = 0.04, "01" = 0.04, "11" = 0.02))
  }
  if (!is.null(fixed_pi)) estimate_pi <- FALSE
  if (!is.null(dirichlet_alpha)) {
    dirichlet_alpha <- normalize_pi(dirichlet_alpha) * sum(dirichlet_alpha)
    if (any(dirichlet_alpha < 1)) {
      stop("dirichlet_alpha entries must be >= 1 for the MAP update.", call. = FALSE)
    }
  }

  trace_loglik <- numeric()
  trace_change <- numeric()
  trace_pi <- list()
  trace_eta1 <- list()
  trace_eta2 <- list()
  monotone <- TRUE
  converged <- FALSE
  eta_diag1 <- eta_diag2 <- NULL

  for (iter in seq_len(max_iter)) {
    log_B1 <- compute_log_gene_bf(d1, genes, eta1)
    log_B2 <- compute_log_gene_bf(d2, genes, eta2)
    q_old <- compute_joint_quantities(log_B1, log_B2, pi)
    tau <- q_old$tau
    w1 <- tau[, "10"] + tau[, "11"]
    w2 <- tau[, "01"] + tau[, "11"]
    rho1 <- compute_rho(d1, eta1)
    rho2 <- compute_rho(d2, eta2)

    pi_new <- pi
    if (estimate_pi) {
      if (is.null(dirichlet_alpha)) {
        pi_new <- colMeans(tau)
      } else {
        numer <- colSums(tau) + dirichlet_alpha - 1
        pi_new <- numer / (length(genes) + sum(dirichlet_alpha) - 4)
      }
      pi_new <- normalize_pi(pi_new)
    }

    eta1_new <- eta1
    eta2_new <- eta2
    if (estimate_eta) {
      upd1 <- update_eta_trait(d1, genes, w1, rho1, eta1)
      upd2 <- update_eta_trait(d2, genes, w2, rho2, eta2)
      eta1_new <- upd1$eta
      eta2_new <- upd2$eta
      eta_diag1 <- upd1$diagnostics
      eta_diag2 <- upd2$diagnostics
    }

    log_B1_new <- compute_log_gene_bf(d1, genes, eta1_new)
    log_B2_new <- compute_log_gene_bf(d2, genes, eta2_new)
    q_new <- compute_joint_quantities(log_B1_new, log_B2_new, pi_new)
    change <- max(
      max(abs(pi_new - pi)),
      if (length(eta1)) max(abs(eta1_new - eta1)) else 0,
      if (length(eta2)) max(abs(eta2_new - eta2)) else 0
    )
    rel_ll <- abs(q_new$loglik - q_old$loglik) / (1 + abs(q_old$loglik))
    if (q_new$loglik + 1e-8 < q_old$loglik) monotone <- FALSE

    trace_loglik[iter] <- q_new$loglik
    trace_change[iter] <- change
    trace_pi[[iter]] <- pi_new
    trace_eta1[[iter]] <- eta1_new
    trace_eta2[[iter]] <- eta2_new

    pi <- pi_new
    eta1 <- eta1_new
    eta2 <- eta2_new
    if (rel_ll < tolerance_loglik && change < tolerance_parameter) {
      converged <- TRUE
      break
    }
  }

  log_B1 <- compute_log_gene_bf(d1, genes, eta1)
  log_B2 <- compute_log_gene_bf(d2, genes, eta2)
  q <- compute_joint_quantities(log_B1, log_B2, pi)
  tau <- q$tau
  colnames(tau) <- paste0("tau_", colnames(tau))
  pp_trait1 <- tau[, "tau_10"] + tau[, "tau_11"]
  pp_trait2 <- tau[, "tau_01"] + tau[, "tau_11"]
  pp_pleio <- tau[, "tau_11"]

  has1 <- genes %in% unique(d1$Gene)
  has2 <- genes %in% unique(d2$Gene)
  nvar1 <- tabulate(match(d1$Gene, genes), nbins = length(genes))
  nvar2 <- tabulate(match(d2$Gene, genes), nbins = length(genes))

  posterior <- data.frame(
    Gene = genes,
    tau,
    PP_trait1 = pp_trait1,
    PP_trait2 = pp_trait2,
    PP_pleiotropy = pp_pleio,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  out <- list(
    call = match.call(),
    trait_names = trait_names,
    posterior = posterior,
    gene_bf = data.frame(
      Gene = genes,
      log_B_trait1 = as.numeric(log_B1),
      log_B_trait2 = as.numeric(log_B2),
      BF_trait1 = exp(pmin(as.numeric(log_B1), 700)),
      BF_trait2 = exp(pmin(as.numeric(log_B2), 700)),
      has_data_trait1 = has1,
      has_data_trait2 = has2,
      n_variants_trait1 = nvar1,
      n_variants_trait2 = nvar2,
      stringsAsFactors = FALSE
    ),
    parameters = list(
      pi = pi,
      delta_trait1 = unname(pi[["10"]] + pi[["11"]]),
      delta_trait2 = unname(pi[["01"]] + pi[["11"]]),
      eta_trait1 = eta1,
      eta_trait2 = eta2
    ),
    bfdr = list(
      trait1 = bfdr_table(genes, pp_trait1, trait_names[1]),
      trait2 = bfdr_table(genes, pp_trait2, trait_names[2]),
      pleiotropy = bfdr_table(genes, pp_pleio, "pleiotropy")
    ),
    traces = list(
      loglik = trace_loglik,
      max_parameter_change = trace_change,
      pi = do.call(rbind, trace_pi),
      eta_trait1 = do.call(rbind, trace_eta1),
      eta_trait2 = do.call(rbind, trace_eta2)
    ),
    variant_info = list(
      trait1 = transform(d1, rho = compute_rho(d1, eta1)),
      trait2 = transform(d2, rho = compute_rho(d2, eta2))
    ),
    diagnostics = list(
      converged = converged,
      n_iter = length(trace_loglik),
      loglik_non_decreasing = monotone,
      eta_trait1 = eta_diag1,
      eta_trait2 = eta_diag2,
      boundary_pi = pi < 1e-8,
      missing_gene_data = data.frame(
        Gene = genes,
        has_data_trait1 = has1,
        has_data_trait2 = has2,
        stringsAsFactors = FALSE
      )
    )
  )
  class(out) <- "mirage_two_trait"
  out
}

print.mirage_two_trait <- function(x, ...) {
  cat("Two-trait MIRAGE fit\n")
  cat("Traits:", paste(x$trait_names, collapse = ", "), "\n")
  cat("Genes:", nrow(x$posterior), "\n")
  cat("Iterations:", x$diagnostics$n_iter, "\n")
  cat("Converged:", x$diagnostics$converged, "\n")
  cat("Pi:", paste(names(x$parameters$pi), signif(x$parameters$pi, 4), sep = "=", collapse = ", "), "\n")
  invisible(x)
}

