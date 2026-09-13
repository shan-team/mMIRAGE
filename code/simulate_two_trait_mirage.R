simulate_two_trait_mirage <- function(
    n_genes = 200,
    variants_per_gene = 8,
    pi = c("00" = 0.85, "10" = 0.05, "01" = 0.05, "11" = 0.05),
    eta_trait1 = c("1" = 0.20, "2" = 0.05),
    eta_trait2 = c("1" = 0.20, "2" = 0.05),
    signal_mean_trait1 = 2.5,
    signal_mean_trait2 = 2.5,
    null_mean = -0.7,
    null_sd = 0.15,
    signal_sd = 0.30,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pi <- normalize_pi(pi)
  genes <- paste0("GENE", seq_len(n_genes))
  states <- sample(names(pi), n_genes, replace = TRUE, prob = pi)
  make_trait <- function(trait, eta, signal_mean, active) {
    rows <- vector("list", n_genes)
    cats <- names(eta)
    for (i in seq_len(n_genes)) {
      category <- sample(cats, variants_per_gene, replace = TRUE)
      z <- if (active[i]) {
        rbinom(variants_per_gene, 1, eta[category])
      } else {
        rep(0, variants_per_gene)
      }
      log_bf <- rnorm(variants_per_gene, null_mean, null_sd)
      log_bf[z == 1] <- rnorm(sum(z == 1), signal_mean, signal_sd)
      rows[[i]] <- data.frame(
        ID = paste0("t", trait, "_", genes[i], "_v", seq_len(variants_per_gene)),
        Gene = genes[i],
        No.case = 0,
        No.contr = 0,
        category = category,
        log_bf = log_bf,
        true_z = z,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }
  active1 <- states %in% c("10", "11")
  active2 <- states %in% c("01", "11")
  list(
    trait1 = make_trait(1, eta_trait1, signal_mean_trait1, active1),
    trait2 = make_trait(2, eta_trait2, signal_mean_trait2, active2),
    truth = data.frame(
      Gene = genes,
      state = states,
      U1 = as.integer(active1),
      U2 = as.integer(active2),
      pleiotropic = as.integer(states == "11"),
      stringsAsFactors = FALSE
    ),
    parameters = list(pi = pi, eta_trait1 = eta_trait1, eta_trait2 = eta_trait2)
  )
}
