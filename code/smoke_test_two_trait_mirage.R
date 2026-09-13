source("code/two_trait_mirage.R")
source("code/simulate_two_trait_mirage.R")

sim <- simulate_two_trait_mirage(
  n_genes = 60,
  variants_per_gene = 4,
  pi = c("00" = 0.70, "10" = 0.10, "01" = 0.10, "11" = 0.10),
  seed = 1
)

fit <- mirage_two_trait(
  sim$trait1,
  sim$trait2,
  n1_trait1 = 10000,
  n0_trait1 = 10000,
  n1_trait2 = 10000,
  n0_trait2 = 10000,
  log_bf_col_trait1 = "log_bf",
  log_bf_col_trait2 = "log_bf",
  max_iter = 500
)

stopifnot(inherits(fit, "mirage_two_trait"))
stopifnot(nrow(fit$posterior) == 60)
stopifnot(all(names(fit$parameters$pi) == c("00", "10", "01", "11")))
stopifnot(abs(sum(fit$parameters$pi) - 1) < 1e-8)
stopifnot(all(fit$posterior$PP_trait1 >= 0 & fit$posterior$PP_trait1 <= 1))
stopifnot(all(fit$posterior$PP_trait2 >= 0 & fit$posterior$PP_trait2 <= 1))
stopifnot(all(fit$posterior$PP_pleiotropy >= 0 & fit$posterior$PP_pleiotropy <= 1))

print(fit)
cat("Smoke test passed.\n")
