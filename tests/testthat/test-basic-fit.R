test_that("package simulation and joint fit run", {
  sim <- simulate_case_control_two_trait(
    n_genes = 80,
    variants_per_gene = 6,
    seed = 11
  )
  fit <- mirage_two_trait(
    sim$trait1,
    sim$trait2,
    n1_trait1 = 3000,
    n0_trait1 = 3000,
    n1_trait2 = 3000,
    n0_trait2 = 3000,
    log_bf_col_trait1 = "log_bf",
    log_bf_col_trait2 = "log_bf",
    pi_init = sim$parameters$pi,
    eta_init_trait1 = sim$parameters$eta_trait1,
    eta_init_trait2 = sim$parameters$eta_trait2,
    max_iter = 50
  )
  expect_s3_class(fit, "mirage_two_trait")
  expect_equal(nrow(fit$posterior), 80)
  expect_equal(sum(fit$parameters$pi), 1, tolerance = 1e-8)
})
