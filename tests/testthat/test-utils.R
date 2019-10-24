context("test-utility_functions")

test_that("meta_retrieve works", {
  es <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  v <- c(.05, .1, .08, .1, .03)
  w <- 1 / v
  d <- tibble::tibble(es, v)
  re_meta <- metafor::rma(data = d, yi = es, vi = v, method = "FE")

  test <- meta_retrieve(re_meta, fixed = TRUE)

  testthat::expect_equal(test$weights, w)
  testthat::expect_equal(as.vector(test$effect), es)
  testthat::expect_equal(test$k, rep(length(es), length(es)))
})

test_that("calc_ppp works", {
  discrepancy <- matrix(data = c(4), nrow = 1, ncol = 1)
  simulated_discrepancy <- data.frame(c(2, 3, 5, 6, 7)) %>%
    dplyr::rename(sim_disc = 1)

  testthat::expect_equal(.4, calc_ppp(discrepancy = discrepancy,
                            simulated_discrepancy = simulated_discrepancy))
})
