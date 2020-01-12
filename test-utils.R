context ("test-utility_functions")

test_that ("meta_retrieve works", {
  es <- c (0.5, 0.6, 0.7, 0.8, 0.9)
  v <- c (.05, .1, .08, .1, .03)
  w <- 1 / v
  d <- tibble::tibble (es, v)
  re_meta <- metafor::rma (data = d, yi = es, vi = v, method = "FE")

  test <- meta_retrieve (re_meta, fixed = TRUE)

  testthat::expect_equal(test$weights, w)
  testthat::expect_equal(as.vector(test$effect), es)
  testthat::expect_equal(test$k, rep(length(es), length(es)))
})

test_that ("calc_ppp works", {
  discrepancy <- 4
  simulated_discrepancy <- c(2, 3, 5, 6, 7)

  testthat::expect_equal (.4, calc_ppp(c(4),
                            simulated_discrepancy = simulated_discrepancy))
})
