context("test-utility_functions")

test_that("meta_retrieve works", {
  library(metafor)
  es <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  v <- c(.05, .1, .08, .1, .03)
  w <- 1 / v
  d <- tibble::tibble(es, v)
  re_meta <- rma(data = d, yi = es, vi = v, method = "FE")

  test <- meta_retrieve(re_meta, fixed = TRUE)

  expect_equal(test$weights, w)
  expect_equal(as.vector(test$effect), es)
  expect_equal(test$k, rep(length(es), length(es)))
})
