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

test_that("ma_monte_p_value works", {
  meta_es <- 0.6
  obs_es <- c(0.4, 0.5, 0.6, 0.7, 0.8)
  obs_w <- c(10, 12, 14, 12, 10)
  sim_med_es <- rep(c(0.4, 0.5, 0.6, 0.7, 0.8), 2000)
  sim_min_es <- rep(c(0.2, 0.3, 0.4, 0.5, 0.6), 2000)
  sim_max_es <- rep(c(.6, .7, .8, .9, 1), 2000)
  sim_mean_w <- rep(c(10, 12, 14, 12, 10), 2000)
  sim_min_w <- rep(c(6, 8, 10, 8, 6), 2000)
  sim_max_w <- rep(c(10, 12, 14, 12, 10), 2000)
  sim_es <- rep(c(.5, .55, .6, .55, .5), 2000)

  smallest_effect_monte_p <- .4
  largest_effect_monte_p <- .4
  mean_weight_monte_p <- .4
  smallest_weight_monte_p <- .8
  largest_weight_monte_p <- .8
  mean_es_monte_p <- .8

  test_results <- tibble::tibble(smallest_effect_monte_p, largest_effect_monte_p,
                                 mean_weight_monte_p, smallest_weight_monte_p,
                                 largest_weight_monte_p, mean_es_monte_p)

  sim <- tibble::tibble(sim_med_es, sim_min_es, sim_max_es, sim_mean_w,
                        sim_min_w, sim_max_w, sim_es) %>%
    dplyr::rename(`Smallest Effect` = sim_min_es,
                  `Largest Effect` = sim_max_es,
                  `Mean Weight` = sim_mean_w,
                  `Smallest Weight` = sim_min_w,
                  `Largest Weight` = sim_max_w,
                  `Mean Effect` = sim_es)
  expect_equal(ma_monte_p_value(sim, obs_es, obs_w, meta_es), test_results)
})

test_that("calc_ppp works", {
  meta_es <- .6
  obs_es <- c(0.4, 0.5, 0.6, 0.7, 0.8)
  obs_w <- c(10, 12, 14, 12, 10)

  sim_mean_es <- rep(c(0.4, 0.5, 0.6, 0.7, 0.8), 2000)
  sim_min_es <- rep(c(0.2, 0.3, 0.4, 0.5, 0.6), 2000)
  sim_max_es <- rep(c(.6, .7, .8, .9, 1), 2000)
  sim_mean_w <- rep(c(10, 12, 14, 12, 10), 2000)
  sim_min_w <- rep(c(6, 8, 10, 8, 6), 2000)
  sim_max_w <- rep(c(10, 12, 14, 12, 10), 2000)
  sim_es <- rep(c(.5, .55, .6, .55, .5), 2000)

  min_effect_ppp_value <- .4
  mean_effect_ppp_value <- .4
  max_effect_ppp_value <- .4
  mean_weight_ppp_value <- .4
  min_weight_ppp_value <- .8
  max_weight_ppp_value <- .8
  es_ppp_value <- .8

  test_results <- tibble::tibble(min_effect_ppp_value, mean_effect_ppp_value,
                                 max_effect_ppp_value,
                                 mean_weight_ppp_value, min_weight_ppp_value,
                                 max_weight_ppp_value, es_ppp_value)

  sim <- tibble::tibble(sim_mean_es, sim_min_es, sim_max_es, sim_mean_w,
                        sim_min_w, sim_max_w, sim_es) %>%
    dplyr::rename(`minEffect` = sim_min_es,
                  `meanEffect` = sim_mean_es,
                  `maxEffect` = sim_max_es,
                  `meanWeight` = sim_mean_w,
                  `minWeight` = sim_min_w,
                  `maxWeight` = sim_max_w,
                  `ES_agg` = sim_es)
  expect_equal(calc_ppp(sim, obs_es, obs_w), test_results)
})
