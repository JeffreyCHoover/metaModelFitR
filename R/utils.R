#' @title Retrieve meta-analysis result components
#'
#' @description This function returns the observed effect sizes, weights, and
#' number of included studies from an `rma` object.
#'
#' @param meta An `rma` object from the metafor package.
#' @param fixed A boolean variable indicating whether the meta-analysis was
#' using a fixed-effect or random-effects model
#'
#' @return results A data frame containing vectors for the study effect sizes,
#' weights, and the number of included studies.
meta_retrieve <- function(meta, fixed = FALSE)
{
  if(is.null(meta)) {
    stop("Meta-analysis results parameter is required.")
  } else if (!("rma" %in% class(meta))) {
    stop("Meta-analysis results must be an object of class `rma`.")
  }

  tau <- dplyr::case_when(fixed == TRUE ~ 0,
                          fixed == FALSE ~ meta$tau2)

  results <- tibble::as_tibble(data.frame(meta$yi, meta$vi)) %>%
    dplyr::rename(effect = meta.yi, var = meta.vi) %>%
    dplyr::mutate(tau = tau,
                  weights = 1 / (var + tau),
                  k = nrow(.)) %>%
    dplyr::select(effect, weights, k)

  return(results)
}

#' @title Simulate meta-analysis results based on an actual meta-analysis
#'
#' @description This function simulates results for a specified number of
#' meta-analyses based on the results from the actual meta-analysis.
#'
#' @param effect The pooled effect size from the meta-analysis.
#' @param effect_sd The standard error of the pooled effect size.
#' @param weight The mean value of the vector of weights for the studies
#' included in the meta-analysis
#' @param weight_sd The standard deviation of the vector of weights for the
#' studies included in the meta-analysis.
#' @param studies The number of included studies in the meta-analysis.
#' @param iter The number of meta-analyses results to simulate.
#'
#' @return results A data frame containing the value of the discrepancy function
#' for each of the simulated meta-analyses.
simulate <- function(effect, effect_sd, weight, weight_sd, studies,
                     iter = 100000)
{
  results <- matrix(nrow = iter, ncol = 6, dimnames = list(c(1:iter),
                                                           c("Smallest Effect",
                                                             "Largest Effect",
                                                             "Mean Weight",
                                                             "Smallest Weight",
                                                             "Largest Weight",
                                                             "Mean Effect")))

  i <- 1
  while(i <= iter)
  {
    j <- 1
    step <- matrix(nrow = studies, ncol = 2, dimnames = list(c(1:studies),
                                                             c("Effect",
                                                               "Weight")))

    while(j <= studies)
    {
      step[j, 1] <- stats::rnorm(1, mean = effect, sd = effect_sd)
      step[j, 2] <- stats::rnorm(1, mean = weight, sd = weight_sd)

      j <- j + 1
    }

    results[i, 1] <- min(step[,1])
    results[i, 2] <- max(step[,1])
    results[i, 3] <- mean(step[,2])
    results[i, 4] <- min(step[,2])
    results[i, 5] <- max(step[,2])
    results[i, 6] <- sum(step[,1] * step[,2]) / sum(step[,2])

    i <- i + 1
  }

  return(as.data.frame(results))
}

#' @title Plot the Monte Carlo Resampling results
#'
#' @description This function plots the results of the meta-analysis to the
#' results of the simualted meta-analyses on the seven specified discrepancy
#' measures.
#'
#' @param simulated The data frame of values from the simulated meta-analyses
#' based on the discrepancy measures.
#' @param observed_effects A vector of observed effect sizes from the
#' meta-analysis.
#' @param observed_weights A vector of observed weights from the
#' meta-analysis.
#' @param meta_es The pooled effect size from the meta-analysis.
#' @param filename A character variable that serves as the stem for all saved
#'`   output files.
plot_simulation <- function(simulated, observed_effects, observed_weights,
                            meta_es, filename)
{
  mean_effect_obs <- meta_es
  min_effect_obs <- min(observed_effects)
  max_effect_obs <- max(observed_effects)
  mean_weight_obs <- mean(observed_weights)
  min_weight_obs <- min(observed_weights)
  max_weight_obs <- max(observed_weights)

  if(!dir.exists(here::here("/figures"))) {
    dir.create(here::here("/figures"))
  }

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-minimum_effect.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")

  print(simulated %>%
          ggplot2::ggplot(ggplot2::aes(x = `Smallest Effect`)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00", bins = 30) +
          ggplot2::geom_vline(xintercept = min_effect_obs) +
          ggplot2::labs(y = "Count"))

  grDevices::dev.off()

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-maximum_effect.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")

  print(simulated %>%
          ggplot2::ggplot(ggplot2::aes(x = `Largest Effect`)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00", bins = 30) +
          ggplot2::geom_vline(xintercept = max_effect_obs) +
          ggplot2::labs(y = "Count"))

  grDevices::dev.off()

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-mean_weight.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")

  print(simulated %>%
          ggplot2::ggplot(ggplot2::aes(x = `Mean Weight`)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00", bins = 30) +
          ggplot2::geom_vline(xintercept = mean_weight_obs) +
          ggplot2::labs(y = "Count"))

  grDevices::dev.off()

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-minimum_weight.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")

  print(simulated %>%
          ggplot2::ggplot(ggplot2::aes(x = `Smallest Weight`)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00", bins = 30) +
          ggplot2::geom_vline(xintercept = min_weight_obs) +
          ggplot2::labs(y = "Count"))

  grDevices::dev.off()

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-maximum_weight.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")

  print(simulated %>%
          ggplot2::ggplot(ggplot2::aes(x = `Largest Weight`)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00", bins = 30) +
          ggplot2::geom_vline(xintercept = max_weight_obs) +
          ggplot2::labs(y = "Count"))

  grDevices::dev.off()

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-mean_effect.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")

  print(simulated %>%
          ggplot2::ggplot(ggplot2::aes(x = `Mean Effect`)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00", bins = 45) +
          ggplot2::geom_vline(xintercept = mean_effect_obs) +
          ggplot2::xlim(min_effect_obs, max_effect_obs) +
          ggplot2::labs(y = "Count", x = "Pooled Effect Size"))

  grDevices::dev.off()
}

#' @title Calculate Monte Carlo p-values for meta-analysis results
#'
#' @description This function calculates the Monte Carlo p-values for the
#' meta-analysis results
#'
#' @param simulated The data frame of values from the simulated meta-analyses
#' based on the discrepancy measures.
#' @param observed_effects A vector of observed effect sizes from the
#' meta-analysis.
#' @param observed_weights A vector of observed weights from the
#' meta-analysis.
#' @param meta_es The pooled effect size from the meta-analysis.
#' @return simulated a vector containing the Monte Carlo p-values for each of
#' the discrepancy measures.
ma_monte_p_value <- function(simulated, observed_effects, observed_weights,
                             meta_es)
{
  simulated <- simulated %>%
    dplyr::mutate(mean_effect = meta_es,
                  smallest_effect = min(observed_effects),
                  largest_effect = max(observed_effects),
                  mean_weight = mean(observed_weights),
                  smallest_weight = min(observed_weights),
                  largest_weight = max(observed_weights),
                  m_e_monte_p = dplyr::case_when(mean_effect > `Mean Effect`~ 1,
                                                 TRUE ~ 0),
                  s_e_monte_p = dplyr::case_when(smallest_effect >
                                                   `Smallest Effect` ~ 1,
                                                 TRUE ~ 0),
                  l_e_monte_p = dplyr::case_when(largest_effect >
                                                   `Largest Effect` ~ 1,
                                                 TRUE ~ 0),
                  w_monte_p = dplyr::case_when(mean_weight > `Mean Weight` ~ 1,
                                               TRUE ~ 0),
                  s_w_monte_p = dplyr::case_when(smallest_weight >
                                                   `Smallest Weight` ~ 1,
                                                 TRUE ~ 0),
                  l_w_monte_p = dplyr::case_when(largest_weight >
                                                   `Largest Weight` ~ 1,
                                                 TRUE ~ 0)) %>%
    dplyr::summarize(mean_es_monte_p = mean(m_e_monte_p),
                     smallest_effect_monte_p = mean(s_e_monte_p),
                     largest_effect_monte_p = mean(l_e_monte_p),
                     mean_weight_monte_p = mean(w_monte_p),
                     smallest_weight_monte_p = mean(s_w_monte_p),
                     largest_weight_monte_p = mean(l_w_monte_p)) %>%
    dplyr::select(smallest_effect_monte_p, largest_effect_monte_p,
                  mean_weight_monte_p, smallest_weight_monte_p,
                  largest_weight_monte_p, mean_es_monte_p)

  return(simulated)
}

#' @title Calculate Posterior Predictive p-values for meta-analysis results
#'
#' @description This function calculates the Posterior Predictive p-values for the
#' meta-analysis results
#'
#' @param mcmc The data frame containing the Markov Chain Monte Carlo values
#' @param observed_effects A vector of observed effect sizes from the
#' meta-analysis.
#' @param observed_weights A vector of observed weights from the
#' meta-analysis.
#' @return mcmc A vector containing the Postrior Predictive p-values for each of
#' the discrepancy measures.
calc_ppp <- function(mcmc, observed_effects, observed_weights)
{
  obtained_es <- sum(observed_effects * observed_weights) /sum(observed_weights)

  mcmc <- mcmc %>%
    dplyr::mutate(smallest_effect = min(observed_effects),
                  largest_effect = max(observed_effects),
                  mean_weight = mean(observed_weights),
                  smallest_weight = min(observed_weights),
                  largest_weight = max(observed_weights),
                  min_effect_ppp = dplyr::case_when(minEffect <
                                                      smallest_effect ~ 1,
                                                    TRUE ~ 0),
                  max_effect_ppp = dplyr::case_when(maxEffect <
                                                      largest_effect ~ 1,
                                                    TRUE ~ 0),
                  min_weight_ppp = dplyr::case_when(minWeight <
                                                      smallest_weight ~ 1,
                                                    TRUE ~ 0),
                  mean_weight_ppp = dplyr::case_when(meanWeight <
                                                       mean_weight ~ 1,
                                                     TRUE ~ 0),
                  max_weight_ppp = dplyr::case_when(maxWeight <
                                                      largest_weight ~ 1,
                                                    TRUE ~ 0),
                  es_ppp = dplyr::case_when(ES_agg < obtained_es ~ 1,
                                            TRUE ~ 0)) %>%
    dplyr::summarize(min_effect_ppp_value = mean(min_effect_ppp),
                     max_effect_ppp_value = mean(max_effect_ppp),
                     min_weight_ppp_value = mean(min_weight_ppp),
                     mean_weight_ppp_value = mean(mean_weight_ppp),
                     max_weight_ppp_value = mean(max_weight_ppp),
                     es_ppp_value = mean(es_ppp)) %>%
    dplyr::select(min_effect_ppp_value,
                  max_effect_ppp_value, min_weight_ppp_value,
                  mean_weight_ppp_value, max_weight_ppp_value, es_ppp_value)

  return(mcmc)
}

#' @title Calculate Posterior Predictive p-values for meta-analysis results
#'
#' @description Graph the Posterior Predictive p-values for the meta-analysis
#' results
#'
#' @param mcmc The data frame containing the Markov Chain Monte Carlo values
#' @param observed_effects A vector of observed effect sizes from the
#' meta-analysis.
#' @param observed_weights A vector of observed weights from the
#' meta-analysis.
#' @param meta An 'rma' object that contains the results of a meta-analysis
#'    performed in the metafor package.
#' @param filename A character variable that serves as the stem for all saved
#'`   output files.
plot_ppmc <- function(mcmc, observed_effects, observed_weights, meta, filename)
{
  if(!dir.exists(here::here("/figures"))) {
    dir.create(here::here("/figures"))
    }

  min_effect_obs <- min(observed_effects)
  max_effect_obs <- max(observed_effects)
  mean_weight_obs <- mean(observed_weights)
  min_weight_obs <- min(observed_weights)
  max_weight_obs <- max(observed_weights)

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-minimum_effect.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")
  print(mcmc %>%
          ggplot2::ggplot(ggplot2::aes(x = minEffect)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00") +
          ggplot2::geom_vline(xintercept = min_effect_obs))
  grDevices::dev.off()

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-maximum_effect.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")
  print(mcmc %>%
          ggplot2::ggplot(ggplot2::aes(x = maxEffect)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00") +
          ggplot2::geom_vline(xintercept = max_effect_obs))
  grDevices::dev.off()

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-mean_weight.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")
  print(mcmc %>%
          ggplot2::ggplot(ggplot2::aes(x = meanWeight)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00") +
          ggplot2::geom_vline(xintercept = mean_weight_obs))
  grDevices::dev.off()

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-minimum_weight.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")
  print(mcmc %>%
          ggplot2::ggplot(ggplot2::aes(x = minWeight)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00") +
          ggplot2::geom_vline(xintercept = min_weight_obs))
  grDevices::dev.off()

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-maximum_weight.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")
  print(mcmc %>%
          ggplot2::ggplot(ggplot2::aes(x = maxWeight)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00") +
          ggplot2::geom_vline(xintercept = max_weight_obs))
  grDevices::dev.off()

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-pooled_effect.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")
  print(mcmc %>%
          ggplot2::ggplot(ggplot2::aes(x = ES_agg)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00") +
          ggplot2::geom_vline(xintercept = meta$b))
  grDevices::dev.off()
}

#' @title Perform diagnostics on the Markov Chain Monte Carlo chains
#'
#' @description The function tests for autocorrelation, parameter values, and
#' density. The results of these tests are saved as jpeg files for easy
#' reference and insertion into manuscripts.
#'
#' @param codaObject The Markov Chain Monte Carlo chains.
#' @param parName A character for the name of the discrepancy measure that is
#' being tested.
#' @param filename A character for the filename root, which is used when saving graphs.
diagMCMC <- function(codaObject, parName = varnames(codaObject)[1], filename) {
  if(!dir.exists(here("/diagnostics"))) {
    dir.create(here("/diagnostics"))
    }

  grDevices::jpeg(filename =
                    here::here(glue::glue("diagnostics/{filename}-{parName}-param.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")
  print(codaObject[[1]] %>%
          tibble::as_tibble() %>%
          dplyr::select(parName) %>%
          dplyr::mutate(chain = 1,
                 row_num = row_number()) %>%
          rbind(codaObject[[2]] %>%
                  as_tibble() %>%
                  select(parName) %>%
                  mutate(chain = 2,
                         row_num = row_number())) %>%
          rbind(codaObject[[3]] %>%
                  as_tibble() %>%
                  select(parName) %>%
                  mutate(chain = 3,
                         row_num = row_number())) %>%
          dplyr::mutate(chain = factor(chain, levels = c(1, 2, 3))) %>%
          ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes_string(x = "row_num", y = parName,
                                        color = "chain", group = "chain"),
                    stat = "identity", show.legend = FALSE) +
          ggplot2::labs(x = "Iteration", y = "Parameter Value"))
  grDevices::dev.off()

  x_matrix = NULL
  y_matrix = NULL

  for(i in 1:length(codaObject)) {
    acf_output <- acf(codaObject[,c(parName)][[i]], plot = FALSE)
    x_matrix <- cbind(x_matrix, acf_output$lag)
    y_matrix <- cbind(y_matrix, acf_output$acf)
  }

  grDevices::jpeg(filename =
                    here::here(glue::glue("diagnostics/{filename}-{parName}-autocorr.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")
  print(x_matrix %>%
          cbind(y_matrix) %>%
          tibble::as_tibble() %>%
          dplyr::mutate(ess =
                          round(coda::effectiveSize(codaObject[,c(parName)]),
                                digits = 2)) %>%
          ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(x = V1, y = V4)) +
          ggplot2::geom_line(ggplot2::aes(x = V2, y = V5)) +
          ggplot2::geom_line(ggplot2::aes(x = V3, y = V6)) +
          ggplot2::geom_hline(yintercept = 0) +
          ggplot2::geom_text(ggplot2::aes(label = glue::glue("ESS: {ess}")), x = Inf,
                             y = Inf, hjust = 1, vjust = 1) +
          ggplot2::labs(x = "Lag", y = "Autocorrelation") +
          ggplot2::scale_y_continuous(breaks = c(0, .1, .2, .3, .4, .5, .6, .7,
                                                 .8, .9, 1)))
  grDevices::dev.off()

  x_matrix = NULL
  y_matrix = NULL

  for(i in 1:length(codaObject)) {
    density_info <- density(codaObject[,c(parName)][[i]])
    x_matrix <- cbind(x_matrix, density_info$x)
    y_matrix <- cbind(y_matrix, density_info$y)
  }

  grDevices::jpeg(filename =
                    here::here(glue::glue("diagnostics/{filename}-{parName}-density.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300,
                  type = "cairo")
  print(x_matrix %>%
          cbind(y_matrix) %>%
          tibble::as_tibble() %>%
          dplyr::mutate(ess =
                          round(coda::effectiveSize(codaObject[,c(parName)]),
                                digits = 2),
                    mcse = round(stats::sd(as.matrix(codaObject[,c(parName)])) /
                                       sqrt(ess), digits = 5)) %>%
          ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(x = V1, y = V4)) +
          ggplot2::geom_line(ggplot2::aes(x = V2, y = V5)) +
          ggplot2::geom_line(ggplot2::aes(x = V3, y = V6)) +
          ggplot2::geom_hline(yintercept = 0) +
          ggplot2::geom_text(ggplot2::aes(label = glue::glue("ESS: {ess}")), x = -Inf,
                             y = Inf, hjust = 0, vjust = 1) +
          ggplot2::geom_text(ggplot2::aes(label = glue::glue("MCSE: {mcse}")), x = Inf,
                             y = Inf, hjust = 1, vjust = 1) +
          ggplot2::labs(x = "Parameter Value", y = "Density") +
          ggplot2::scale_y_continuous(breaks = c(0, .1, .2, .3, .4, .5, .6, .7,
                                                 .8, .9, 1)))
  grDevices::dev.off()
}
