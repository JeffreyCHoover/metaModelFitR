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
  meta.yi <- NULL
  meta.vi <- NULL
  var <- NULL
  effect <- NULL
  weights <- NULL
  k <- NULL
  . <- NULL

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
  results <- matrix(nrow = iter, ncol = 1, dimnames = list(c(1:iter),
                                                           c("Mean Effect")))

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
      #step[j, 2] <- stats::rnorm(1, mean = weight, sd = weight_sd)
      step[j, 2] <- stats::runif(1, min = min(weight), max(weight))

      j <- j + 1
    }

    results[i, 1] <- sum(step[,1] * step[,2]) / sum(step[,2])

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
  `Mean Effect` <- NULL

  mean_effect_obs <- meta_es

  if(!dir.exists(here::here("/figures"))) {
    dir.create(here::here("/figures"))
  }

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
  m_e_monte_p <- NULL
  mean_es_monte_p <- NULL

  simulated <- simulated %>%
    dplyr::mutate(mean_effect = meta_es,
                  m_e_monte_p = dplyr::case_when(mean_effect > `Mean Effect`~ 1,
                                                 TRUE ~ 0)) %>%
    dplyr::summarize(mean_es_monte_p = mean(m_e_monte_p)) %>%
    dplyr::select(mean_es_monte_p) %>%
    dplyr::rename(es_monte_p = mean_es_monte_p)

  return(simulated)
}

#' @title Calculate Posterior Predictive p-values for meta-analysis results
#'
#' @description This function calculates the Posterior Predictive p-values for the
#' meta-analysis results
#'
#' @param mcmc The data frame containing the Markov Chain Monte Carlo values
#' @param discrepancy A vector of observed discrepancies between the obtained
#' study-level effects.
#' @param simulated_discrepancy A vector of the simulated discrepancies of the
#' simulated study-level effects.
#' @return mcmc A vector containing the Postrior Predictive p-values for each of
#' the discrepancy measures.
calc_ppp <- function(mcmc, discrepancy, simulated_discrepancy)# observed_effects, observed_weights)
{
  #obtained_es <- sum(observed_effects * observed_weights) /sum(observed_weights)
  es_ppp <- NULL

  temp <- simulated_discrepancy %>%
    dplyr::mutate(es_ppp = dplyr::case_when(sim_disc < discrepancy[1,1] ~ 1,
                                            TRUE ~ 0)) %>%
    dplyr::summarize(es_ppp_value = mean(es_ppp)) %>%
    dplyr::pull()

  return(round(temp, 4))
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
#' @param discrepancy A vector of observed discrepancies between the obtained
#' study-level effects.
#' @param simulated_discrepancy A vector of the simulated discrepancies of the
#' simulated study-level effects.
plot_ppmc <- function(mcmc, observed_effects, observed_weights, meta, filename,
                      discrepancy, simulated_discrepancy)
{
  sim_disc <- NULL

  if(!dir.exists(here::here("/figures"))) {
    dir.create(here::here("/figures"))
    }

  grDevices::jpeg(filename =
                    here::here(glue::glue("figures/{filename}-pooled_effect.jpeg")),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 100,
                  type = "cairo")
  print(simulated_discrepancy %>%
          ggplot2::ggplot(ggplot2::aes(x = sim_disc)) +
          ggplot2::geom_histogram(color = "black", fill = "#D55E00",
                                  binwidth = 50) +
          ggplot2::geom_vline(xintercept = discrepancy[1,1], size = 1.5))

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
diagMCMC <- function(codaObject, parName, filename) {
  chain <- NULL
  V1 <- NULL
  V2 <- NULL
  V3 <- NULL
  V4 <- NULL
  V5 <- NULL
  V6 <- NULL
  ess <- NULL

  if(!dir.exists(here::here("/diagnostics"))) {
    dir.create(here::here("/diagnostics"))
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
                 row_num = dplyr::row_number()) %>%
          rbind(codaObject[[2]] %>%
                  tibble::as_tibble() %>%
                  dplyr::select(parName) %>%
                  dplyr::mutate(chain = 2,
                         row_num = dplyr::row_number())) %>%
          rbind(codaObject[[3]] %>%
                  tibble::as_tibble() %>%
                  dplyr::select(parName) %>%
                  dplyr::mutate(chain = 3,
                         row_num = dplyr::row_number())) %>%
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
    acf_output <- stats::acf(codaObject[,c(parName)][[i]], plot = FALSE)
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
    density_info <- stats::density(codaObject[,c(parName)][[i]])
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
