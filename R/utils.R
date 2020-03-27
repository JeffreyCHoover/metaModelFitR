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
meta_retrieve <- function(meta, fixed = FALSE) {
  meta.yi <- NULL
  meta.vi <- NULL
  var <- NULL
  effect <- NULL
  weights <- NULL
  k <- NULL
  . <- NULL

  if (is.null(meta)) {
    stop ("Meta-analysis results parameter is required.")
  } else if (!("rma" %in% class(meta))) {
    stop ("Meta-analysis results must be an object of class `rma`.")
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

#' @title Calculate Posterior Predictive p-values for meta-analysis results
#'
#' @description This function calculates the Posterior Predictive p-values for the
#' meta-analysis results
#'
#' @param discrepancy A numeric value of observed discrepancies between the
#' obtained study-level effects.
#' @param simulated_discrepancy A vector of the simulated discrepancies of the
#' simulated study-level effects.
#'
#' @return mcmc A vector containing the Postrior Predictive p-values for each of
#' the discrepancy measures.
calc_ppp <- function(discrepancy, simulated_discrepancy) {
  es_ppp <- NULL
  name <- NULL

  temp <- simulated_discrepancy %>%
    tibble::enframe() %>%
    dplyr::rename(sim_disc = 2) %>%
    dplyr::select(-name) %>%
    dplyr::mutate(disc = discrepancy,
                  es_ppp = dplyr::case_when(.$sim_disc < disc ~ 1,
                                            TRUE ~ 0)) %>%
    dplyr::summarize(es_ppp_value = mean(es_ppp)) %>%
    dplyr::pull()

  return(round(temp, 4))
}

#' @title Retrieve standard errors from confidence intervals
#'
#' @description This function returns the standard error from a reported
#' confidence interval.
#'
#' @param estimate A `numeric` object for the point estimate of the confidence
#' interval.
#' @param lower A `numeric` object for the lower bound of a confidence interval.
#' @param upper A `numeric` object for the upper bound of a confidence interval.
#' @param ci A `numeric` object for the percentage of the confidence interval.
#' @param corr A `boolean` object for whether the confidence interval is for a
#' correlation coefficient.
#'
#' @return se A numeric object containing the standard error of the reported
#' confidence interval.
se_retrieve <- function(estimate, lower = NA_real_, upper = NA_real_, ci = .95,
                        corr = FALSE) {
  if(is.na(lower) & is.na(upper)) {
    warning("Either the lower or upper bound of the confidence interval must be entered.")
  } else if(!is.numeric(lower) & !is.numeric(upper)) {
    warning("Enter a numeric value for either the lower or upper bound of the confidence interval.")
  } else if(!is.numeric(estimate)) {
    warning("Enter a numeric value for the point estimate")
  } else if(!is.numeric(ci)) {
    warning("Enter a numeric value for the confidence interval percentage.")
  } else if((ci <= 0) | (ci >= 1)){
    warning("Enter a value between 0 and 1 for the confidence interval percentage.")
  }

  mult <- qnorm(ci + ((1 - ci) / 2))

  if(!corr) {
    if(!is.na(upper)) {
      se = (upper - estimate) / mult
    } else {
      se = (estimate - lower) / mult
    }
    return(se)
  } else {
    z_est = .5 * log((1 + estimate) / (1 - estimate))
    if(!is.na(upper)) {
      z_upper = .5 * log((1 + upper) / (1 - upper))
      z_se = (z_upper - z_est) / mult
    } else {
      z_lower = .5 * log((1 + lower) / (1 - lower))
      z_se = (z_est - z_lower) / mult
    }
    return(z_se)
  }
}
