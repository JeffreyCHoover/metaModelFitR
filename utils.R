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
