#' @title Monte Carlo Resampling for Meta-Analysis Models
#'
#' @description This function performs Monte Carlo resampling on meta-analysis
#' results from the metafor package to assess model fit. This function returns
#' Monte Carlo p-values for the following discrepancy measures: mean study-level
#' effects, minimum study-level effects, maximum study-level effects, mean
#' study-level weights, minimum study-level weights, maximum study-level effects
#' and pooled effect size. Graphical results are also saved for each of these
#' discrepancy functions.
#'
#' @param meta An 'rma' object that contains the results of a meta-analysis
#'   performed in the metafor package.
#' @param fixed Logical variable specifying whether the meta-analysis results
#'   were using a fixed-effect model (default) or a random-effects model.
#' @param fileName A character variable that serves as the stem for all saved
#'`   output files.
#'
#' @return monte-p_values A vector containing the Monte Carlo p-values for each
#'   of the discrepancy measures.
#'
#' @examples
#' \dontrun{
#' monte_carlo_ma(meta = re_meta, fixed = FALSE, fileName = "re-meta_monte-carlo")
#' }
#'
#' @export

monte_carlo_ma <- function(meta, fixed = FALSE, fileName)
{
  options(warn = -1)

  if(is.null(meta)) {
    stop("Meta-analysis results parameter is required.")
  } else if (!("rma" %in% class(meta))) {
      stop("Meta-analysis results must be an object of class `rma`.")
    }

  myData <- meta_retrieve(meta = meta, fixed = fixed)

  monte <- simulate(effect = meta$b,
                    effect_sd = meta$se,
                    weight = mean(myData$weights),
                    weight_sd = stats::sd(myData$weights),
                    studies = myData$k[1],
                    iter = 10000)

  calc_monte <- plot_simulation(simulated = monte,
                               observed_effects = myData$effect,
                               observed_weights = myData$weights,
                               meta_es = meta$b,
                               fileName)

  monte_p_values <- ma_monte_p_value(simulated = monte,
                                 observed_effects = myData$effect,
                                 observed_weights = myData$weights,
                                 meta_es = meta$b)

  return(monte_p_values)
}
