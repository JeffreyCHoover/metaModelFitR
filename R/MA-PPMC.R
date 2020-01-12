#' @title Posterior Predictive Model Check for Meta-Analysis Models
#' @description This function performs a posterior predictive model check on
#'   meta-analysis results from the metafor package to assess model fit. This
#'   function returns posterior predictive p-values (PPP-values) for the
#'   following discrepancy measures: mean study-level effects, minimum
#'   study-level effects, maximum study-level effects, mean study-level
#'   weights, minimum study-level weights, maximum study-level effects and
#'   pooled effect size. Graphical results are also saved for each of these
#'   discrepancy functions.
#'
#' @param meta An 'rma' object that contains the results of a meta-analysis
#'   performed in the metafor package.
#' @param fixed Logical variable specifying whether the meta-analysis results
#'   were using a fixed-effect model (default) or a random-effects model.
#' @param iter A double specifying the number of iterations to be used in the
#'   posterior predictive model check.
#'
#' @return ppp A vector containing the PPP-values.
#'
#' @examples {
#' \dontrun{
#' library(Rcpp)
#' library(metafor)
#' d = c(.2, .3, .25, .3, .4, .15)
#' w = c(250, 300, 400, 275, 250, 200)
#' v = 1 / w
#'
#' re_meta = metafor::rma(yi = d, vi = v, method = "REML")
#'
#' ppmc_ma(meta = re_meta, fixed = FALSE)
#' }
#' }
#'
#' @importFrom magrittr %>%
#' @export
ppmc_ma <- function(meta, fixed = FALSE) {
  if (is.null(meta)) {
    stop ("Meta-analysis results parameter is required.")
  } else if (!("rma" %in% class(meta))) {
    stop ("Meta-analysis results must be an object of class `rma`.")
  }

  es <- as.numeric(meta$b)

  meta_data <- list(
    n = length(meta$yi),
    tau = max(.01, (meta$tau2)),
    es = es,
    es_sd = meta$se * sqrt(length(meta$yi)),
    w_mean = mean(1 / (meta$vi + meta$tau2)),
    w_sd = stats::sd(1 / (meta$vi + meta$tau2))
  )

  model_string <- "
  data {
    int<lower=0> n;
    real<lower=0> tau;
    real es;
    real w_mean;
    real<lower=0> w_sd;
    real<lower=0> es_sd;
  }
  parameters {
    real theta[n];
    real y[n];
    real weights[n];
  }
  model {
    weights ~ normal(w_mean, w_sd * 2);
    theta ~ normal(es, tau);
    y ~ normal(theta, es_sd);
  }
  "

  fit1 <- rstan::stan(
    model_code = model_string,
    data = meta_data,
    chains = 3,
    warmup = 1000,
    iter = 17000,
    cores = 2,
    refresh = 0
  )

  mcmc <- as.matrix(fit1)[, -ncol(as.matrix(fit1))]
  mcmc <- mcmc[, -c(1:length(meta$yi))]

  effects <- mcmc[, 1:(ncol (mcmc) / 2)]
  weights <- mcmc[, ( (ncol (mcmc) / 2) + 1):ncol (mcmc)]

  num <- rowSums (effects * weights)
  denom <- rowSums (weights)

  es_agg <- num / denom

  disc <- sum ( ( (meta$yi[1:length (meta$yi)] - rep (meta$beta[1],
                                                 length (meta$yi))) ^ 2))

  sim_disc <- vector (length = nrow(mcmc))

  for (ii in 1:nrow (mcmc)) {
    sim_disc[ii] <- sum ( (effects[ii, ] - rep (es_agg[ii],
                                             length (effects[ii, ]))) ^ 2)
  }

  ppp <- calc_ppp (discrepancy = disc,
                   simulated_discrepancy = sim_disc)
  return (ppp)
}
