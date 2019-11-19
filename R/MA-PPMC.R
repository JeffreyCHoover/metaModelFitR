#' @title Posterior Predictive Model Check for Meta-Analysis Models
#' @description This function performs a posterior predictive model check on
#'   meta-analysis results from the metafor package to assess model fit. This
#'   function returns posterior predictive p-values (PPP-values) for the
#'   following discrepancy measures: mean study-level effects, minimum
#'   study-level effects, maximum study-level effects, mean study-level
#'   weights, minimum study-level weights, maximum study-level effects and
#'   pooled effect size. Graphical results are also saved for each of these
#'   discrepancy functions.
#' @param fileName A character variable that serves as the stem for all saved
#'   output files.
#' @param meta An 'rma' object that contains the results of a meta-analysis
#'   performed in the metafor package.
#' @param fixed Logical variable specifying whether the meta-analysis results
#'   were using a fixed-effect model (default) or a random-effects model.
#' @return ppp A vector containing the PPP-values.
#' @examples
#' \dontrun{
#' ppmc_ma(fileName = "re-meta_ppmc", meta = re_meta, fixed = FALSE)
#' }
#' @importFrom magrittr %>%
#' @export
ppmc_ma <- function(fileName = NULL, meta, fixed = FALSE)
{
  meta.vi <- NULL
  meta.yi <- NULL

  if(is.null(meta)) {
    stop("Meta-analysis results parameter is required.")
  } else if (!("rma" %in% class(meta))) {
    stop("Meta-analysis results must be an object of class `rma`.")
  }

  grDevices::graphics.off()

  theData = data.frame(meta$yi, meta$vi) %>%
    dplyr::mutate(weights = 1 / meta.vi) %>%
    dplyr::rename(effects = meta.yi)
  fileNameRoot = fileName

  if(fixed) {
  dataList = list(
    n = length(meta$yi),
    es = meta$b,
    es_se = meta$se,
    w = 1 / meta$vi
    )

  modelString = "
  model {
    w_sd <- sd(w)
    w_mean <- mean(w)
    w_sum <- sum(w)

    # priors
    for(i in 1:n) {
      weights[i] ~ dunif(min(w), max(w))
      y[i] ~ dnorm(es, es_se)
    }

    for(i in 1:n) {
      ES_ag[i] <- sum(y[i] * weights[i]) / sum(weights[i])
    }

    ES_agg <- mean(ES_ag)

    minEffect <- min(y)
    maxEffect <- max(y)
  }
  "

  } else {
    dataList = list(
      n = length(meta$yi),
      tau = meta$tau2,
      es = meta$b,
      es_v = meta$vi,
      w = 1 / (meta$vi + meta$tau2))

    modelString = "
  model {
    w_sd <- sd(w)
    es_sd <- mean(es_v)
    w_mean <- mean(w)
    tau_d <- max(.1, tau)

    # priors
    for(i in 1:n) {
      weights[i] ~ dunif(min(w), max(w))
      theta[i] ~ dnorm(es, tau_d)
      y[i] ~ dnorm(theta[i], es_sd)
    }

    for(i in 1:n) {
      ES_ag[i] <- sum(y[i] * weights[i]) / sum(weights[i])
    }

    ES_agg <- mean(ES_ag)

    minEffect <- min(y)
    maxEffect <- max(y)
  }
  "
  }

  writeLines(modelString , con="TEMPmodel.txt")

  parameters = c("y", "ES_agg")
  adaptSteps = 1000
  burnInSteps = 1000
  numSavedSteps=20000
  thinSteps=20
  nChains = 3
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

  jagsModel = rjags::jags.model( "TEMPmodel.txt" , data=dataList ,
                          n.chains=nChains , n.adapt=adaptSteps )
  cat( "Burning in the MCMC chain...\n" )
  stats::update( jagsModel , n.iter=burnInSteps )

  cat( "Sampling final MCMC chain...\n" )
  codaSamples = rjags::coda.samples( jagsModel , variable.names=parameters ,
                              n.iter=nIter , thin=thinSteps )
  if ( !is.null(fileNameRoot) ) {
    save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata",sep="") )
  }

  mcmcMat = as.matrix(codaSamples)

  diagMCMC( codaObject=codaSamples , parName="ES_agg", filename=fileNameRoot)

  disc <- sum(abs(meta$yi - rep(meta$b, length(meta$yi)))) / abs(meta$b)

  sim_data <- as.data.frame(mcmcMat) %>%
    dplyr::select(2:7)

  sim_disc <- as.data.frame(rowSums(abs(sim_data - mcmcMat[,1]))/
                              abs(mcmcMat[,1])) %>%
    dplyr::rename(sim_disc = 1)

  ppp <- calc_ppp(discrepancy = disc,
                  simulated_discrepancy = sim_disc)
  return(ppp)
}
