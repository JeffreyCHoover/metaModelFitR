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
#' @export
ppmc_ma <- function(fileName, meta, fixed = FALSE)
{
  if(is.null(meta)) {
    stop("Meta-analysis results parameter is required.")
  } else if (!("rma" %in% class(meta))) {
    stop("Meta-analysis results must be an object of class `rma`.")
  }

  grDevices::graphics.off() # This closes all of R's graphics windows.

  theData = data.frame(meta$yi, meta$vi) %>%
    mutate(weights = 1 / meta.vi) %>%
    rename(effects = meta.yi)
  fileNameRoot = fileName # For output file names.

  # Package the data for JAGS:
  if(fixed) {
  dataList = list(
    n = length(meta$yi),
    es = meta$b,
    es_se = meta$se,
    w = 1 / meta$vi
    )

  # Define the JAGS model:
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
  " # close quote for modelString

  } else {
    dataList = list(
      n = length(meta$yi),
      tau = meta$tau2,
      es = meta$b,
      es_v = meta$vi,
      w = 1 / (meta$vi + meta$tau2))

    # Define the JAGS model:
    modelString = "
  model {
    w_sd <- sd(w)
    es_sd <- mean(es_v)
    w_mean <- mean(w)
    tau_d <- tau#max(1, tau)

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
  " # close quote for modelString
  }

  writeLines(modelString , con="TEMPmodel.txt")

  # Run the chains:
  parameters = c("y", "ES_agg")
  adaptSteps = 1000
  burnInSteps = 1000
  numSavedSteps=20000
  thinSteps=20
  nChains = 3
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

  jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , #inits=initsList ,
                          n.chains=nChains , n.adapt=adaptSteps )
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel , variable.names=parameters ,
                              n.iter=nIter , thin=thinSteps )
  if ( !is.null(fileNameRoot) ) {
    save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata",sep="") )
  }

  #save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )
  mcmcMat = as.matrix(codaSamples)

  # Examine the chains:
  # Convergence diagnostics:
  diagMCMC( codaObject=codaSamples , parName="ES_agg", filename=fileNameRoot)

  disc <- sum(abs(meta$yi - rep(meta$b, length(meta$yi)))) / abs(meta$b)

  sim_data <- as.data.frame(mcmcMat) %>%
    select(2:7)

  sim_disc <- as.data.frame(rowSums(abs(sim_data - mcmcMat[,1]))/
                              abs(mcmcMat[,1])) %>%
    rename(sim_disc = 1)

  # # Posterior descriptives:
  ppp <- calc_ppp(mcmc = as.data.frame(mcmcMat),
                  discrepancy = disc,
                  simulated_discrepancy = sim_disc)

  #graphs <- plot_ppmc(mcmc = as_tibble(mcmcMat),
  #                    discrepancy = disc,
  #                    simulated_discrepancy = sim_disc,
  #                    observed_effects = meta$yi,
  #                    observed_weights = 1 / meta$vi,
  #                    meta = meta,
  #                    filename = fileNameRoot)

  return(ppp)
}
