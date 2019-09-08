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
    tau = 0.1,
    es = meta$b,
    es_se = meta$se,
    w = 1 / meta$vi,
    es_list = meta$yi
    )
  } else {
    dataList = list(
      n = length(meta$yi),
      tau = 1,
      es = meta$b,
      es_se = meta$se,
      w = 1 / meta$vi,
      es_list = meta$yi)
  }

  # Define the JAGS model:
  modelString = "
  model {
    w_sd <- sd(w)
    es_sd <- sd(es_list)
    w_mean <- mean(w)
    w_min <- min(w)
    w_max <- max(w)

    # priors
    tau_d ~ dunif(0.05, tau)
    for(i in 1:n) {
      weights[i] ~ dnorm(w_mean, w_sd)
      effects[i] ~ dnorm(es, es_sd)
      theta[i] ~ dnorm(es, tau_d)
      y[i] ~ dnorm(theta[i], es_se)
    }

    #ES <- sum(FinalEffect * FinalWeight) / sum(FinalWeight)
    for(i in 1:n) {
      ES_ag[i] <- sum(y[i] * weights[i]) / sum(weights[i])
    }

    ES_agg <- mean(ES_ag)

    minWeight <- min(weights)
    maxWeight <- max(weights)
    meanWeight <- mean(weights)
    minEffect <- min(y)
    maxEffect <- max(y)
  }
  " # close quote for modelString
  writeLines(modelString , con="TEMPmodel.txt")

  # Run the chains:
  parameters = c("minEffect", "maxEffect", "minWeight",
                 "meanWeight", "maxWeight", "ES_agg")
  adaptSteps = 1000
  burnInSteps = 1000
  numSavedSteps=20000
  thinSteps=20
  nChains = 3

  runJagsOut <- runjags::run.jags( method=c("rjags","parallel")[2] ,
                                  model="TEMPmodel.txt" ,
                                  monitor=parameters ,
                                  data=dataList ,
                                  #inits=initsList ,
                                  n.chains=nChains ,
                                  adapt=adaptSteps ,
                                  burnin=burnInSteps ,
                                  sample=ceiling(numSavedSteps/nChains) ,
                                  thin=thinSteps ,
                                  summarise=FALSE ,
                                  plots=FALSE )
  codaSamples = coda::as.mcmc.list( runJagsOut )

  save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )
  mcmcMat = as.matrix(codaSamples)

  # Examine the chains:
  # Convergence diagnostics:
  diagMCMC( codaObject=codaSamples , parName="minEffect", filename=fileNameRoot)
  #diagMCMC( codaObject=codaSamples , parName="meanEffect", filename=fileNameRoot)
  diagMCMC( codaObject=codaSamples , parName="maxEffect", filename=fileNameRoot)
  diagMCMC( codaObject=codaSamples , parName="minWeight", filename=fileNameRoot)
  diagMCMC( codaObject=codaSamples , parName="meanWeight", filename=fileNameRoot)
  diagMCMC( codaObject=codaSamples , parName="maxWeight", filename=fileNameRoot)
  diagMCMC( codaObject=codaSamples , parName="ES_agg", filename=fileNameRoot)

  # # Posterior descriptives:

  cexMain=1.25
  cexMax = 2.5
  cexMin = 1.0
  #cexSlope1 = (cexMax-cexMin)/(max(theData$effect)-min(theData$effect))
  #cexSlope2 = (cexMax-cexMin)/(max(theData$weights)-min(theData$weights))
  #cexInter1 = cexMax - cexSlope1*max(theData$effect)
  #cexInter2 = cexMax - cexSlope2*max(theData$weights)

  #------------------------------------------------------------------------------

  #openGraph(height=2.75,width=5)
  #layout(matrix(1:2,ncol=2))
  #par( mar=c(3.0,0.5,3.0,0.5) , mgp=c(2.0,0.7,0) )
  #plotPost( mcmcMat[,"minEffect"] , xlab=bquote("Minimum Effect"),
  #          xlim = c(min(mcmcMat[, "minEffect"]), max(mcmcMat[, "minEffect"])),
  #          main="Minimum Effect" , cex.main=cexMain,
  #          compVal = min(theData$effect),
  #          ROPE = c(min(theData$effect) - .1, min(theData$effect) + .1))
  #plotPost( mcmcMat[,"multTSD"] , xlab=bquote(sigma*" of "*mu[j]) ,
  #          main="SD of Multiplier for Treatment" ,  cex.main=cexMain  )
  #graphName = paste0(fileNameRoot,"-minEffect")
  #saveGraph( file=graphName , type="eps" )
  #saveGraph( file=graphName , type="png" )

  #------------------------------------------------------------------------------

  #openGraph(height=2.75,width=5)
  #layout(matrix(1:2,ncol=2))
  #par( mar=c(3.0,0.5,3.0,0.5) , mgp=c(2.0,0.7,0) )
  #plotPost( mcmcMat[,"meanEffect"] , xlab=bquote("Mean Effect"),
  #          xlim = c(min(mcmcMat[,"meanEffect"]), max(mcmcMat[,"meanEffect"])),
  #          main="Mean Effect" , cex.main=cexMain,
  #          compVal = mean(theData$effect),
  #          ROPE = c(mean(theData$effect) - .1, mean(theData$effect) + .1))
  #plotPost( mcmcMat[,"multTSD"] , xlab=bquote(sigma*" of "*mu[j]) ,
  #          main="SD of Multiplier for Treatment" ,  cex.main=cexMain  )
  #graphName = paste0(fileNameRoot,"-meanEffect")
  #saveGraph( file=graphName , type="eps" )
  #saveGraph( file=graphName , type="png" )

  #------------------------------------------------------------------------------

  #openGraph(height=2.75,width=5)
  #layout(matrix(1:2,ncol=2))
  #par( mar=c(3.0,0.5,3.0,0.5) , mgp=c(2.0,0.7,0) )
  #plotPost( mcmcMat[,"maxEffect"] , xlab=bquote("Maximum Effect"),
  #          xlim = c(min(mcmcMat[,"maxEffect"]), max(mcmcMat[,"maxEffect"])),
  #          main="Maximum Effect" , cex.main=cexMain,
  #          compVal = max(theData$effect),
  #          ROPE = c(max(theData$effect) - .1, max(theData$effect) + .1))
  #plotPost( mcmcMat[,"multTSD"] , xlab=bquote(sigma*" of "*mu[j]) ,
  #          main="SD of Multiplier for Treatment" ,  cex.main=cexMain  )
  #graphName = paste0(fileNameRoot,"-maxEffect")
  #saveGraph( file=graphName , type="eps" )
  #saveGraph( file=graphName , type="png" )

  #------------------------------------------------------------------------------

  #openGraph(height=2.75,width=5)
  #layout(matrix(1:2,ncol=2))
  #par( mar=c(3.0,0.5,3.0,0.5) , mgp=c(2.0,0.7,0) )
  #plotPost( mcmcMat[,"minWeight"] , xlab=bquote("Minimum Weight"),
  #          xlim = c(-2, max(mcmcMat[, "minWeight"])),
  #          main="Minimum Weight" , cex.main=cexMain,
  #          compVal = min(theData$weights),
  #          ROPE = c(min(theData$weights) - 5, min(theData$weights) + 5))
  #plotPost( mcmcMat[,"multTSD"] , xlab=bquote(sigma*" of "*mu[j]) ,
  #          main="SD of Multiplier for Treatment" ,  cex.main=cexMain  )
  #graphName = paste0(fileNameRoot,"-minWeight")
  #saveGraph( file=graphName , type="eps" )
  #saveGraph( file=graphName , type="png" )

  #------------------------------------------------------------------------------

  #openGraph(height=2.75,width=5)
  #layout(matrix(1:2,ncol=2))
  #par( mar=c(3.0,0.5,3.0,0.5) , mgp=c(2.0,0.7,0) )
  #plotPost( mcmcMat[,"meanWeight"] , xlab=bquote("Mean Weight"),
  #          xlim = c(-2, max(mcmcMat[,"meanWeight"])),
  #          main="Mean Weight" , cex.main=cexMain,
  #          compVal = mean(theData$weights),
  #          ROPE = c(mean(theData$weights) - 5, mean(theData$weights) + 5))
  #plotPost( mcmcMat[,"multTSD"] , xlab=bquote(sigma*" of "*mu[j]) ,
  #          main="SD of Multiplier for Treatment" ,  cex.main=cexMain  )
  #graphName = paste0(fileNameRoot,"-meanWeight")
  #saveGraph( file=graphName , type="eps" )
  #saveGraph( file=graphName , type="png" )

  #------------------------------------------------------------------------------

  #openGraph(height=2.75,width=5)
  #layout(matrix(1:2,ncol=2))
  #par( mar=c(3.0,0.5,3.0,0.5) , mgp=c(2.0,0.7,0) )
  #plotPost( mcmcMat[,"maxWeight"] , xlab=bquote("Maximum Weight"),
  #          xlim = c(-2, max(mcmcMat[,"maxWeight"])),
  #          main="Maximum Weight" , cex.main=cexMain,
  #          compVal = max(theData),
  #          ROPE = c(max(theData$weights) - 5, max(theData$weights) + 5))
  #plotPost( mcmcMat[,"multTSD"] , xlab=bquote(sigma*" of "*mu[j]) ,
  #          main="SD of Multiplier for Treatment" ,  cex.main=cexMain  )
  #graphName = paste0(fileNameRoot,"-maxWeight")
  #saveGraph( file=graphName , type="eps" )
  #saveGraph( file=graphName , type="png" )

  #------------------------------------------------------------------------------

  #openGraph(height=2.75,width=5)
  #layout(matrix(1:2,ncol=2))
  #par( mar=c(3.0,0.5,3.0,0.5) , mgp=c(2.0,0.7,0) )
  #plotPost( mcmcMat[,"ES_agg"] , xlab=bquote("Aggregate Effect Size"),
  #          xlim = c(-10, 10),
  #          main="Aggregate Effect Size" , cex.main=cexMain,
  #          compVal = max(meta$b),
  #          ROPE = c(max(meta$b) - .1, max(meta$b) + .1))
  #plotPost( mcmcMat[,"multTSD"] , xlab=bquote(sigma*" of "*mu[j]) ,
  #          main="SD of Multiplier for Treatment" ,  cex.main=cexMain  )
  #graphName = paste0(fileNameRoot,"-ES_agg")
  #saveGraph( file=graphName , type="eps" )
  #saveGraph( file=graphName , type="png" )

  #------------------------------------------------------------------------------

  ppp <- calc_ppp(mcmc = as.data.frame(mcmcMat),
                  observed_effects = theData$effect,
                  observed_weights = theData$weights)

  graphs <- plot_ppmc(mcmc = as_tibble(mcmcMat),
                      observed_effects = theData$effect,
                      observed_weights = theData$weights,
                      meta = meta,
                      filename = fileNameRoot)

  print(ppp)
  return(ppp)
}
