% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calc_ppp}
\alias{calc_ppp}
\title{Calculate Posterior Predictive p-values for meta-analysis results}
\usage{
calc_ppp(discrepancy, simulated_discrepancy)
}
\arguments{
\item{discrepancy}{A numeric value of observed discrepancies between the
obtained study-level effects.}

\item{simulated_discrepancy}{A vector of the simulated discrepancies of the
simulated study-level effects.}
}
\value{
mcmc A vector containing the Postrior Predictive p-values for each of
the discrepancy measures.
}
\description{
This function calculates the Posterior Predictive p-values for the
meta-analysis results
}
