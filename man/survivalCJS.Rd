\name{Survival (CJS)}
\alias{survCJS}
\alias{survCJSaj}
\alias{BsurvCJS}
\title{
Survival from recapture data
}
\description{
Calculation of apparent survival (accounting for recapture probability) from mark-recapture data, with time-dependent phi or p, possibly with covariates. Function \code{survCHSaj} allows for different survival parameters for juveniles and adults.
}
\usage{
survCJS(DH, model=list(phi~1, p~1), data=NULL, freq=1, ci = 0.95)

survCJSaj(DHj, DHa=NULL, model=list(phiJ~1, phiA~1, p~1), data=NULL,
    freqj=1, freqa=1, ci = 0.95)

BsurvCJS(DH, model=list(phi~1, p~1), data = NULL, freq=1,
    numSavedSteps=1e4, thinSteps=1, burnInSteps = 1e3) 
}
\arguments{
  \item{DH}{
a 1/0 matrix with detection histories with a row for each animal captured and a column for each capture occasion.
}
  \item{model}{
a list of formulae symbolically defining a linear predictor for each parameter in terms of covariates.
}
  \item{data}{
a data frame with a row for each survival interval / recapture occasion and columns for each of the covariates used to estimate phi or p.
}
  \item{freq}{
a scalar or a vector of length \code{nrow(DH)} with the frequency of each detection history.
}
  \item{DHj, DHa}{
detection history matrices for animals marked a juveniles and adults respectively; DHa should be NULL if no animals were marked as adults.
}
  \item{freqj, freqa}{
frequencies of each detection history in DHj and DHa; freqa is ignored if DHa = NULL.
}
  \item{ci}{
the required confidence interval.
}
  \item{numSavedSteps}{
the number of MCMC observations to be returned.
}
  \item{thinSteps}{
thinning rate. If set to n > 1, one value in n steps of the MCMC chain is returned. This is useful if autocorrelation is high.
}
  \item{burnInSteps}{
number of steps to discard as burn-in at the beginning of the chain.
}
}

\value{
\code{survCJS} returns an object of class \code{wiqid}, a list with elements:

  \item{call}{The call used to produce the results}
  \item{beta }{Estimates of the coefficients in the linear predictors for phi and p.}
  \item{beta.vcv }{The variance-covariance matrix for the beta estimates.}
  \item{real }{Back-transformed estimates of phi and p for each interval / occasion. }
  \item{logLik}{a vector with elements for log(likelihood), number of parameters, and effective sample size. If parameters \emph{and their SEs} cannot be estimated, the first element should be \code{NA}.}
  
There are \code{print}, \code{logLik}, and \code{nobs} methods for class \code{wiqid}.

\code{BsurvCJS} returns an object of class \code{Bwiqid}, a data frame with columns for each p and psi value containing the series of MCMC samples, and attributes for \code{effectiveSize} and \code{defaultPlot}.
}
\references{
Lebreton, J-D; K P Burnham; J Clobert; D R Anderson. 1992. Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. \emph{Ecological Monographs} 62:67-118.
}
\author{
Mike Meredith
}

\section{Benchmarks}{
Output of \code{survCJS} has been checked against program MARK with the dipper data set: coefficients are not the same as MARK uses models without an intercept, but the real values agree to 3 decimal places.
}

\examples{
# New examples to do; need a suitable data set, or add dippers to wiqid.

}
