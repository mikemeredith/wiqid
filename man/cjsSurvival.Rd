\name{Survival (CJS)}
\alias{survCJS}
\alias{survCJSaj}
\alias{BsurvCJS}
\alias{ch2mArray}
\title{
Survival from recapture data
}
\description{
Calculation of apparent survival (accounting for recapture probability) from recapture data formated as an m-array, with time-dependent phi or p, possibly with covariates. \code{ch2mArray} converts capture history matrices to m-arrays.
}
\usage{
survCJS(mArray, model=list(phi~1, p~1), data = NULL, ci = 0.95)

survCJSaj(mArray, mArrayJ, model=list(phiA~1, phiJ~1, p~1), data=NULL, ci = 0.95)

BsurvCJS(mArray, model=list(phi~1, p~1), data = NULL,
    numSavedSteps=1e4, thinSteps=1, burnInSteps = 1e3) 

ch2mArray(CH, freq = 1)
}
\arguments{
  \item{mArray}{
matrix giving the recapture history in m-array format, with zeros (not NAs) in the lower triangle.
}
  \item{mArrayJ}{
for a model with different survival probabilities for juveniles, \code{mArray} refers to animals released as adults and \code{mArrayJ} to those released as juveniles.
}
  \item{ci}{
the required confidence interval.
}
  \item{model}{
a list of formulae symbolically defining a linear predictor for each parameter in terms of covariates.
}
  \item{data}{
a data frame with a row for each survival interval / recapture occasion and columns for each of the covariates used to estimate phi or p.
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
  \item{CH}{
a 1/0 matrix with capture histories with a row for each animal captured and a column for each capture occasion.
}
  \item{freq}{
a scalar or a vector of length \code{nrow(CH)} with the frequency of each capture history.
}
}
\details{
For n REcapture occasions (ie. n + 1 capture occasions total), the m-array has:

  = n rows corresponding to marked animals caught and released on each occasion except for the last (those animals will not be recatured anyway).

  = n columns with the number of animals recaptured on each occasion except the first, PLUS a last column with the number never recaptured.

 The lower triangle of the matrix is all zeros.


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

\code{ch2mArray} returns an m-array.
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
# The classic Lebreton et al dipper data set:
dippers <- matrix(c(
  11,2,0,0,0,0,9,
  0,24,1,0,0,0,35,
  0,0,34,2,0,0,42,
  0,0,0,45,1,2,32,
  0,0,0,0,51,0,37,
  0,0,0,0,0,52,46),nrow=6,byrow=TRUE)

survCJS(dippers)  # default is a phi(.) p(.) model
survCJS(dippers, phi ~ time) 

BsurvCJS(dippers)  # default is also a phi(.) p(.) model
BsurvCJS(dippers, phi ~ time)   # using in-built time covariate

dd <- data.frame(flood = c(0, 1, 1, 0, 0, 0))
survCJS(dippers, list(phi~flood, p~time), data=dd)
BsurvCJS(dippers, list(phi ~ flood, p ~ time), data=dd) 


# The following model is overparameterised: phi[6] and p[6] cannot
#    be estimated separately and the variance-covariance matrix
#    cannot be calculated:
survCJS(dippers, list(phi~time, p~time))

}
