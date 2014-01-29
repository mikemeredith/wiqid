\name{Bayesian SCR}
\alias{Bsecr0}


\title{
Spatially explicit capture-recapture (secr) density estimation using MCMC
}
\description{
Functions to estimate density from mark-recapture data using MCMC methods and JAGS. 
}
\usage{
Bsecr0(capthist, buffer = 100, start = NULL, nAug = NA,
                 nChains = 3, numSavedSteps=1e4, thinSteps=1, burnInSteps = 0, priorOnly=FALSE) 
}
\arguments{
  \item{capthist}{
a \code{capthist} object as defined in package \code{secr} including capture data and detector (trap) layout }
  \item{buffer}{
scalar mask buffer radius (default 100 m)
}
  \item{start}{
an optional object of class \code{secr}, ie, output from the \code{secr.fit} function in package \code{secr}; objects of other classes are silently ignored.
}
  \item{nAug}{
number of individuals in the augmented population; if NA, a suitable default is chosen based on the object passed to \code{start} or a preliminary run of \code{secr.fit}.
}
  \item{nChains}{
the number of MCMC chains to run.
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
  \item{priorOnly}{
if TRUE, the function produces random draws from the appropriate \emph{prior} distributions, with a warning.
}
}
\details{
\code{Bsecr0} implements an intercept-only model (D ~ 1, g0 ~ 1, sigma ~ 1). 

}
\value{
Returns an object of class \code{Bwiqid}, data frame with one column for each parameter, ie. D, lam0 and sigma.

There are print and plot methods for \code{Bwiqid}.
}

\references{
Borchers & Efford (2008) Spatially explicit maximum likelihood methods for capture-recapture studies \emph{Biometrics} 64, 377-385

Royle & Dorazio (2008) \emph{Hierarchical modeling and inference in ecology}. Academic Press
}

\author{
Mike Meredith
}
\seealso{
The function \code{secr.fit} in package \code{secr}.
}
\examples{
\donttest{
# The stoats data set in 'secr'
require(secr)
data(stoatDNA)
tmp <- Bsecr0(stoatCH, buffer=1000)
plot(tmp)

}
}
