\name{Bayesian SCR}
\alias{Bsecr0}


\title{
Spatially explicit capture-recapture (secr) density estimation using MCMC
}
\description{
Functions to estimate density from mark-recapture data using MCMC methods and JAGS. 
}
\usage{
Bsecr0(capthist, buffer = 100, nAug = NA, 
         numSavedSteps=1e3, thinSteps=1, burnInSteps = 0) 
}
\arguments{
  \item{capthist}{
a \code{capthist} object as defined in package \code{secr} including capture data and detector (trap) layout }
  \item{buffer}{
scalar mask buffer radius (default 100 m)
}
  \item{nAug}{
number of individuals in the augmented population; if NA, a suitable default is chosen based on a preliminary run of \code{secr.fit}.
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
\details{
\code{Bsecr0} implements an intercept-only model (D ~ 1, g0 ~ 1, sigma ~ 1). 

}
\value{
Returns a data frame with one column for each parameter, ie. D, lam0 and sigma.

}

\references{
Efford ...

Royle et al ...
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
colMeans(tmp)

}
}
