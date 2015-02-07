\name{Bayesian regression}
\alias{Bregress}
\alias{Blogit}


\title{
Bayesian Gaussian and logistic regression
}
\description{
These are wrappers for the \pkg{MCMCpack} functions \code{\link{MCMCregress}} and  \code{\link{MCMClogit}}, with the same calling syntax as other Bayesian functions in the \pkg{wiqid} package, and returning an object of class \code{Bwiqid}. 

}
\usage{

Bregress(formula, data = NULL, priors=list(),
              n.iter=11000, n.burnin=1000, ...)
              
Blogit(formula, data = NULL, priors=list(),
              n.iter=11000, n.burnin=1000, ...)
              
}
\arguments{
  \item{formula}{
a symbolic description of the model to be fitted; see \code{\link{formula}}. 
}
  \item{data}{
a data frame containing the variables in the model. 
}
  \item{priors}{
a list with elements specifying the priors for the mean and variance; see Details. 
}
  \item{n.iter}{
number of values per chain, including burn-in.
}
  \item{n.burnin}{
number of values per chain to discard as burn-in.
}
  \item{\dots}{
other named arguments passed to \code{MCMCregress} or \code{MCMClogit}.
}
}

\details{
These are wrappers for the corresponding \code{MCMC*} functions.

The coefficients are given normal (Gaussian) priors with means \code{b0} and precision \code{B0}. \code{b0} can be specified as a scalar, when the same value will be used for all coefficients, or as a vector of length K, where K is the number of coefficients, including the intercept; the default is 0. \code{B0} can be specified as a scalar, or as a K x K matrix; the default, \code{B0 = 0}, is an improper flat prior.

The residual variance is given an inverse Gamma prior: \code{c0} is 2 x shape parameter and \code{d0} is 2 x scale parameter. This corresponds to \code{c0} observations with precision \code{d0}.

The default is equivalent to \code{prior = list(b0=0, B0=0, c0=0.001, d0=0.001)}.
}

\value{
Returns an object of class \code{Bwiqid}, which is a data frame with a column for each parameter in the model.

There are \code{print} and \code{plot} methods for class \code{Bwiqid}, as well as \link{diagnostic plots}.
}

\author{
Wrappers by Mike Meredith.
}

\examples{
# 

}
