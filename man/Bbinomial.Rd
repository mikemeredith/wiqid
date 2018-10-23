\name{Bayesian binomial simulation}
\alias{Bbinom}
\alias{Bbinomial}


\title{
Bayesian analysis of binomial data
}
\description{
Simulates a sample from the posterior for a binomial likelihood and beta prior.
}
\usage{
Bbinom(y, n, priors=NULL, sample=10000)
Bbinomial(y, n, priors=NULL, sample=100000)
}
\arguments{
  \item{y}{
the number of successes
}
  \item{n}{
the number of trials
}
  \item{priors}{
an optional list with elements specifying the priors for the mode and concentration of the beta prior distribution; see Details.
}
  \item{sample}{
the number of MCMC observations to be returned.
}
}

\details{
The function generates a vector of random draws from the posterior distribution of the probability of success. It uses conjugacy to determine the parameters of the posterior beta distribution, and draws independent samples from this.

A prior can be specified with the \code{priors} argument. A beta prior is used, specified by mode, \code{mode}, and concentration, \code{conc}.

When \code{priors = NULL} (the default), a uniform prior corresponding to beta(1, 1) is used.
}

\value{
Returns an object of class \code{Bwiqid}, which is a data frame with a column for each parameter in the model.

There are \code{print} and \code{plot} methods for class \code{Bwiqid}, as well as \link{diagnostic plots}.
}

\author{
Mike Meredith.
}

\examples{
# Generate a sample from a binomial distribution, maybe the number of sites
#   where a species was detected:
n <- 10  # number of trials (sites visited)
( y <- rbinom(1, n, 0.75) ) # number of successes (sites where species detected)
Bbinomial(y, n)  # with uniform prior
Bbinomial(y, n, priors=list(mode=0.4, conc=5))  # with informative prior
}

