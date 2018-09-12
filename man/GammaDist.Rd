\name{GammaDist}
\alias{GammaDist}
\alias{dgamma2}
\alias{pgamma2}
\alias{qgamma2}
\alias{rgamma2}
\alias{getGammaPar}
\title{
The Gamma Distribution
}
\description{
Density, distribution function, quantile function and random generation for the Gamma distribution with parameters \code{mean} and \code{sd}. These are wrappers for \code{stats::dgamma}, etc. \code{getGammaPar} returns the shape and rate parameters.
}
\usage{
dgamma2(x, mean, sd)
pgamma2(q, mean, sd, lower.tail=TRUE, log.p=FALSE)
qgamma2(p, mean, sd, lower.tail=TRUE, log.p=FALSE)
rgamma2(n, mean, sd)
getGammaPar(mean, sd)
}
\arguments{
  \item{x}{
vector of parameter values
}
  \item{q}{
vector of quantiles
}
  \item{p}{
vector of provabilities
}
  \item{n}{
number of random draws required.
}
  \item{mean}{
mean of the t-distribution
}
  \item{sd}{
standard deviation of the t-distribution
}
  \item{lower.tail}{
logical; if TRUE (default), cumulative probabilities up to x, otherwise, above x.
}
  \item{log.p}{
logical; if TRUE, probabilities p are given as log(p).
}
}
\value{
\code{dgamma2} gives the density, \code{pgamma2} gives the distribution function, \code{qgamma2} gives the quantile function, and \code{rgamma2} generates random deviates.

code{getGammaPar} returns a named vector with the shape and rate parameters corresponding to \code{mean} and \code{sd}. }
\author{
Mike Meredith
}

\seealso{
See the \pkg{stats} functions \code{\link{dgamma}}, \code{\link{pgamma}}, \code{\link{qgamma}}, \code{\link{rgamma}}.
}
\examples{
## to do
}
\keyword{"distribution"}