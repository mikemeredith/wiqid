\name{BetaDist}
\alias{BetaDist}
\alias{dbeta2}
\alias{pbeta2}
\alias{qbeta2}
\alias{rbeta2}
\alias{getBetaPar}

\title{
The Beta Distribution
}
\description{
Density, distribution function, quantile function and random generation for the Beta distribution with parameters \code{mean} and \code{sd}. These are wrappers for \code{stats::dbeta}, etc. \code{getBetaPar} returns the shape and rate parameters.
}
\usage{
dbeta2(x, mean, sd)
pbeta2(q, mean, sd, lower.tail=TRUE, log.p=FALSE)
qbeta2(p, mean, sd, lower.tail=TRUE, log.p=FALSE)
rbeta2(n, mean, sd)
getBetaPar(mean, sd)
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
\code{dbeta2} gives the density, \code{pbeta2} gives the distribution function, \code{qbeta2} gives the quantile function, and \code{rbeta2} generates random deviates.

code{getBetaPar} returns a named vector with the shape and rate parameters corresponding to \code{mean} and \code{sd}.
}
\author{
Mike Meredith; getBetaPars function by "Max" at CrossValidated \url{https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance}
}

\seealso{
See the \pkg{stats} functions \code{\link{dbeta}}, \code{\link{pbeta}}, \code{\link{qbeta}}, \code{\link{rbeta}}.
}
\examples{
## to do
}
\keyword{"distribution"}