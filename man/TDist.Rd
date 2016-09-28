\name{TDist}
\alias{TDist}
\alias{dt2}
\alias{pt2}
\alias{qt2}
\alias{rt2}
\title{
The Student t Distribution
}
\description{
Density, distribution function, quantile function and random generation for the t distribution with df degrees of freedom, using mean and sd instead of a t-statistic. These are wrappers for stats::dt, etc.
}
\usage{
dt2(x, mean, sd, df)
pt2(x, mean, sd, df, lower.tail=TRUE, log.p=FALSE)
qt2(p, mean, sd, df, lower.tail=TRUE, log.p=FALSE)
rt2(n, mean, sd, df)
}
\arguments{
  \item{x}{
vector of parameter values
}
  \item{mean}{
mean of the t-distribution
}
  \item{sd}{
standard deviation of the t-distribution
}
  \item{df}{
degrees of freedom
}
  \item{lower.tail}{
logical; if TRUE (default), cumulative probabilities up to x, otherwise, above x.
}
  \item{log.p}{
logical; if TRUE, probabilities p are given as log(p).
}
  \item{p}{
probability.
}
  \item{n}{
number of random draws required.
}
}
\value{
\code{dt} gives the density, \code{pt} gives the cumulative probability, \code{qt} gives the quantile function, and \code{rt} generates random deviates.
}
\author{
Mike Meredith
}

\seealso{
See the \pkg{stats} functions \code{\link{dt}}, \code{\link{pt}}, \code{\link{qt}}, \code{\link{rt}}.
}
\examples{
## to do
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"distribution"}