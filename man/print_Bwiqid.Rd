\name{print.Bwiqid}
\alias{print.Bwiqid}
\title{
Print method for objects of class 'Bwiqid'
}
\description{
Prints summary statistics for the parameters and several MCMC diagnostic measures for objects of class \code{Bwiqid}.
}
\usage{
\method{print}{Bwiqid}(x, digits=4, ...)  
    
}
\arguments{
  \item{x}{
an object of class \code{Bwiqid}.
}
  \item{digits}{
the number of digits to print.
}
  \item{\dots}{
further arguments for the print function.
}
}
\details{
The function prints a table with a row for each parameter \emph{after} removing duplicated rows. Duplication usually arises because a covariate has only a few unique values.

There are columns for each of the following summary statistics ...

 \tabular{ll}{
  \code{mean} \tab the mean of each MCMC chain. \cr
  \code{sd} \tab  the standard deviation of each MCMC chain. \cr
  \code{median} \tab the median of each MCMC chain. \cr
  \code{HDIlo} and \code{HDIup} \tab the lower and upper values of a 95\% Highest Density Interval CrI for each MCMC chain. }

... and for some or all of the following diagnostics, depending on the MCMC engine used for fitting the model:

 \tabular{ll}{
  \code{n.eff} \tab  the effective sample sizes for the parameters adjusted for autocorrelation; for stable estimates of credible intervals this should be at least 10,000. See \code{\link[coda]{effectiveSize}}.\cr
  \code{MCerror} \tab the Monte Carlo errors for the parameters, expressed as a percentage of the standard deviation. Values less than 5\% are acceptable. \cr
  \code{Rhat} \tab the with potential scale reduction factors for the parameters, which is 1 on convergence and should be < 1.05 for all parameters. See \code{\link[coda]{gelman.diag}}. }

}
\value{
Returns \code{x} invisibly.
}
\author{
Mike Meredith.
}

\examples{
# See examples for dippers.
}
\keyword{print}
\keyword{methods}
