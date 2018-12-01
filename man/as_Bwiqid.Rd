\name{Bwiqid-class}
\alias{Bwiqid-class}
\alias{as.Bwiqid}
\alias{as.Bwiqid.default}
\alias{as.Bwiqid.data.frame}
\alias{as.Bwiqid.mcmc.list}
\alias{as.Bwiqid.bugs}
\alias{as.Bwiqid.rjags}
\alias{as.Bwiqid.runjags}

\title{
Conversion to class Bwiqid
}
\description{
Convert output containing MCMC chains to the class \code{Bwiqid}. The function is generic, with methods for a range of input objects.
}
\usage{
as.Bwiqid(object, ...)

\method{as.Bwiqid}{default}(object, ...)

\method{as.Bwiqid}{data.frame}(object, header, defaultPlot,
  n.chains=1, Rhat=TRUE, n.eff=TRUE, ...)

\method{as.Bwiqid}{mcmc.list}(object, header, defaultPlot, ...)

\method{as.Bwiqid}{bugs}(object, header, defaultPlot, ...)

\method{as.Bwiqid}{rjags}(object, header, defaultPlot, ...)

\method{as.Bwiqid}{runjags}(object, header, defaultPlot, ...)

}

\arguments{
  \item{object}{
an object containing the MCMC chains; see Details.
}
  \item{header}{
an optional character object with the text to be displayed by \code{print.Bwiqid}
}
  \item{defaultPlot}{
an optional character object specifying which column should be plotted by default.
}
  \item{n.chains}{
integer, the number of chains included in the object.
}
  \item{Rhat}{
either logical or a vector of Rhat values of length equal to the number of columns in \code{object}; if \code{Rhat=TRUE}, \code{\link{simpleRhat}} will be called to provide values.
}
  \item{n.eff}{
either logical or a vector of effective sample sizes of length equal to the number of columns in \code{object}; if \code{n.eff=TRUE}, \code{\link{effectiveSize}} will be called to provide values.
}
  \item{\dots}{
named parameters to be passed to other methods.
}
}

\details{
Several R packages provide interfaces with WinBUGS, OpenBUGS or JAGS, and have their own classes and methods for printing, plotting, etc. \code{as.Bwiqid} allows all of these to be converted to a common class.
}

\value{
An object of class \code{Bwiqid}. This is a data frame with a column for the MCMC chain for each parameter. Column names are stripped of square brackets, and commas are replaced with dots, so (eg) \code{p[1, 1]} becomes \code{p1.1}.

Any of the following attributes may be present:

  \tabular{ll}{
  \code{header} \tab text to be displayed as the first line when the object is printed. \cr
  \code{call} \tab the original function call. \cr
  \code{n.chains} \tab  the number of chains. \cr
  \code{n.eff} \tab  a vector of effective sample sizes for the parameters. \cr
  \code{MCerror} \tab a vector with Monte Carlo errors for the parameters. \cr
  \code{Rhat} \tab a vector with potential scale reduction factors for the parameters. \cr
  \code{defaultPlot} \tab the name of the parameter to be plotted as the default by \code{plot.Bwiqid}. \cr
  \code{timetaken} \tab an object of class \code{difftime} with the time taken for the MCMC run. }

}

\author{
Mike Meredith.
}


\keyword{methods}

