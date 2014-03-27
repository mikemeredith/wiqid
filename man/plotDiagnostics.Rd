\name{diagnostic plots}
\alias{diagnostic plots}
\alias{tracePlot}
\alias{densityPlot}
\alias{acfPlot}
\title{
Diagnostic graphics for MCMC output
}
\description{
Display trace plots, density plots or autocorrelation plots for the chains in the MCMC output. Each chain is plotted with a different colour.
}
\usage{
tracePlot(object, ask=TRUE, ...)

densityPlot(object, ask=TRUE, ...)

acfPlot(object, lag.max=NULL, ask=TRUE, ...)

}
\arguments{
  \item{object}{
An object of class \code{Bwiqid}.
}
  \item{ask}{
If \code{TRUE}, the user will be prompted before the next page of output is displayed.
}
  \item{lag.max}{
Maximum lag at which to calculate the acf; see \code{\link{acf}}.
}
  \item{\dots}{
Additional graphical parameters.
}
}
\details{
}
\value{
Return NULL invisibly. Used for their plotting side-effects.
}
\author{
Mike Meredith
}

\seealso{
\code{\link{plotPost}} for a histogram and summary statistics.
}
\examples{
# TO DO. Need a data set which is Bwiqid output

}
\keyword{hplot}