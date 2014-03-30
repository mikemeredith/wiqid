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
# Create a fake Bwiqid object:
fake <- data.frame(mu = rnorm(3000), sigma=rlnorm(3000))
class(fake) <- c("Bwiqid", "data.frame")
attr(fake, "n.chains") <- 3
fake

tracePlot(fake)
densityPlot(fake, col=1:5, lwd=2)
acfPlot(fake, lag.max=10)

}
\keyword{hplot}