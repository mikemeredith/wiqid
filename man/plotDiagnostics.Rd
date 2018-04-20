\name{diagnostic plots}
\alias{diagnostic plots}
\alias{diagPlot}
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

diagPlot(object, which, ask=TRUE, maxRow=4, RhatBad=1.05, ...)

tracePlot(object, ask=TRUE, ...)

densityPlot(object, ask=TRUE, ...)

acfPlot(object, lag.max=NULL, ask=TRUE, ...)

}
\arguments{
  \item{object}{
An object of class \code{Bwiqid}.
}
  \item{which}{
An optional vector of parameter names or numbers.
}
  \item{ask}{
If \code{TRUE}, the user will be prompted before the next page of output is displayed.
}
  \item{maxRow}{
Maximum number of rows to display in one window; each row consists of a trace plot
 and a density plot for one parameter.
}
  \item{RhatBad}{
Threshold for Rhat; parameters with \code{Rhat > RhatBad} are highlighted in red.
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
\donttest{
# Create a fake Bwiqid object:
fake <- data.frame(
  mu0 = rnorm(3000),        # normal, mean zero
  mu10 = rnorm(3000, 10),   # normal, mean 10
  sigma=rlnorm(3000),       # non-negative, skewed
  prob = rbeta(3000, 2,2),  # probability, central mode
  prob0 = rbeta(3000, 1,2), # probability, mode = 0
  N = rpois(3000, 20),      # large integers (no zeros)
  n = rpois(3000, 2),       # small integers (some zeros)
  const1 = rep(1, 3000),    # all values = 1
  const3.2 = rep(3.2, 3000))# all values the same but not integer
class(fake) <- c("Bwiqid", "data.frame")
attr(fake, "n.chains") <- 3
attr(fake, "Rhat") <- c(1, 1.01, 1.05, 1.1, 1.2, 1, NA, 1, NA)
attr(fake, "n.eff") <- c(2345, NA, 457, 63, 1234, 324, 543, 1, 1)

fake
diagPlot(fake)
diagPlot(fake, which=3:6)
diagPlot(fake, which="prob0")

par(mfrow=c(3,3))
tracePlot(fake, ask=FALSE)
densityPlot(fake, col=1:5, lwd=2, ask=FALSE)
acfPlot(fake, lag.max=10, ask=FALSE)
par(mfrow=c(1,1))
}
}
\keyword{hplot}