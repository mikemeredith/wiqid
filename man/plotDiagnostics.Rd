\name{diagnostic plots}
\alias{diagnostic plots}
\alias{diagPlot}
\alias{tracePlot}
\alias{densityPlot}
\alias{acfPlot}
\alias{crosscorrPlot}
\title{
Diagnostic graphics for MCMC output
}
\description{
Display trace plots, density plots or autocorrelation plots for the chains in the MCMC output. Each chain is plotted with a different colour.
}
\usage{

diagPlot(object, which, howMany, ask=TRUE, maxRow=4, RhatBad=1.05, ...)

tracePlot(object, ask=TRUE, ...)

densityPlot(object, ask=TRUE, ...)

acfPlot(object, lag.max=NULL, ask=TRUE, ...)

crosscorrPlot(object, which, ...)

}
\arguments{
  \item{object}{
An object of class \code{Bwiqid}.
}
  \item{which}{
An optional vector of parameter names or numbers.
}
  \item{howMany}{
How many iterations to plot; if negative, the iterations at the end of the chains will be plotted.
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
Additional graphical parameters passed to \code{\link{plot.default}} or for \code{crosscorrPlot} to \code{\link{corrplot}}.
}
}

\value{
\code{crosscorrPlot} returns the correlation matrix invisibly; other functions return NULL invisibly. Used for their plotting side-effects.
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
  mu0 = rnorm(3000),         # normal, mean zero
  mu10 = rnorm(3000, rep(9:11, each=1000), 1),
                             # normal, mean 10, but poor mixing
  sigma=rlnorm(3000),        # non-negative, skewed
  prob = rbeta(3000, 2, 2),  # probability, central mode
  prob0 = rbeta(3000, 1,2),  # probability, mode = 0
  N = rpois(3000, rep(c(24, 18, 18), each=1000)),
                             # large integers (no zeros), poor mixing
  n = rpois(3000, 2),        # small integers (some zeros)
  const1 = rep(1, 3000),     # all values = 1
  const3.2 = rep(10/3, 3000))# all values the same but not integer
class(fake) <- c("Bwiqid", "data.frame")
attr(fake, "n.chains") <- 3
attr(fake, "Rhat") <- c(1, 1.3, 1.05, 0.98, 1, 1.2, NA, 1, NA)
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

crosscorrPlot(fake)
}
}
\keyword{hplot}