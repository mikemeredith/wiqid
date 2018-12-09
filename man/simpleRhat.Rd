\name{simpleRhat}
\alias{simpleRhat}

\title{
The Brooks-Gelman-Rubin (BGR) convergence diagnostic
}
\description{
An 'interval' estimator of the 'potential scale reduction factor' (Rhat) for MCMC output. Similar to the function \code{\link{gelman.diag}} in \pkg{coda}, but much faster when thousands of parameters are involved and will not cause R to crash.
}

\usage{
simpleRhat(object, n.chains, burnin=0)
}

\arguments{
  \item{object}{
a vector, matrix or data frame with MCMC output in columns
}
  \item{n.chains}{
scalar integer, the number of chains concatenated to form the columns of \code{x}; multiple chains are required to calculate Rhat; if \code{n.chains} is missing and \code{object} has an attribute \code{n.chains}, that value will be used.
}
  \item{burnin}{
scalar, between 0 and 0.9; the proportion of values at the start of each chain to discard as burn-in.
}
}

\details{
Calculates the Gelman-Rubin convergence statistic, as modified by Brooks and Gelman (1998). Following the WinBUGS User Manual (Spiegelhalter et al, 2003), we use the width of the central 80\% interval as the measure of width. Rhat is the ratio of the width of the pooled chains to the mean width of the individual chains. At convergence, this should be 1; values less than 1.1 or 1.05 are considered satisfactory.
}

\value{
A named vector with the Rhat values.
}

\references{
Brooks, S.P. & Gelman, A. (1998) General methods for monitoring convergence of iterative simulations. \emph{Journal of Computational and Graphical Statistics}, 7, 434-455.

Spiegelhalter, Thomas, Best & Lunn (2003) WinBUGS User Manual Version 1.4, on line \href{https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/manual14.pdf}{here}.
}

\author{
Mike Meredith
}


\examples{
# Get some output to use
data(salamanders)
y <- rowSums(salamanders)
( out <- BoccSS0(y, 5) )

simpleRhat(out)
}
