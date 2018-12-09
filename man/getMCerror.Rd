\name{getMCerror}
\alias{getMCerror}

\title{
MCMC error using the batch method
}
\description{
The Monte Carlo error is a result of using (pseudo-)random numbers to generate the MCMC chain: multiple runs with different random number sequences produce different chains with slightly different summary statistics. A simple way to estimate MC error is to divide values into batches and calculate the variance of the batch means. See Lunn et al (2013, p77) for details. There is a trade-off between number of batches and batch size, and we make these approximately equal (only approximately as the number of batches must be a multiple of the number of chains). The function gives the same values as \code{coda::batchSE} with an appropriate \code{batchSize}.
}

\usage{
getMCerror(object, n.chains, SDpc=FALSE)
}

\arguments{
  \item{object}{
a vector, matrix or data frame with MCMC output in columns
}
  \item{n.chains}{
scalar integer, the number of chains concatenated to form the columns of \code{x}.
}
  \item{SDpc}{
logical, if TRUE the MC error is expressed as a percentage of the standard deviation of the corresponding posterior.
}
}

\value{
If SDpc is FALSE (the default), a named vector with the estimates of MC error. If TRUE, the MC error as a percentage of the standard deviation of the posterior chain. A value <5\% of SD is adequate for most purposes, but <1.5\% is needed to properly estimate tail probabilities (Lunn et al 2013, p78-79).
}

\references{
Lunn, D., Jackson, C., Best, N., Thomas, A., & Spiegelhalter, D. (2013) \emph{The BUGS book: a practical introduction to Bayesian analysis}, Chapman and Hall.

Roberts, G.O. (1996). Markov chain concepts related to sampling algorithms. In \emph{Markov Chain Monte Carlo in practice} (eds W.R. Gilks, D.J. Spiegelhalter & S. Richardson). Chapman & Hall, London.
}

\author{
Mike Meredith
}


\examples{
# Get some output to use
data(salamanders)
y <- rowSums(salamanders)
( out <- BoccSS0(y, 5) )

getMCerror(out)
getMCerror(out, SDpc=TRUE)
}
