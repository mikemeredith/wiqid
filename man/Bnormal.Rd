\name{Bayesian normal estimation}
\alias{Bnormal}


\title{
Bayesian modelling of a normal (Gaussian) distribution
}
\description{
Bayesian estimation of mean and standard deviation of a normal distribution based on a sample using Gibbs sampling implemented in R.

}
\usage{
Bnormal(x, priors=list(), n.chains=3, n.iter=10100, n.burnin=100)
 
}
\arguments{
  \item{x}{
a vector (length > 1) with observed sample values; missing values not allowed.
}
  \item{priors}{
a list with elements specifying the priors for the mean and variance; see Details. Note that the default prior for the mean = 0, which will often be inappropriate. If \code{priors = NULL}, improper flat priors are used.
}
  \item{n.chains}{
number of MCMC chains to run.
}
  \item{n.iter}{
number of values per chain, including burn-in.
}
  \item{n.burnin}{
number of values per chain to discard as burn-in.
}
}

\details{
The function uses a Gibbs sampler, sampling alternately from the conjugate posterior distributions for the mean with known variance, and the variance with known mean: see Gill (2008:74-79). This allows for separate specification of the priors for the mean and standard deviation of the Gaussian distribution, with values in the \code{priors} argument list. See Ntzoufras (2009:72).

The prior for the mean is a normal distribution with mean \code{mMean} and standard deviation \code{sMean}. Note that \code{sMean} is \emph{not} the standard deviation of the modelled distribution, but reflects the 'width' of the prior for the mean. The default is \code{mMean = 0} and \code{sMean = 100}. A warning is issued if the sample mean is outside the range \code{mMean +- sMean}.

The prior for the variance is an inverse gamma function with parameters \code{aSigma} and \code{bSigma}. The default is the uninformative prior, \code{aSigma = bSigma = 0.001}.
}

\value{
Returns an object of class \code{Bwiqid}, which is a data frame with a column for each parameter in the model.

There are \code{print} and \code{plot} methods for class \code{Bwiqid}, as well as \link{diagnostic plots}.
}

\references{

Ntzoufras 2009. \emph{Bayesian modeling using WinBUGS}. Wiley, Hoboken NJ.

Gill 2008. \emph{Bayesian methods: A social and behavioral sciences approach}. 2nd edition, Chapman and Hall/CRC.
}

\author{
Mike Meredith, after Ntzoufras (2009) and Brian Neelon's web site.
}

\examples{
# Generate a sample from a normal distribution, maybe the head-body length of a
#   carnivore in meters:
HB <- rnorm(10, 0.9, 0.15)
Bnormal(HB)
# The default prior for the mean, Normal(0, 100^2), works for this, but not if
#   we measure in mm:
HBmm <- HB * 1000
Bnormal(HBmm)  # warning
Bnormal(HBmm, priors=list(mMean=1000, sMean=500)) # more appropriate prior


}
