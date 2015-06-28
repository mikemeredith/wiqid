\name{Bayesian normal estimation}
\alias{Bnormal}
\alias{Bnormal2}


\title{
Bayesian modelling of a normal (Gaussian) distribution
}
\description{
Bayesian estimation of centre (\eqn{\mu}) and scale (\eqn{\sigma}) of a normal distribution based on a sample. \code{Bnormal} uses a gamma prior on the precision, \eqn{\tau = 1/\sigma^2}, while \code{Bnormal2} applies a gamma prior to \eqn{\sigma}.
}
\usage{
Bnormal(y, priors=NULL,
    numChains=3, numSavedSteps=50000, burnInSteps=100)

Bnormal2(y, priors=NULL, doPriorsOnly=FALSE,
    numSavedSteps=1e5, thinSteps=1, burnInSteps = 1000,
    verbose=TRUE, rnd.seed=NULL, parallel=NULL)

}
\arguments{
  \item{y}{
a vector (length > 1) with observed sample values; missing values not allowed.
}
  \item{priors}{
an optional list with elements specifying the priors for the centre and scale; see Details.
}
 \item{doPriorsOnly}{
if TRUE, \code{Bnormal2} returns MCMC chains representing the prior distributions, \emph{not} the posterior distributions for your data set.
}
  \item{numChains}{
the number of MCMC chains to run.
}
  \item{numSavedSteps}{
the number of MCMC observations to be returned.
}
  \item{thinSteps}{
thinning rate. If set to n > 1, n steps of the MCMC chain are calculated for each one returned. This is useful if autocorrelation is high.
}
  \item{burnInSteps}{
number of steps to discard as burn-in at the beginning of each chain.
}
   \item{verbose}{
if FALSE, output to the R Console is suppressed.
}
  \item{rnd.seed}{
a positive integer (or NULL): the seed for the random number generator, used to obtain reproducible samples if required.
}
  \item{parallel}{
if NULL or TRUE and > 3 cores are available, the MCMC chains are run in parallel. (If TRUE and < 4 cores are available, a warning is given.)
}
}

\details{
The function generates vectors of random draws from the posterior distributions of the population centre (\eqn{\mu}) and scale (\eqn{\sigma}). \code{Bnormal} uses a Gibbs sampler inplemented in R, while \code{Bnormal2} uses JAGS (Plummer 2003).

Priors for all parameters can be specified by including elements in the \code{priors} list. For both functions, \eqn{\mu} has a normal prior, with mean \code{muMean} and standard deviation \code{muSD}. For \code{Bnormal}, a gamma prior is used for the precision, \eqn{\tau = 1\\sigma^2}, with parameters specified by \code{tauShape} and \code{tauRate}. For \code{Bnormal2}, a gamma prior is placed on \eqn{\sigma}, with parameters specified by mode, \code{sigmaMode}, and SD, \code{sigmaSD}.

When \code{priors = NULL} (the default), \code{Bnormal} uses improper flat priors for both \eqn{\mu} and \eqn{\tau}, while \code{Bnormal2} uses a broad normal prior (muMean = mean(y), muSD = sd(y)*5) for \eqn{\mu} and a uniform prior on (sd(y) / 1000, sd(y) * 1000) for \eqn{\sigma}.
}

\value{
Returns an object of class \code{Bwiqid}, which is a data frame with a column for each parameter in the model.

There are \code{print} and \code{plot} methods for class \code{Bwiqid}, as well as \link{diagnostic plots}.
}

\references{
Kruschke, J K. 2013. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} 142(2):573-603. doi: 10.1037/a0029146

Plummer, Martyn (2003). JAGS: A Program for Analysis of Bayesian Graphical Models Using Gibbs Sampling, \emph{Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003)}, March 20-22, Vienna, Austria. ISSN 1609-395X
}

\author{
Mike Meredith, \code{Bnormal} based on code by Brian Neelon, \code{Bnormal2} adapted from code by John Kruschke.
}

\examples{
# Generate a sample from a normal distribution, maybe the head-body length of a
#   carnivore in mm:
HB <- rnorm(10, 900, 15)
Bnormal(HB)  # with improper flat priors for mu and tau
Bnormal(HB, priors=list(muMean=1000, muSD=200))
Bnormal(HB, priors=list(muMean=1, muSD=0.2)) # a silly prior produces a warning.
\donttest{
Bnormal2(HB)  # with broad normal prior for mu, uniform for sigma
Bnormal2(HB, priors=list(muMean=1000, muSD=200, sigmaMode=20, sigmaSD=10))
} 
}
