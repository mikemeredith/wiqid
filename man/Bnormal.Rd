\name{Bayesian normal estimation}
\alias{Bnormal}


\title{
Bayesian modelling of a normal (Gaussian) distribution
}
\description{
Bayesian estimation of mean and standard deviation of a normal distribution based on a sample using JAGS.

}
\usage{
Bnormal(y, priors=NULL, doPriorsOnly=FALSE,
    numSavedSteps=1e5, thinSteps=1, burnInSteps = 1000,
    verbose=TRUE, rnd.seed=NULL, parallel=NULL)
 
}
\arguments{
  \item{y}{
a vector (length > 1) with observed sample values; missing values not allowed.
}
  \item{priors}{
an optional list with elements specifying the priors for the mean and variance; see Details.
}
 \item{doPriorsOnly}{
if TRUE, \code{Bnormal} returns MCMC chains representing the prior distributions, \emph{not} the posterior distributions for your data set.
}
  \item{numSavedSteps}{
the number of MCMC observations to be returned.
}
  \item{thinSteps}{
thinning rate. If set to n > 1, n steps of the MCMC chain are calculated for each one returned. This is useful if autocorrelation is high.
}
  \item{burnInSteps}{
number of steps to discard as burn-in at the beginning of the chain.
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
The function generates vectors of random draws from the posterior distributions of the population mean (\eqn{\mu}) and standard deviation (\eqn{\sigma}). The procedure uses a Bayesian MCMC process implemented in JAGS (Plummer 2003).

Priors for all parameters can be specified by including elements in the \code{priors} list:\cr
\eqn{\mu} : the population mean has a normal prior, with mean \code{muM} and standard deviation \code{muSD};\cr
\eqn{\sigma} : the population standard deviation has gamma priors, with mode \code{sigmaMode} and standard deviation \code{sigmaSD}.

The default priors correspond to \code{list(muM = mean(y), muSD = sd(y)*5, sigmaMode = sd(y), sigmaSD = sd(y)*5)}.

If \code{priors = NULL}, vague priors as described by Kruschke (2013) are used for all parameters.
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
Mike Meredith, adapted from code in the \code{BEST} package.
}

\examples{
# Generate a sample from a normal distribution, maybe the head-body length of a
#   carnivore in mm:
HB <- rnorm(10, 900, 0.15)
Bnormal(HB)  # with broad default priors
Bnormal(HB, priors=list(muM=1000, muSD=200)) # an appropriate prior
Bnormal(HB, priors=list(muM=1, muSD=0.2))    # a silly prior produces a warning.

}
