\name{Bayesian Occupancy Single Season}
\alias{BoccSS}
\alias{BoccSS0}


\title{
Bayesian single-season occupancy modelling
}
\description{
Functions to estimate occupancy from detection/non-detection data for a single season using a Gibbs sampler coded in R or JAGS.

\code{BoccSS0} runs a model in R without covariates, and allows priors to be specified as beta distributions for probability of occupancy and probability of detection.

\code{BoccSS} runs a model in R allowing for covariates, using a probit link and conjugate normal priors, which can be specified as mean and covariance.
}
\usage{
BoccSS0(y, n, psiPrior=c(1,1), pPrior=c(1,1),
                    n.chains=3, n.iter=10100, n.burnin=100)

BoccSS(DH, model=NULL, data=NULL, priors=list(),
                    n.chains=3, n.iter=11000, n.burnin=1000, parallel,
                    seed=NULL, doWAIC=FALSE)
}
\arguments{
  \item{y}{
a vector with the number of detections at each site.
}
  \item{n}{
a scalar or vector with the number of visits (survey occasions) at each site.
}
  \item{psiPrior, pPrior}{
parameters for beta distributions to be used as priors for psi and p.
}
  \item{DH}{
a 1/0/NA matrix (or data frame) of detection histories, sites x occasions.
}
  \item{model}{
a list of formulae symbolically defining a linear predictor for each parameter in terms of covariates. If NULL, an intercept-only model is used, ie, psi(.) p(.).
}
  \item{data}{
a data frame containing the variables in the model. For \code{occSStime}, a data frame with a row for each survey occasion; otherwise, a row for each site. Each site covariate has one column. Each survey covariate has one column for each occasion, and the column name must end with the occasion number (without leading zeros); eg, \code{Cov1, Cov2, ..., Cov15}. All covariates should be included in \code{data}, otherwise they will be sought in enclosing environments, which may not produce what you want -- and they won't be standardised.
}
  \item{priors}{
a list with elements for prior mean and variance for coefficients; see Details. If NULL, improper flat priors are used.
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
  \item{parallel}{
logical; if TRUE \emph{and} \code{n.chains} > 1 \emph{and} available cores (as returned by \code{parallel::detectCores}) > 2, chains will be run in parallel. If missing, chains will be run in parallel if \code{n.chains} < available cores.
}
  \item{doWAIC}{
logical; if TRUE, the Watanabe-Akaike Information Criterion is calculated. NOTE: THIS FEATURE IS STILL EXPERIMENTAL.
}
  \item{seed}{
for reproducible results; note that parallel and sequential methods use different random number generators, so will give different results with the same seed.
}

}
\details{
\code{BoccSS0} implements a simple model with one parameter for probability of occupancy and one for probability of detection, ie. a \code{psi(.) p(.)} model, using a Gibbs sampler implemented in R.

Independent beta distributions are used as priors for \code{BoccSS0}, as specified by \code{psiPrior} and \code{pPrior}. The defaults, \code{c(1, 1)}, correspond to uniform priors on the probabilities of occupancy and detection.

\code{BoccSS} uses a probit link to model occupancy and detection as a function of covariates (Dorazio and Rodriguez 2011); most software uses a logistic (logit) link.
See \link{Links}.
Coefficients on the probit scale are about half the size of the equivalent on the logit scale.

Priors for \code{BoccSS} are listed in the \code{priors} argument, which may contain elements:

\code{muPsi} and \code{muP} : the means for occupancy and detection coefficients respectively. This may be a vector with one value for each coefficient, including the intercept, or a scalar, which will be used for all. The default is 0.

\code{sigmaPsi} and \code{sigmaP} : the (co)variance for occupancy and detection coefficients respectively. This may be (1) a vector with one value for each coefficient, including the intercept, which represents the variance, assuming independence, or (2) a scalar, which will be used for all, or (3) a variance-covariance matrix. The default is 100.


When specifying priors, note that numerical covariates are standardized internally before fitting the model. For an intercept-only model, a prior of Normal(0, 1) on the probit scale implies a Uniform(0, 1) or Beta(1, 1) prior on the probability scale.

If you are unsure of the order of predictors, do a short run and check the output, or pass unusable values (eg, \code{muPsi=numeric(100)}) and check the error message.
}

\value{
Returns an object of class \code{Bwiqid}, which is a data frame with a column for each parameter in the model.


There are \code{print} and \code{plot} methods for class \code{Bwiqid}, as well as \link{diagnostic plots}.
}

\references{

MacKenzie, D I; J D Nichols; A J Royle; K H Pollock; L L Bailey; J E Hines 2006. \emph{Occupancy estimation and modeling : inferring patterns and dynamics of species occurrence}. Elsevier Publishing.

Dorazio and Rodriguez. 2012. A Gibbs sampler for Bayesian analysis of site-occupancy data. \emph{Methods in Ecology and Evolution}, 3, 1093-1098
}

\author{
Mike Meredith. \code{BoccSS} uses the Gibbs sampler described by Dorazio and Rodriguez (2012).
}
\seealso{
See the examples for the \code{\link{weta}} data set.
}
\examples{
# The blue ridge salamanders data from MacKenzie et al (2006) p99:
data(salamanders)
y <- rowSums(salamanders)
n <- rowSums(!is.na(salamanders))

tmp <- BoccSS0(y, n)
tmp
occSS0(y, n)  # for comparison
plot(tmp)


}
