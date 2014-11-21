\name{Bayesian Occupancy Single Season}
\alias{BoccSS0}
\alias{BoccSS0_JAGS}


\title{
Single-season occupancy estimation using MCMC
}
\description{
Functions to estimate occupancy from detection/non-detection data for a single season using a Gibbs sampler coded in R or JAGS. 
}
\usage{
BoccSS0(y, n, psiPrior=c(1,1), pPrior=c(1,1), 
                    n.chains=3, n.iter=10100, n.burnin=100)
 
BoccSS0_JAGS(y, n, priorOnly=FALSE, ...)

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
  \item{n.chains}{
number of MCMC chains to run.
}
  \item{n.iter}{
number of values per chain, including burn-in.
}
  \item{n.burnin}{
number of values per chain to discard as burn-in.
}
  \item{priorOnly}{
if TRUE, the function produces random draws from the appropriate \emph{prior} distributions, with a warning.
}
  \item{...}{
further arguments passed to \code{\link{autorun.jags}}.
}
}
\details{
\code{BoccSS0} and \code{BoccSS0_JAGS} implement a simple model with one parameter for probability of occupancy and one for probability of detection, ie. a \code{psi(.) p(.)} model. The first uses a Gibbs sampler encoded in R, the second calls JAGS.

}
\value{
Returns an object of class \code{Bwiqid}, which is a data frame with a column for each parameter in the model.


There are \code{print} and \code{plot} methods for class \code{Bwiqid}.
}

\references{
MacKenzie, D I; J D Nichols; G B Lachman; S Droege; J A Royle; C A Langtimm. 2002. Estimating site occupancy rates when detection probabilities are less than one. \emph{Ecology} 83:2248-2255.

MacKenzie, D I; J D Nichols; A J Royle; K H Pollock; L L Bailey; J E Hines 2006. \emph{Occupancy estimation and modeling : inferring patterns and dynamics of species occurrence}. Elsevier Publishing.

Hines, J. E. (2006). PRESENCE - Software to estimate patch occupancy and related parameters. SGS-PWRC. http://www.mbr-pwrc.usgs.gov/software/presence.html.

Hines, J.E., Nichols, J.D., Royle, J.A., MacKenzie, D.I., Gopalaswamy, A.M., Kumar, N.S., Karanth, K.U. (2010) Tigers on trails: occupancy modeling for cluster sampling. \emph{Ecological Applications}, 20, 1456-1466.

}

\author{
Mike Meredith
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
