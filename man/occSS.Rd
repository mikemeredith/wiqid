\name{Occupancy Single Season}
\alias{occSS0}
\alias{occSS.t}
\alias{occSS.T}
\alias{occSSrn}
\alias{occSScovSite}
\alias{occSScov}

\title{
Single-season occupancy estimation
}
\description{
Functions to estimate occupancy from detection/non-detection data for a single season. 
}
\usage{
occSS0(y, n)

occSS.t(DH, time = TRUE)

occSS.T(DH)

occSScovSite(DH, psi = ~1, p = ~1, data = NULL)

occSScov(DH, psi = ~1, p = ~1, data = NULL)

occSSrn(y, n)
}
\arguments{
  \item{y}{
a vector with the number of detections at each site.
}
  \item{n}{
a vector with the number of visits (survey occasions) at each site.
}
 \item{DH}{
a 1/0/NA matrix of detection histories, sites x occasions
}
  \item{time}{
logical; if FALSE, only one value for the probability of detection is estimated.
}
  \item{psi}{
a one-sided formula symbolically defining a linear predictor for \code{psi} in terms of site covariates.
}
  \item{p}{
a one-sided formula symbolically defining a linear predictor for \code{p} in terms of site covariates.
}
  \item{data}{
an optional data frame or list containing the variables in the model. If survey covariates are included this must be a list, with each survey covariate as a sites x surveys matrix.
}
}
\details{
\code{occSS0} implements a simple model with one parameter for probability of occupancy and one for probability of detection, ie. a \code{psi(.) p(.)} model.

\code{occSS.t} with \code{time = TRUE} estimates a separate detection probability for each survey occasion, ie. a \code{psi(.) p(t)} model. With \code{time = FALSE}, it is equivalent to occSS0.

\code{occSS.T} estimates a logit-linear trend in detection probability, producing estimates of the intercept and slope.

\code{occSScovSite} allows for psi or p to be modelled as a logistic function of site covariates, while \code{occSScov} extends this to allow for survey covariates to be included in the model for p. \code{occSScov} has a built in \code{.time} covariate which can be used for modelling p.

occSSrn implements the Royle-Nichols model for abundance-induced heterogeneity in detection probability.
}
\value{

##TODO: Change return values so that all return a list.

A matrix with columns for the point estimate and 95\% confidence interval and a row for each of the parameters estimated, with occupancy in the first row.

The matrix has an attribute with the AIC of the model.
}
\references{
MacKenzie, D I; J D Nichols; G B Lachman; S Droege; J A Royle; C A Langtimm. 2002. Estimating site occupancy rates when detection probabilities are less than one. \emph{Ecology} 83:2248-2255.

MacKenzie, D I; J D Nichols; A J Royle; K H Pollock; L L Bailey; J E Hines 2006. \emph{Occupancy estimation and modeling : inferring patterns and dynamics of species occurrence}. Elsevier Publishing.

Royle, J A; J D Nichols. 2003. Estimating abundance from repeated presence-absence data or point counts. \emph{Ecology} 84:777-790.
}
\author{
Mike Meredith
}
\note{
%%  ~~further notes~~
}
\seealso{
See the examples for the \code{\link{weta}} data set.
}
\examples{
# The blue ridge salamanders data from MacKenzie et al (2006) p99:
BRS <- matrix(0, 39, 5)
BRS[c(4,10,13,14,41,42,43,82,83,84,85,86,87,89,90,92,118,121,128,
        129,130,131,153,154,157,167,168,169,194,195)] <- 1
y <- rowSums(BRS)
n <- rep(5, 39)

occSS0(y, n)
occSS.t(BRS, time = FALSE)
occSS.t(BRS, time = TRUE)
occSS.T(BRS)
occSSrn(y, n)

}
