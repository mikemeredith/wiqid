\name{Occupancy Single Season}
\alias{occSS0}
\alias{occSS.t}
\alias{occSS.T}
\alias{occSScovSite}
\alias{occSScov}
\alias{print.occupancy}

\title{
Single-season occupancy estimation
}
\description{
Functions to estimate occupancy from detection/non-detection data for a single season. 
}
\usage{
occSS0(y, n)

occSS.t(DH)

occSS.T(DH)

occSScovSite(y, n, psi = ~1, p = ~1, data = NULL)

occSScov(DH, psi = ~1, p = ~1, data = NULL)

\method{print}{occupancy}(x, digits = 4, ...)
}
\arguments{
  \item{y}{
a vector with the number of detections at each site.
}
  \item{n}{
a vector with the number of visits (survey occasions) at each site.
}
 \item{DH}{
a 1/0/NA matrix (or data frame) of detection histories, sites x occasions
}
  \item{psi}{
a one-sided formula symbolically defining a linear predictor for \code{psi} in terms of site covariates.
}
  \item{p}{
a one-sided formula symbolically defining a linear predictor for \code{p} in terms of site covariates and (for \code{occSScov}) survey covariates.
}
  \item{data}{
an optional data frame or list containing the variables in the model. If survey covariates are included this must be a list - see details.
}
  \item{x}{
an object of class \code{occupancy} to be printed.
}
  \item{digits}{
the number of digits to display.
}
  \item{...}{
other arguments passed to the print command
}
}
\details{
\code{occSS0} implements a simple model with one parameter for probability of occupancy and one for probability of detection, ie. a \code{psi(.) p(.)} model.

\code{occSS.t} estimates a separate detection probability for each survey occasion, ie. a \code{psi(.) p(t)} model.

\code{occSS.T} estimates a logit-linear trend in detection probability, producing estimates of the intercept and slope.

\code{occSScovSite} allows for psi or p to be modelled as a logistic function of site covariates, while \code{occSScov} extends this to allow for survey covariates to be included in the model for p. \code{occSScov} has a built in \code{.time} covariate which can be used for modelling p.

Survey covariates for use by \code{occSScov} should be included in the \code{data} list either as a sites x surveys matrix or as a vector of length sites*surveys, beginning with the observations for all the sites for the first survey. The latter is useful when the covariate is a factor, as matrices of factors don't generally work in R. Using \code{as.factor(covariate_matrix)} produces the right vector.

For speed, use the simplest function which will cope with your model. For example, you can run psi(.) p(.) models in \code{occSScovSite} or \code{occSScov}, but \code{occSS0} is much faster.
}
\value{
Returns an object of class \code{occupancy}, which is a list with the following elements:

  \item{call}{The call used to produce the results}
  \item{beta}{Values of the coefficients of the terms in the linear predictors, with standard errors and 95\% conidence intervals.}
  \item{real}{Estimates of occupancy and probability of detection on the real scale, with 95\% confidence intervals.}
  \item{AIC}{The AIC of the fitted model.}
}
\section{Benchmarks}{Output has been checked against output from PRESENCE (Hines 2006) v.5.5 for the \code{\link{salamanders}} and \code{\link{weta}} data sets. Real values are mostly the same to 4 decimal places, though there is occasionally a discrepancy of 0.0001. AICs are the same.
}
\references{
MacKenzie, D I; J D Nichols; G B Lachman; S Droege; J A Royle; C A Langtimm. 2002. Estimating site occupancy rates when detection probabilities are less than one. \emph{Ecology} 83:2248-2255.

MacKenzie, D I; J D Nichols; A J Royle; K H Pollock; L L Bailey; J E Hines 2006. \emph{Occupancy estimation and modeling : inferring patterns and dynamics of species occurrence}. Elsevier Publishing.

Hines, J. E. (2006). PRESENCE - Software to estimate patch occupancy and related parameters. SGS-PWRC. \url{http://www.mbr-pwrc.usgs.gov/software/presence.html}.
}
\author{
Mike Meredith
}
\note{
%%  ~~further notes~~
}
\seealso{
\code{\link{occSSrn}} for the Royle-Nichols model. See the examples for the \code{\link{weta}} data set.
}
\examples{
# The blue ridge salamanders data from MacKenzie et al (2006) p99:
data(salamanders)
y <- rowSums(salamanders)
n <- rowSums(!is.na(salamanders))

occSS0(y, n)
occSS.t(salamanders)
occSS.T(salamanders)
}
