\name{Occupancy Single Season}
\alias{occSS0}
\alias{occSSrn}
\alias{occSStime}
\alias{occSScovSite}
\alias{occSScov}
\alias{occSSautocorA}

\title{
Single-season occupancy estimation
}
\description{
Functions to estimate occupancy from detection/non-detection data for a single season. 
}
\usage{
occSS0(y, n, ci=0.95)

occSSrn(y, n, ci=0.95)

occSStime(DH, model=p~1, data=NULL, ci=0.95, plot=TRUE)

occSSautocorA(DH, ci=0.95)

occSScovSite(y, n, model=list(psi~1, p~1), data = NULL, ci=0.95)

occSScov(DH, model=list(psi~1, p~1), data = NULL, ci=0.95)
}
\arguments{
  \item{y}{
a vector with the number of detections at each site.
}
  \item{n}{
a scalar or vector with the number of visits (survey occasions) at each site.
}
  \item{ci}{
the confidence interval to use.
}
 \item{DH}{
a 1/0/NA matrix (or data frame) of detection histories, sites x occasions
}
  \item{model}{
a list of formulae symbolically defining a linear predictor for each parameter in terms of covariates.
}
  \item{data}{
a data frame containing the variables in the model. For \code{occSStime}, a data frame with a row for each survey occasion; otherwise, a row for each site. Each site covariate has one column. Each survey covariate has one column for each occasion, and the column name must end with the occasion number (without leading zeros); eg, \code{Cov1, Cov2, ..., Cov15}. All covariates should be included in \code{data}, otherwise they will be sought in enclosing environments, which may not produce what you want -- and they won't be standardised.
}
  \item{plot}{
if TRUE (default), draws a plot of probability of detection vs time.
}
}
\details{
\code{occSS0} implements a simple model with one parameter for probability of occupancy and one for probability of detection, ie. a \code{psi(.) p(.)} model.

\code{occSSrn} implements the Royle-Nichols method for estimation of site occupancy allowing for abundance-induced heterogeneity in detection probability.

\code{occSStime} allows for a time-varying covariate that is the same across all sites, eg, moon-phase. A categorical time variable \code{.time} and a time trend \code{.Time} are built-in. A plot of detection probability vs time is produced if \code{plot=TRUE}.

\code{occSSautocorA} allows for autocorrelation among surveys: the probability of detection on a survey depends on detection/nondetection on the previous survey (Hines et al 2010, Appendix A).

\code{occSScovSite} allows for psi or p to be modelled as a logistic function of site covariates, while \code{occSScov} extends this to allow for survey covariates to be included in the model for p. \code{occSScov} has a built in \code{.time} covariate which can be used for modelling p with time as a fixed effect, and \code{.Time} for a linear or quadratic trend.

Numeric covariates in \code{data} are standardised to facilitate convergence. This applies to binary covariates coded as 1/0; if this is not what you want, code these as TRUE/FALSE or as factors. 

For speed, use the simplest function which will cope with your model. For example, you can run psi(.) p(.) models in \code{occSScovSite} or \code{occSScov}, but \code{occSS0} is much faster.
}
\value{
Returns an object of class \code{wiqid}, which is a list with the following elements:

  \item{call}{The call used to produce the results}
  \item{beta}{Values of the coefficients of the terms in the linear predictors, with standard errors and confidence intervals.}
  \item{beta.vcv }{The variance-covariance matrix for the beta estimates.}
  \item{real}{Estimates of occupancy and probability of detection on the real scale, with confidence intervals; \code{occSSrn} gives occupancy, lambda and r.}
  \item{logLik}{a vector with elements for log(likelihood), number of parameters, and effective sample size. If parameters \emph{and their SEs} cannot be estimated, the first element should be \code{NA}.}

There are \code{print}, \code{logLik}, and \code{nobs} methods for class \code{wiqid}.
}

\section{Benchmarks}{Output has been checked against output from PRESENCE (Hines 2006) v.5.5 for the \code{\link{salamanders}} and \code{\link{weta}} data sets. Real values are mostly the same to 4 decimal places, though there is occasionally a discrepancy of 0.0001. AICs are the same.
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

occSS0(y, n)
occSSrn(y, n)
occSStime(salamanders, p ~ .time)  # time as a fixed effect
occSStime(salamanders, p ~ .Time + I(.Time^2))  # a quadratic time effect
}
