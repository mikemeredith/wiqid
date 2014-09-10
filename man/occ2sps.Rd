\name{occ2sps}
\alias{occ2sps}

\title{
Single-season two-species occupancy estimation
}

\description{
Estimates occupancy and probability of detection for two species, where one (dominant) species affects the occupancy or detection of the other (subordinate) species (see Richmond et al, 2010). The model has the following parameters:
  \tabular{ll}{
  psiA \tab probability of occupancy of species A \cr
  psiBA \tab probability of occupancy of B if A is present \cr
  psiBa \tab probability of occupancy of B if A is absent \cr
  pA \tab probability of detection of species A if B is absent \cr
  pB \tab probability of detection of species B if A is absent \cr
  rA \tab probability of detection of species A if both are present \cr
  rBA \tab probability of detection of species B if both are present and A was detected \cr
  rBa \tab probability of detection of species B if both are present and A was not detected}
}
\usage{
occ2sps(DHA, DHB, modelSpec = 111, ci = 0.95)
}
\arguments{
  \item{DHA}{
a 1/0/NA matrix (or data frame) of detection histories, sites x occasions, for the dominant species.
}
  \item{DHB}{
detection histories for the subordinate species in the same format as DHA.
}
  \item{modelSpec}{
a 3-digit number, where the first digit specifies the occupancy model (1: psiBA = psiBa, 2: psiBA != psiBa); the second, detection probability of species A (1: pA = rA, 2: pA != rA); and the third, detection probability of species B (1: pB = rBA = rBa,
  2: pB != rBA = rBa, 3: pB = rBa != rBA, 4: pB != rBA != rBa).
}
  \item{ci}{
the confidence interval to use.
}
}

\value{
Returns an object of class \code{wiqid}, which is a list with the following elements:

  \item{call}{The call used to produce the results}
  \item{beta}{Estimates of the parameters on the logit scale, with standard errors and confidence intervals.}
  \item{beta.vcv }{The variance-covariance matrix for the beta estimates.}
  \item{real}{Estimates of the parameters on the real scale, with confidence intervals}
  \item{logLik}{a vector with elements for log(likelihood), number of parameters, and effective sample size. If parameters \emph{and their SEs} cannot be estimated, the first element should be \code{NA}.}
}

\references{

Richmond, Hines, Bessinger (2010) Two-species occupancy models: a new parameterization applied to co-occurrence of secretive rails. \emph{Ecological Applications} 20(7):2036-2046

MacKenzie, D I; J D Nichols; A J Royle; K H Pollock; L L Bailey; J E Hines 2006. \emph{Occupancy estimation and modeling : inferring patterns and dynamics of species occurrence}. Elsevier Publishing.
}
\author{
Mike Meredith
}

\seealso{
See \code{\link{occSS}} for single-season single-species occupancy estimation.
}
\examples{
## To do
}

