\name{occSSrn}
\alias{occSSrn}

\title{
Royle-Nichols model for occupancy estimation
}
\description{
Inplements the Royle-Nichols method for estimation of site occupancy allowing for abundance-induced heterogeneity in detection probability.
}
\usage{
occSSrn(y, n)
}
\arguments{
  \item{y}{
a vector with the number of detections at each site.
}
  \item{n}{
a vector with the number of visits (survey occasions) at each site.
}
}

\value{

A matrix with columns for the point estimate and 95\% confidence interval and a row for each of the parameters estimated, with occupancy in the first row.

The matrix has an attribute with the AIC of the model.
}
\references{
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
See \code{link{occSS0}} for other single-season estimators.
}
\examples{
# The blue ridge salamanders data from MacKenzie et al (2006) p99:
y <- rep(0:4, c(21,12,1,4,1))
n <- rep(5, 39)
occSSrn(y, n)
}
