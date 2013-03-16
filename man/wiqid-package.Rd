\name{wiqid-package}
\alias{wiqid-package}
\alias{wiqid}
\docType{package}
\title{
Estimation functions for wildlife simulations
}
\description{
Quick and dirty functions to estimate occupancy, survival, abundance, etc. for wildlife populations during simulations and bootstraps.
}
\details{
\tabular{ll}{
Package: \tab wiqid\cr
Type: \tab Package\cr
Version: \tab 0.02\cr
Date: \tab 2013-03-16\cr
License: \tab GPL-3\cr
}
There are a number of sophisticated programs for the analysis of wildlife data, producing estimates of occupancy, survival, abundance, or density. But if you are generating hundreds of simulated data sets in R, you need fast, bare-bones functions which can be run from R.

So far the package has functions for single-season occupancy (cf PRESENCE), closed-population mark-recapture studies and CJS survival models from m-arrays (cf MARK), and species richness from species x sample matrices (cf EstimateS).

Simulations and bootstraps often generate weird data sets, eg. capture histories with no captures. These functions do not throw errors or give warnings if the data are weird, but return NAs if estimates cannot be calculated. Errors may still occur if the data are impossible, eg. 6 detections in 5 occasions.

}
\author{
Mike Meredith

Maintainer: Mike Meredith <mmeredith@wcs.org>
}
\references{
Williams, B K; J D Nichols; M J Conroy 2002. \emph{Analysis and management of animal populations: modeling, estimation, and decision making.} Academic Press, San Diego CA.
}
\keyword{ package }
\seealso{

}
