\name{railSims}
\alias{railSims}
\docType{data}
\title{
Simulated detection/non-detection data for two species of rails
}
\description{
A data set for single-season two-species occupancy modeling. See \code{\link{occ2sps}} for details of these kinds of models.
}
\usage{data("railSims")}
\format{
  A data frame with detection (1) vs non-detection data for 2 species at 160 sites on three occasions.
  \describe{
    \item{A1, A2, A3}{detection histories for the dominant species on 3 occasions}
    \item{B1, B2, B3}{detection histories for the subordinate species on 3 occasions}
  }
}
\details{
The data come from a simulated scenario with the following parameters:
  \tabular{lll}{
  psiA  \tab = 0.39 = \tab probability of occupancy of species A \cr
  psiBA \tab = 0.77 = \tab probability of occupancy of B if A is present \cr
  psiBa \tab = 0.39 = \tab probability of occupancy of B if A is absent \cr
  pA \tab = 0.74 = \tab probability of detection of species A if B is absent \cr
  pB \tab = 0.80 = \tab probability of detection of species B if A is absent \cr
  rA \tab = 0.60 = \tab probability of detection of species A if both are present \cr
  rBA \tab = 0.88 = \tab probability of detection of species B if both are present and A was detected \cr
  rBa \tab = 0.79 = \tab probability of detection of species B if both are present and A was not detected}

These values approximate those reported by Richmond et al (2010).
}
\source{
Simulated data
}
\references{
Richmond, O.M.W., Hines, J.E., & Beissinger, S.R. (2010) Two-species occupancy models: a new parameterization applied to co-occurrence of secretive rails. Ecological Applications, 20, 2036-2046.
}
\examples{
data(railSims)
# Separate the two detection histories
DHA <- railSims[, 1:3]
DHB <- railSims[, 4:6]
# Default model (no interaction)
occ2sps(DHA, DHB)
# Model with full interaction
occ2sps(DHA, DHB, modelSpec=224)
}
\keyword{datasets}
