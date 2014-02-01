\name{weta}
\alias{weta}
\docType{data}
\title{
Detection data for weta in gorse bushes
}
\description{
Results of an occupancy survey to see if presence/absence of weta in 72 gorse bushes is affected by browsing of the bushes by goats. Probability of detection may be different for each observer, and observer ID is also recorded.
}
\usage{data(weta)}
\format{
  A data frame with 72 observations on the following 11 variables.
  \describe{
    \item{\code{D1, D2, D3, D4, D5}}{a numeric vector for each of 5 daily suveys, showing whether weta were detected (1) or not detected (0); NA if the bush was not surveyed on the relevant day. }
    \item{\code{Browsed}}{a logical vector indicating whether the bush was browsed (TRUE) or not browsed (FALSE)}
    \item{\code{ObsD1, ObsD2, ObsD3, ObsD4, ObsD5}}{a factor with levels \code{A} \code{B} \code{C}, indicating which observer carried out each survey.}
  }
}
\details{
The data are provided as a data frame, such as would result from reading in data for a text file. Further formatting is needed before using these for analysis: see the examples. 
}
\source{
Discussed in MacKenzie et al (2006) p116. Data distributed with PRESENCE.
}
\references{
MacKenzie, D I; J D Nichols; A J Royle; K H Pollock; L L Bailey; J E Hines 2006. \emph{Occupancy estimation and modeling : inferring patterns and dynamics of species occurrence}. Elsevier Publishing.
}
\examples{
data(weta)

DH <- weta[, 1:5]  # extract the detection history:
occSScov(DH)
occSStime(DH, p~.time)
# With covariates
occSScov(DH, list(psi ~ Browsed, p ~ ObsD), data=weta)
occSScov(DH, list(psi ~ Browsed, p ~ Browsed), data=weta)

}
\keyword{datasets}
