\name{weta}
\alias{weta}
\docType{data}
\title{
Detection data for weta in gorse bushes
}
\description{
Rewults of an occupancy survey to see if presence/absence of weta in 72 gorse bushes is affected by browsing of the bushes by goats. Probability of detection may be different for each observer, and observer ID is also recorded.
}
\usage{data(weta)}
\format{
  A data frame with 72 observations on the following 11 variables.
  \describe{
    \item{\code{D1 - D5}}{a numeric vector for each of 5 daily suveys, showing whether weta were detected (1) or not detected (0); NA if the bush was not surveyed on the relevant day. }
    \item{\code{Browsed}}{a numeric vector indicating whether the bush was browsed (1) or not browsed (0)}
    \item{\code{ObsD1-obsD5}}{a factor with levels \code{A} \code{B} \code{C}, indicating which observer carried out each survey.}
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
# Extract the detection history and summarize:
DH <- weta[, 1:5]
y <- rowSums(DH, na.rm=TRUE)
n <- rowSums(!is.na(DH))
weta1 <- occSS0(y, n)
weta2 <- occSS.t(DH)

# Extract site covariate and use this:
weta.covs <- weta[, 6, drop=FALSE]
  # drop=FALSE ensures that it remains as a data frame with 1 column.
weta3 <- occSScovSite(DH, psi= ~ Browsed, data=weta.covs)

# Extract the obs covariate and create the list of covariates
weta.obs <- weta[, 7:11]
weta.dat <- as.list(weta.covs)
weta.dat$obs <- as.matrix(weta.obs)
weta4 <- occSScov(DH, psi= ~ Browsed, p= ~ obs, data=weta.dat)

}
\keyword{datasets}
