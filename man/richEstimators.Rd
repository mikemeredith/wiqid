\name{Species richness estimators}
\alias{Species richness estimators}
\alias{richACE}
\alias{richICE}
\alias{richChao1}
\alias{richChao2}
\alias{richJack1}
\alias{richJack2}
\alias{richBoot}
\alias{richMM}

\title{
Species richness estimators
}
\description{
Functions to estimate species richness, based on samples from one or more surveys (quadrats, sites, occasions, ...) as included in EstimateS. See Details for individual estimators.
}
\usage{
richACE(cntVec, threshold = 10)
richICE(incMat, threshold = 10)
richChao1(cntVec, correct = FALSE)
richChao2(incMat, correct = FALSE)
richJack1(incMat)
richJack2(incMat)
richBoot(incMat)
richMM(incMat)
}
\arguments{
  \item{cntVec}{
a vector of species counts (abundances) with one element for each species. A matrix of counts, species x sites, may also be provided.
}
  \item{incMat}{
a 1/0 matrix of species incidence (presence), species x sites. A matrix of counts may also be provided.
}
  \item{threshold}{
the definition of rare or infrequent species: species with \code{threshold} or smaller counts (ACE) or incidences (ICE) are rare or infrequent.
}
  \item{correct}{
if TRUE, bias-corrected estimates are calculated.
}
}
\details{
\code{richACE} and \code{richICE} calculate Anne Chao's Abundance-based and Incidence-based Coverage Estimators of species richness respectively (Chao et al, 2000). 

\code{richChao1} and \code{richChao2} calculate Anne Chao's Chao 1 (abundance-based) and Chao 2 (incidence-based) estimators (Chao 1984, 1987). 

\code{richJack1} and \code{richJack2} calculate first and second order jackknife estimators of species richness (Smith & van Belle, 1984). 

\code{richBoot} calculates a bootstrap estimator of species richness (Smith & van Belle, 1984). 

\code{richMM} calculates the assymptotic species richness from a Michaelis-Menten curve fitted to the species rarefaction curve (Colwell et al. 2004).

See the EstimateS Users Guide at \url{http://viceroy.eeb.uconn.edu/estimates/EstimateSPages/EstSUsersGuide/EstimateSUsersGuide.htm} for the equations used.
}
\value{
\code{richChao1} and \code{richChao2} return a vector with a point estimate, upper and lower 95\% confidence limits, and standard deviation.

The other functions return a scalar.
}
\references{
Chao, A. 1984. Non-parametric estimation of the number of classes in a population. \emph{Scandinavian Journal of Statistics} 11, 265-270.

Chao, A. 1987. Estimating the population size for capture-recapture data with unequal capture probabilities. \emph{Biometrics} 43:783-791. 

Chao, A., W.-H. Hwang, Y.-C. Chen, and C.-Y. Kuo. 2000. Estimating the number of shared species in two communities. \emph{Statistica Sinica} 10:227-246.

Colwell, R. K., C. X. Mao, & J. Chang. 2004. Interpolating, extrapolating, and comparing incidence-based species accumulation curves. \emph{Ecology} 85, 2717-2727.

Smith, E.P. & van Belle, G. 1984. Nonparametric estimation of species richness. \emph{Biometrics} 40, 119-129.

}
\author{
Mike Meredith
}
\note{
%%  ~~further notes~~
}
\section{Benchmarks}{
Output from these functions has been checked against EstimateS for the \code{seedbank} and \code{killarneyBirds} data sets.
} 

\seealso{
\code{\link{richRarefy}} for rarefaction curves, and \code{\link{richCurve}} for a function to give richness estimates for sub-sets of samples.
}
\examples{
data(seedbank)

richACE(seedbank)

}
