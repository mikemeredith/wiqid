\name{Closed Captures}
\alias{closedCapM0}
\alias{closedCapMh2}
\alias{closedCapMhJK}

\title{
Analysis of mark-recapture data for closed populations
}
\description{
Functions to analyse the classical models for closed populations without covariates. Currently M0, Mh2 (a 2-mixture model), and the jackknife estimator for Mh are available.
}
\usage{
closedCapM0(freq, n.occ = length(freq))

closedCapMh2(freq, n.occ = length(freq))

closedCapMhJK(freq, n.occ = length(freq))
}
\arguments{
  \item{freq}{
a vector of capture frequencies; trailing zeros are not required provided n.occ is given.
}
  \item{n.occ}{
the total number of capture occasions.
}
}
\details{
Model M0 assumes all animals have the same capture probability.

Model Mh2 and the jackknife estimator allow for heterogeneity in capture probability.
}
\value{
A matrix with columns for the point estimate and 95\% confidence interval, and rows for each of the parameters estimated.

The matrix has an attribute with the AIC, though this will always be NA for the jackknife estimator.

}
\references{
Basic work on mark-recapture for closed populations is in:

Otis, D L; K P Burnham; G C White; D R Anderson. 1978. Statistical inference from capture data on closed animal populations. \emph{Wildlife Monographs} 62:1-135.

White, G C; D R Anderson; K P Burnham; D L Otis. 1982. \emph{Capture-recapture and removal methods for sampling closed populations}. Los Alamos National Laboratory, Los Alamos NM.

The jackknife estimator is described in:

Burnham, K P; W S Overton. 1979. Robust estimation of population size when capture probabilities vary among animals. \emph{Ecology} 60:927-936.

Rexstad, E; K Burnham 1992. \emph{User's guide for interactive program CAPTURE}. USGS Patuxent.

Rabbit data in the examples are from:

Edwards, W.R. and L.L. Eberhardt 1967. Estimating cottontail abundance from live trapping data. \emph{Journal of Wildlife Management} 31:87-96. 
}
\author{
Mike Meredith
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Data from White et al (1982):
freq1 <- c(50, 46, 35, 24, 14, 5, 0) # there were 7 capture occasions
closedCapM0(freq1)
closedCapMh2(freq1)
closedCapMhJK(freq1)

# Rabbit data from Edwards and Eberhart (1967)
freq2 <- c(43, 16, 8, 6, 0, 2, 1) ; t2 <- 18 
closedCapM0(freq2, t2)
closedCapMh2(freq2, t2)
closedCapMhJK(freq2, t2)

# Kanha tiger data
freq3 <- c(10,6,6,2,2) ; t3 <- 10
closedCapM0(freq3, t3)
closedCapMh2(freq3, t3)
closedCapMhJK(freq3, t3)

}
