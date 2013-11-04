\name{cjs0}
\alias{cjs0}
\alias{cjsTime}
\alias{ch2mArray}
\title{
Survival from recapture data
}
\description{
Calculation of apparent survival (accounting for recapture probability) from recapture data formated as an m-array. \code{cjs0} runs a phi(.)p(.) model. \code{cjsTime} allows for time-dependent phi or p, possibly with covariates. \code{ch2mArray} converts capture history matrices to m-arrays.
}
\usage{
cjs0(mArray, ci = 0.95)

cjsTime(mArray, phi = ~1, p = ~1, data = NULL, ci = 0.95)

ch2mArray(CH)
}
\arguments{
  \item{mArray}{
matrix giving the recapture history in m-array format, with zeros (not NAs) in the lower triangle.
}
  \item{ci}{
the required confidence interval.
}
  \item{phi}{
a one-sided formula defining a linear predictor for apparent survival (phi).
}
  \item{p}{
a one-sided formula defining a linear predictor for recapture probability (p).
}
  \item{data}{
a data frame with a row for each survival interval / recapture occasion and columns for each of the covariates used to estimate phi or p.
}
  \item{CH}{
a 1/0 matrix with capture histories with a row for each animal captured and a column for each capture occasion.
}
}
\details{
For n REcapture occasions (ie. n + 1 capture occasions total), the m-array has:

  = n rows corresponding to marked animals caught and released on each occasion except for the last (those animals will not be recatured anyway).

  = n columns with the number of animals recaptured on each occasion except the first, PLUS a last column with the number never recaptured.

 The lower triangle of the matrix is all zeros.


}
\value{
\code{cjs0} returns a matrix with estimates of phi and p in rows and columns for the estimate and upper and lower confidence intervals, and AIC as an attribute.

\code{cjsTime} returns a list:
  \item{beta }{Estimates of the coefficients in the linear predictors for phi and p.}
  \item{real }{Back-transformed estimates of phi and p for each interval / occasion. }
  \item{AIC }{AIC for the model.}

\code{ch2mArray} returns an m-array.
}
\references{
Lebreton, J-D; K P Burnham; J Clobert; D R Anderson. 1992. Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. \emph{Ecological Monographs} 62:67-118.
}
\author{
Mike Meredith
}

\section{Benchmarks}{
Output has been checked against program MARK with the dipper data set: coefficients are not the same as MARK uses models without an intercept, but the real values agree to 3 decimal places.
}

\examples{
# The classic Lebreton et al dipper data set:
dippers <- matrix(c(
  11,2,0,0,0,0,9,
  0,24,1,0,0,0,35,
  0,0,34,2,0,0,42,
  0,0,0,45,1,2,32,
  0,0,0,0,51,0,37,
  0,0,0,0,0,52,46),nrow=6,byrow=TRUE)

cjs0(dippers)  # phi(.) p(.) model
cjsTime(dippers)  # default is also a phi(.) p(.) model
cjsTime(dippers, phi = ~time)  # using in-built time covariate
dd <- data.frame(flood = c(0, 1, 1, 0, 0, 0))
cjsTime(dippers, phi = ~flood, p = ~time, data=dd)

# The following model is overparameterised: phi[6] and p[6] cannot
#    be estimated separately and the variance-covariance matrix
#    cannot be calculated:
cjsTime(dippers, phi = ~time, p = ~time)

}
