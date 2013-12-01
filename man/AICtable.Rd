\name{AICtable}
\alias{AICtable}
\title{
Make a table for AIC or other criterion
}
\description{
Takes the output from a call to AIC or AICc returns a data frame with model likelihoods and model weights.
}
\usage{
AICtable(x)
}
\arguments{
  \item{x}{
A data frame with the second column being AIC-type values, as returned by \code{AIC} or \code{AICc}.
}
}
\value{
Returns the original data frame with the following extra columns
  \item{Delta }{The difference between each value in x and min(x)}
  \item{ModelLik }{The model likelihood}
  \item{ModelWt }{The model weight}

It the original data frame had row names, the rows will be sorted by increasing values of AIC/AICc.
}
\references{
Burnham and Anderson (2002) \emph{Model selection and multimodel inference: a practical information-theoretic approach}, 2 edn. Springer-Verlag.
}
\author{
Mike Meredith
}

\examples{
data(KanhaTigers)
freq3 <- table(rowSums(KanhaTigers))
t3 <- ncol(KanhaTigers)
mods <- NULL
mods$M0 <- closedCapM0(freq3, t3)
mods$Mh2 <- closedCapMh2(freq3, t3)
mods$MhJK <- closedCapMhJK(freq3, t3)
mods$Mt <- closedCapMt(KanhaTigers)
AICc <- sapply(mods, AICc)
AICtable(AICc)
# MhJK does not use likelihood maximisation
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
