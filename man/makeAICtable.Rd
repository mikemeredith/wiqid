\name{makeAICtable}
\alias{makeAICtable}
\title{
Make a table for AIC or other criterion
}
\description{
Takes a (named) vector of values for model AIC's or other criteria (AICc, QAIC, BIC, etc) and returns a data frame with model likelihoods and model weights.
}
\usage{
makeAICtable(x)
}
\arguments{
  \item{x}{
A vector of AIC-type values
}
}
\value{
Returns a data frame with the following columns
  \item{x }{The original vector; the column name will be derived from the name of the object provided.}
  \item{Delta }{The difference between each value in x and min(x)}
  \item{ModelLik }{The model likelihood}
  \item{ModelWt }{The model weight}

It the original vector was named, the names will be used as the row names of the data frame and the rows will be sorted by increasing values of x.
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
makeAICtable(AICc)
# MhJK does not use likelihood maximisation
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
