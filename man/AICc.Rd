\name{AICc}
\alias{AICc}
\title{
Akaike's Information Criterion with small-sample correction - AICc
}
\description{
Akaike's Information Criterion with small-sample correction - AICc
}
\usage{
AICc(object, ..., nobs)
}

\arguments{
  \item{object}{
a fitted model object for which there exists a logLik method to extract the corresponding log-likelihood, number of parameters, and number of observations.
}
  \item{\dots}{
optionally more fitted model objects.
}
  \item{nobs}{
the value to use for the effective sample size; overrides the value contained in the model object(s).
}
}
\details{
This is an S3 generic, with a default method which calls \code{\link{logLik}}, and should work with any class that has a \code{logLik} method.
}
\value{
If just one object is provided, the corresponding AICc.

If multiple objects are provided, a data frame with rows corresponding to the objects and columns representing the number of parameters in the model (df) and the AICc.

The result will be \code{Inf} for overparameterised models, ie. when \code{df >= nobs - 1}.
}
\references{
Burnham, K P; D R Anderson 2002. \emph{Model selection and multimodel inference: a practical information-theoretic approach}. Springer-Verlag. 
}
\author{
Essentially the same as \code{\link{AIC}} in package \code{stats}. Modified to return AICc by Mike Meredith. 
}
\note{
For some data types, including occupancy data, there is debate on the appropriate effective sample size to use.
}

\seealso{
\code{\link{AIC}}.
}
\examples{
data(salamanders)
mt <- occSS.t(salamanders)
mT <- occSS.T(salamanders, plot=FALSE)
AIC(mt, mT)
AICc(mt, mT)
nobs(mt)
AICc(mt, mT, nobs=10)
}