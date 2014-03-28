\name{window.Bwiqid}
\alias{window.Bwiqid}
\title{
Time windows for 'Bwiqid' objects
}
\description{
\code{window.Bwiqid} is a method for \code{Bwiqid} objects for the \code{window} generic function.

It extracts the subset of the observations between \code{start} and \code{end}. Setting \code{thin = k} selects every \code{k}th observation starting from \code{start}.
}
\usage{
\method{window}{Bwiqid}(x, start=NULL, end=NULL, thin=1, ...)
    
}
\arguments{
  \item{x}{
an object of class \code{Bwiqid}.
}
  \item{start}{
the first observation to retain.
}
  \item{end}{
the last observation to retain.
}
  \item{thin}{
the interval between retained observations.
}
  \item{\dots}{
further arguments for the window function.
}
}

\value{
Returns an object of class \code{Bwiqid} with the subset of observations.

Convergence diagnostics included as attributes in \code{x} are not retained, as they are unlikely to be valid for the subset. 
}
\author{
Mike Meredith.
}

\examples{
# TO DO
}
\keyword{methods}
