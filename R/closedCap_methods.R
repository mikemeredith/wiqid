# Methods for class closedCap, ie. output for closedCap* family of functions.

print.closedCap <- function(x, digits=4, ...)  {
  cat("Call: ")
  print(x$call)
  cat("\nReal values:\n")
  print(x$real, digits=digits, ...)
  cat("\nAIC:", AIC(x), "\n")
}

logLik.closedCap <- function(object, ...)  {
  tmp <- as.vector(object$logLik)
  ll <- tmp[1]
  attr(ll, 'df') <- tmp[2]
  attr(ll, 'nobs') <- tmp[3]
  class(ll) <- "logLik"
  return(ll)
}

nobs.closedCap <- function(object, ...)  {
  nobs <- as.numeric(object$logLik)[3]
  return(if(is.null(nobs)) NA_integer_ else nobs)
}
