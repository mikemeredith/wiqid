# Print method for class occupancy, ie. output for occSSxxx family of functions.

print.occupancy <- function(x, digits=4, ...)  {
  cat("Call: ")
  print(x$call)
  cat("\nReal values (duplicates omitted):\n")
  print(unique(x$real), digits=digits, ...)
  cat("\nAIC:", AIC(x), "\n")
}

logLik.occupancy <- function(object, ...)  {
  tmp <- as.vector(object$logLik)
  ll <- tmp[1]
  attr(ll, 'df') <- tmp[2]
  attr(ll, 'nobs') <- tmp[3]
  class(ll) <- "logLik"
  return(ll)
}

nobs.occupancy <- function(object, ...)  {
  nobs <- as.numeric(object$logLik)[3]
  return(if(is.null(nobs)) NA_integer_ else nobs)
}
