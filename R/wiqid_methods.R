# Print method for class wigid, ie. output for wiqid MLE functions.

print.wiqid <- function(x, digits=4, ...)  {
  cat("Call: ")
  print(x$call)
 if(anyDuplicated(x$real) > 0) {
    cat("\nReal values (duplicates omitted):\n")
    print(unique(x$real), digits=digits, ...)
  } else {
    cat("\nReal values:\n")
    print(x$real, digits=digits, ...)
  }
  cat("\nAIC:", AIC(x), "\n")
}

logLik.wiqid <- function(object, ...)  {
  tmp <- as.vector(object$logLik)
  ll <- tmp[1]
  attr(ll, 'df') <- tmp[2]
  attr(ll, 'nobs') <- tmp[3]
  class(ll) <- "logLik"
  return(ll)
}

nobs.wiqid <- function(object, ...)  {
  nobs <- as.numeric(object$logLik)[3]
  return(if(is.null(nobs)) NA_integer_ else nobs)
}
