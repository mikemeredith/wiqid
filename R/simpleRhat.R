
# simplified Gelman convergence diagnostic using Brooks & Gelman's "interval" method.

# See Brooks & Gelman (1998) General methods for monitoring convergence of iterative simulations. J Computational and Graphical Statistics, 7, 434-455. p. 441

# This follows WinBUGS in using the central 80% interval as the measure of width (WinBUGS manual p.27).

simpleRhat <- function(object, n.chains, burnin=0) {

  width <- function(y)
    diff(quantile(y, c(0.1, 0.9)))

  x <- as.matrix(object)
  if(missing(n.chains) || !is.numeric(n.chains)) {
    n.chains <- attr(object, "n.chains")
    if(is.null(n.chains))
     stop("Please supply a value for 'n.chains'!", call.=FALSE)
  }
  if(n.chains < 2)
    stop("More than 1 chain is needed to calculate Rhat.", call.=FALSE)
  if(nrow(x) %% n.chains > 0)
    stop("The number of rows in 'x' must be a multiple of 'n.chains'.", call.=FALSE)
  if(burnin < 0 || burnin > 0.9)
    stop("'burnin' must be between 0 and 0.9.", call.=FALSE)

  n.par <- ncol(x)                     # number of parameters
  parNames <- colnames(x)
  n.iter <- nrow(x)                    # total sample size
  ipc <- n.iter / n.chains             # iterations per chain (T)

  dim(x) <- c(ipc, n.chains, n.par)    # separate the chains
  if(burnin) {
    discard <- round(ipc*burnin)       # iters to discard as burn-in
    x <- x[-(1:discard), , ]           # discard unwanted iters
    ipc <- ipc - discard
  }

  W0 <- apply(x, 2:3, width)           # width of individual chains
  W <- colMeans(W0)
  
  dim(x) <- c(ipc * n.chains, n.par)   # combine the chains (after burn-in)
  B <- apply(x, 2, width)              # width of pooled chains

  Rhat <- B / W
  names(Rhat) <- parNames
  return(Rhat)
}


# An error-catching wrapper for coda::effectiveSize
safeNeff <- function(x) {
  # x is a data frame or matrix with a column for each parameter
  safe1 <- function(v) {
    tmp <- try(coda::effectiveSize(v), silent=TRUE)
    if(inherits(tmp, "try-error"))
      return(NA)
    return(tmp)
  }
  apply(x, 2, safe1)
}

