
# simplified Gelman convergence diagnostic using Brooks & Gelman's "interval" method.

# See Brooks & Gelman (1998) General methods for monitoring convergence of iterative simulations. J Computational and Graphical Statistics, 7, 434-455. p. 441

# This follows WinBUGS in using the central 80% interval as the measure of width (WinBUGS manual p.27).


# Calculate Rhat for a vector; called by simpleRhat; not exported

simpleRhat1 <- function(x, n.chains, burnin=0) {
  if(!is.numeric(x))
    return(NA_real_)
  width <- function(y)
    diff(quantile(y, c(0.1, 0.9)))
  xmat <- matrix(x, ncol=n.chains)
  if(burnin) {
    discard <- round(nrow(xmat)*burnin)
    xmat <- xmat[-(1:discard), ]
  }
  unname(width(xmat) / mean(apply(xmat, 2, width)))
}
#.......................................................................

#
simpleRhat <- function(x, n.chains, burnin=0) {
  if(missing(n.chains) || !is.numeric(n.chains))
    stop("Please supply a value for 'n.chains'!", call.=FALSE)
  if(n.chains < 2)
    stop("More than 1 chain is needed to calculate Rhat.", call.=FALSE)
  if(burnin < 0 || burnin > 0.9)
    stop("'burnin' must be between 0 and 0.9.", call.=FALSE)
  if(is.matrix(x) || is.data.frame(x)) {
    if(nrow(x) %% n.chains > 0)
      stop("The number of rows in 'x' must be a multiple of 'n.chains'.", call.=FALSE)
    return(apply(x, 2, simpleRhat1, n.chains=n.chains, burnin=burnin))
  } else if(is.vector(x)) {
    if(length(x) %% n.chains > 0)
      stop("The length of 'x' must be a multiple of 'n.chains'.", call.=FALSE)
    return(simpleRhat1(x, n.chains=n.chains, burnin=burnin))
  } else {
    stop("'x' must be a numeric vector, matrix or data frame.", call.=FALSE)
  }
}

