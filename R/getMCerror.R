
# Calculation of MCMC error

# See Lunn et al 2013, The BUGS Book, p77.

# Lunn et al want number of batches = batch size = sqrt(chain length), but only
#   apply that to a single chain; they comment that with multiple chains, number of
#   batches will be larger.
# If there are many chains (eg, 20), this means number of batches >> batch size.
# To achieve number of batches = batch size, we use sqrt(total iterations),
#   corrected to be a multiple of number of chains.

# It gives the same result as coda::batchSE with an appropriate choice of batchSize.

getMCerror <- function(object, n.chains, SDpc=FALSE) {
  x <- as.matrix(object)                    # make sure it's a matrix
  # oldx <- x
  if(missing(n.chains) || !is.numeric(n.chains)) {
    n.chains <- attr(object, "n.chains")
    if(is.null(n.chains))
     stop("Please supply a value for 'n.chains'!", call.=FALSE)
  }
  if(nrow(x) %% n.chains > 0)
    stop("The number of rows in 'x' must be a multiple of 'n.chains'.", call.=FALSE)
  n.par <- ncol(x)                     # number of parameters
  parNames <- colnames(x)
  n.iter <- nrow(x)                    # total number of draws
  ipc <- n.iter / n.chains             # iterations per chain (T in Lunn et al)

  bpc <- sqrt(n.iter) %/% n.chains     # batches per chain (Q)
  bsize <- floor(nrow(x)/n.chains/bpc) # batch size (a)
  ni <- bpc * bsize                    # number of iters per chain to use

  dim(x) <- c(ipc, n.chains, n.par)        # separate the chains
  x <- x[1:ni, , ]                         # discard unwanted iters
  dim(x) <- c(bsize, bpc, n.chains, n.par) # separate the batches
  bm <- apply(x, 2:4, mean)                # get batch means
  dim(bm) <- c(bpc*n.chains, n.par)        # Combine bm's across chains
  sqdev <- (sweep(bm, 2, colMeans(bm), "-"))^2  # Get squared deviation
  SD <- sqrt(colSums(sqdev) * bsize / (n.chains*bpc - 1)) # sqrt(rho)
  MCE <- SD / sqrt(n.iter)
  names(MCE) <- parNames
  if(SDpc) {
    return(100 * MCE / apply(object, 2, sd))
  } else {
    return(MCE)
  }
}
