
# Bayesian modelling of normal distribution with Gibbs sampler
# ============================================================

Bnormal1 <- function(x, priors=list(),
                    n.chains=3, n.iter=10100, n.burnin=100) {
                    
  startTime <- Sys.time()
  
  # Data summaries
  n <- length(x)
  stopifnot(n > 1)
  x.bar <- mean(x)

  if(!is.null(priors)) {
    priorsDefault <- list(mMean=0, sMean=100, aSigma=0.001, bSigma=0.001)
    priors <- replace (priorsDefault, names(priors), priors)
    if(x.bar > priors$mMean + priors$sMean || x.bar < priors$mMean - priors$sMean)
      warning("Sample mean is outside the prior range mMean \u00B1 sMean.")
    m0 <- priors$mMean
    t0 <- 1/(priors$sMean)^2
    a <- priors$aSigma
    b <- priors$bSigma
  } else {
    m0 <- t0 <- a <- b <- 0
  }
  aBit <- a + n / 2 # This doesn't change

 
  # Starting values
  # tauStart <- sd(x) * 2^runif(n.chains, -1, 1)

  # Objects to hold results
  chainList <- vector('list', n.chains)
  chain <- matrix(nrow=n.iter, ncol=2) # will hold output
  colnames(chain) <- c("mu", "sigma")

  for(ch in 1:n.chains) {
    tau <- 1  # starting values
    for (t in 1:n.iter){  
      # Draw mu from conjugate posterior with known sigma
      v <- 1 / (tau * n + t0)
      m <- v * (tau * sum(x) + t0 * m0)
      mu <-  rnorm(1, m, sqrt(v))
      # Draw tau from conjugate posterior with known mean
      tau <- rgamma(1, aBit, b + sum((x - mu)^2)/2) 
      chain[t, ] <- c(mu, sqrt(1/tau))
    } 
    chainList[[ch]] <- mcmc(chain[(n.burnin+1):n.iter, ])
  }
  MCMC <- mcmc.list(chainList)
  Rhat <- try(gelman.diag(MCMC, autoburnin=FALSE)$psrf[, 1], silent=TRUE)
  if(inherits(Rhat, "try-error") || !all(is.finite(Rhat)))
    Rhat <- NULL
  
  out <- as.Bwiqid(MCMC,
      header = "Model fitted in R with a Gibbs sampler",
      defaultPlot = names(MCMC)[1])
  attr(out, "call") <- match.call()
  attr(out, "n.chains") <- n.chains
  attr(out, "n.eff") <- effectiveSize(MCMC)
  attr(out, "Rhat") <- Rhat
  attr(out, "timetaken") <- Sys.time() - startTime
  return(out)
}
