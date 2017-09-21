
# Bayesian modelling of normal distribution with Gibbs sampler
# ============================================================

# This version allows for a gamma prior for the precision, tau,
#  though the default is an improper flat prior.
# The gamma prior is specified by shape and rate.

Bnormal <- function(y, priors=NULL,
                    chains=3, sample=10000, burnin=100) {

  startTime <- Sys.time()

  # Data summaries
  if(!is.numeric(y))
    stop("'y' must be a numeric vector.")
  n <- length(y)
  if(n < 2)
    stop("'y' must contain at least 2 values")
  y.bar <- mean(y)

  m0 <- t0 <- a <- b <- 0
  if(!is.null(priors$muMean) && !is.null(priors$muSD)) {
    m0 <- priors$muMean
    t0 <- 1/(priors$muSD)^2
    if(y.bar > priors$muMean + priors$muSD || y.bar < priors$muMean - priors$muSD)
      warning("Sample mean is outside the prior range mMean \u00B1 sMean.")
  }
  if(!is.null(priors$tauShape) && !is.null(priors$tauRate)) {
    a <- priors$tauShape
    b <- priors$tauRate
  }
  aBit <- a + n / 2 # This doesn't change

  # Objects to hold results
  n.iter <- sample + burnin
  chainList <- vector('list', chains)
  chain <- matrix(nrow=n.iter, ncol=2) # will hold output
  colnames(chain) <- c("mu", "sigma")

  for(ch in 1:chains) {
    tau <- 1  # starting values
    for (t in 1:n.iter){
      # Draw mu from conjugate posterior with known sigma
      v <- 1 / (tau * n + t0)
      m <- v * (tau * sum(y) + t0 * m0)
      mu <-  rnorm(1, m, sqrt(v))
      # Draw tau from conjugate posterior with known mean
      tau <- rgamma(1, aBit, b + sum((y - mu)^2)/2)
      chain[t, ] <- c(mu, sqrt(1/tau))
    }
    chainList[[ch]] <- mcmc(chain[(burnin+1):n.iter, ])
  }
  MCMC <- mcmc.list(chainList)
  Rhat <- try(gelman.diag(MCMC, autoburnin=FALSE)$psrf[, 1], silent=TRUE)
  if(inherits(Rhat, "try-error") || !all(is.finite(Rhat)))
    Rhat <- NULL

  out <- as.Bwiqid(MCMC,
      header = "Model fitted in R with a Gibbs sampler",
      defaultPlot = names(MCMC)[1])
  attr(out, "call") <- match.call()
  attr(out, "n.chains") <- chains
  attr(out, "n.eff") <- effectiveSize(MCMC)
  attr(out, "Rhat") <- Rhat
  attr(out, "timetaken") <- Sys.time() - startTime
  return(out)
}
