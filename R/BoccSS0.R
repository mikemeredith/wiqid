
# Bayes version of single-season occupancy model with no covariates

# This version uses a Gibbs sampler coded in R

BoccSS0 <- function(y, n, psiPrior=c(1,1), pPrior=c(1,1), 
                    n.chains=3, n.iter=10100, n.burnin=100) {
  startTime <- Sys.time()
  nSites <- length(y)
  if(length(n) == 1)
    n <- rep(n, nSites)
  stopifnot(length(n) == length(y))
  stopifnot(all(n >= y))
  known <- y > 0  # sites known to be occupied
  detected <- sum(y)

  # talk <- round(n.iter / 10) ### 'talk'ing increased time by a factor of 10!
  chain <- matrix(nrow=n.iter, ncol=2)
  colnames(chain) <- c("psi", "p")
  chain[1, ] <- rep(0.5, 2)
  chainList <- vector('list', n.chains)
  for(ch in 1:n.chains)  {
    # cat("\nRunning chain", ch)
    for(i in 2:n.iter) {
      # if(i %% talk == 0)
        # cat(".") ; flush.console()
      # 1. Calculate prob(occupied | y = 0), draw new z vector
      psi.y0 <- (chain[i-1,1] * (1 - chain[i-1,2])^n) /
                    (chain[i-1,1] * (1 - chain[i-1,2])^n + (1 - chain[i-1,1]))
      z <- ifelse(known, 1, rbinom(nSites, 1, psi.y0))
      # 2. Update psi from beta(occupied+1, unoccupied+1)
      chain[i,1] <- rbeta(1, sum(z) + psiPrior[1], sum(z == 0) + psiPrior[2])
      # 3. Update p from beta(detected+1, undetected+1) for occupied sites only.
      chain[i,2] <- rbeta(1, detected + pPrior[1], sum(n[z == 1]) - detected + pPrior[2])
    }
    chainList[[ch]] <- mcmc(chain[(n.burnin+1):n.iter, ], start=n.burnin+1, end=n.iter)
  }
  # cat("\nDone.\n")
  # Diagnostics
  MCMC <- mcmc.list(chainList)
  Rhat <- try(gelman.diag(MCMC, autoburnin=FALSE)$psrf[, 1], silent=TRUE)
  if(inherits(Rhat, "try-error") || !is.finite(Rhat))
    Rhat <- NULL
  
  out <- as.Bwiqid(MCMC,
      header = "Model fitted in R with a Gibbs sampler",
      defaultPlot = "psi")
  attr(out, "call") <- match.call()
  attr(out, "n.chains") <- n.chains
  attr(out, "n.eff") <- effectiveSize(MCMC)
  attr(out, "Rhat") <- Rhat
  attr(out, "timetaken") <- Sys.time() - startTime
  return(out)
}
