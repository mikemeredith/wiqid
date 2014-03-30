
# Bayes version of single-season occupancy model with no covariates

BoccSS0 <-
function(y, n, priorOnly=FALSE, ...) {
  # n is a vector with the number of occasions at each site.
  # y is a vector with the number of detections at each site.
  # ci is the required confidence interval.
  
  stopifnot(testjags(silent=TRUE)$JAGS.found)
  if (detectCores() > 3)
    runjags.options("method"="parallel")

  if (priorOnly)
    warning("The prior distributions will be produced, not the posterior distributions!")

  # Call occSS0 to do sanity checks and get starting values
  start <- occSS0(y=y, n=n)$real[, 1]
  if(length(n) == 1)
    n <- rep(n, length(y))

  # Define the model
  modeltext <- "
    model {

      # Priors
      psi ~ dunif(0, 1)
      p ~ dunif(0, 1)

      # Likelihood
      for (i in 1:R) {
        z[i] ~ dbern(psi)
        y[i] ~ dbinom(z[i] * p, n[i])
      }
    } "

  # Prepare the bits:
  jagsData <- list(n=n, R=length(y))
  if(!priorOnly)
    jagsData$y <- y
  inits <- list(psi=start[1], p=start[2], z=(y > 0)*1)
  wanted <- c("psi", "p")
  
  # Run the model:
  resB <- autorun.jags(modeltext, wanted, jagsData, n.chains=3, inits, ...)
  
  return(as.Bwiqid(resB, 
      header = "Model fitted in JAGS with runjags::autorun.jags",
      defaultPlot = "psi"))
}
