
# Bayes version of single-season occupancy model with no covariates

BoccSS0 <-
function(y, n, numSavedSteps=1e4, thinSteps=1, burnInSteps = 1e3, priorOnly=FALSE) {
  # n is a vector with the number of occasions at each site.
  # y is a vector with the number of detections at each site.
  # ci is the required confidence interval.
  
  if (priorOnly)
    warning("The prior distributions will be produced, not the posterior distributions!")

  # Call occSS0 to do sanity checks and get starting values
  start <- occSS0(y=y, n=n)$real[, 1]
  if(length(n) == 1)
    n <- rep(n, length(y))

  # Define the model
  modelFile <- tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".txt")
  model <- "
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
  writeLines(model, con=modelFile)

  # Prepare the bits:
  jagsData <- list(n=n, R=length(y))
  if(!priorOnly)
    jagsData$y <- y
  inits <- list(psi=start[1], p=start[2], z=(y > 0)*1)
  wanted <- c("psi", "p")
  
  # Create the model and run:
  jm <-jags.model(modelFile, jagsData, inits,
            n.chains = 1, n.adapt=1000, quiet=FALSE)
  update(jm, n.iter=burnInSteps)
  codaSamples <- coda.samples(jm, variable.names=wanted, 
                        n.iter= numSavedSteps * thinSteps, thin=thinSteps)
  
  out <- as.data.frame(as.matrix(codaSamples))
  class(out) <- c("Bwiqid", class(out))
  attr(out, "n.eff") <- effectiveSize(codaSamples)
  attr(out, "data") <- list(y = y, n = n)
  attr(out, "defaultPlot") <- "psi"
  return(out)
}
