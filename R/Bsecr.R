
# Bayesian version of secr to work with stoats data

Bsecr0 <- function(capthist, buffer = 100, start=NULL, nAug = NA,
                    nChains=3, numSavedSteps=1e4, thinSteps=1, burnInSteps = 0,
                    priorOnly=FALSE) {
  stopifnot(inherits(capthist, "capthist"))
  if (priorOnly)
    warning("The prior distributions will be produced, not the posterior distributions!")

  traps <- traps(capthist)
  J <- nrow(traps)
  xl <- min(traps$x) - buffer
  xu <- max(traps$x) + buffer
  yl <- min(traps$y) - buffer
  yu <- max(traps$y) + buffer
  A <- (xu-xl)*(yu-yl) / 1e4 # ha

  # Get starting values, etc from secr.fit
  if(!is.null(start) && inherits(start, "secr")) {
    mle.res <- predict(start)
  } else {
    cat("Running secr.fit to get starting values...") ; flush.console()
    mle.res <- predict(secr.fit(capthist, buffer=buffer, trace=FALSE))
    cat("done\n") ; flush.console()
  }
  if(is.na(nAug))
    nAug <- ceiling(1.5 * mle.res[1, 5] * A)
  # maxSig2 <- (1.1 * mle.res[3, 5]) ^ 2
  maxSig2 <- (2.2 * mle.res[3, 5]) ^ 2
  lam0start <- mle.res[2,2]
  # sigma2start <- mle.res[3,2]^2 # sigma2 = 2 * sigma^2
  sigma2start <- 2 * mle.res[3,2]^2
  psistart  <- (mle.res[1,2] * A) / nAug
    
  # Convert capture histories into an Animals x Traps matrix
  nInd <- dim(capthist)[1]
  nOcc <- dim(capthist)[2]
  yMat <- matrix(0, nAug, J)
  if(length(dim(capthist)) == 3) {
    yMat[1:nInd,] <- apply(capthist, c(1,3), sum)
  } else {
    for(i in 1:nInd)
      for(j in 1:nOcc)
        yMat[i, capthist[i, j]] <- yMat[i, capthist[i, j]] + 1
  }
  
  # Get initial locations of animals
  SX <- SY <- numeric(nInd)
  for(i in 1:nInd) {
    where <- colMeans(traps[which(yMat[i, ] > 0), , drop=FALSE])
    SX[i] <- where[1]
    SY[i] <- where[2]
  }

  # Define the model
  modelFile <- tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".txt")
  # modelFile <- "model.tst"
  model <- "
    model {
    sigma2 ~ dunif(0, maxSig2)      # need to set good max
    sigma <- sqrt(sigma2 / 2)
    lam0 ~ dgamma(0.1, 0.1)
    psi ~ dunif(0, 1)
    for (i in 1:M){         #loop through the augmented population
      z[i] ~ dbern(psi)     #state of individual i (real or imaginary)
      SX[i] ~ dunif(xl, xu) #priors for the activity centre for each individual
      SY[i] ~ dunif(yl, yu) #  xl, yl = lower coordinate; xu, yu = upper value
      for(j in 1:J) {       #loop through the J camera trap locations
         Dsq[i,j] <- pow(SX[i]-trapmat[j,1], 2) + pow(SY[i]-trapmat[j,2],2)
         g[i,j] <- exp(-Dsq[i,j]/sigma2)
         lambda[i,j] <- K * g[i,j] * lam0 * z[i]
         y[i,j] ~ dpois(lambda[i,j])
      }
    }
    N <- sum(z[1:M]) # derive number (check against M)
    D <- N / A       # derive density
  }   "
  writeLines(model, con=modelFile)

   # organise the data:
  jagsData <- list(M = nAug, xl=xl, xu=xu, yl=yl, yu=yu, J=J, trapmat=traps, K=nOcc,
                    A = A, maxSig2 = maxSig2)
  if (!priorOnly)
    jagsData$y <- yMat
  inits <- function() {list(z=rep(1, nAug), 
                        SX=c(SX, runif(nAug-nInd, xl, xu)),
                        SY=c(SY, runif(nAug-nInd, yl, yu)),
                        sigma2=sigma2start, lam0=lam0start,
                        psi=psistart)}
  wanted <- c("D", "lam0", "sigma")
  # Create the model and run:
  jm <-jags.model(modelFile, jagsData, inits,
            n.chains = nChains, n.adapt=500, quiet=FALSE)
  if(burnInSteps > 0)
    update(jm, n.iter=burnInSteps)
  codaSamples <- coda.samples(jm, variable.names=wanted, 
                        n.iter= ceiling(numSavedSteps * thinSteps / nChains), thin=thinSteps)
  # gelman.diag(codaSamples)
  # effectiveSize(codaSamples)

  out <- as.data.frame(as.matrix(codaSamples))
  # names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  attr(out, "n.eff") <- effectiveSize(codaSamples)
  attr(out, "Rhat") <- gelman.diag(codaSamples)$psrf[, 1]
  attr(out, "defaultPlot") <- "D"
  # attr(out, "data") <- ???
  return(out)
}
