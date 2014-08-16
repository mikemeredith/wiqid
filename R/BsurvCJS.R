

# Bayesian version of CJS models

# This uses runjags

BsurvCJS <- function(DH, model=list(phi~1, p~1), data=NULL, freq=1,
    priorOnly=FALSE, ...) {
  # phi(t) p(t) model or models with time covariates for Cormack-Joly-Seber
  # estimation of apparent survival.
  # ** DH is detection history matrix/data frame, animals x occasions.
  # ** freq is vector of frequencies for each detection history
  # ** model is a list of 2-sided formulae for psi and p; can also be a single
  #   2-sided formula, eg, model = psi ~ habitat.
  # ** data a data frame with the covariates.
  # ** ci is required confidence interval.

  stopifnot(testjags(silent=TRUE)$JAGS.found)
  if (detectCores() > 3)
    runjags.options("method"="parallel")
  runjags.options("rng.warning"=FALSE)
  
  # Sanity checks:
  if (priorOnly)
    warning("The prior distributions will be produced, not the posterior distributions!")
    
  ni <- ncol(DH) - 1  # number of survival intervals and REcapture occasions
  stopifnot(is.null(data) || nrow(data) == ni)

  # Convert detection history to m-array to facilitate use of multinomial likelihood
  mArray <- ch2mArray(CH=DH, freq=freq)

  # Standardise the model:
  model <- stdModel(model, defaultModel=list(phi=~1, p=~1))

  # Standardize the data
  dataList <- stddata(data, NULL)
  dataList$.Time <- as.vector(scale(1:ni)) /2
  dataList$.time <- as.factor(1:ni)

  # Set up model matrices
  phiDf <- selectCovars(model$phi, dataList, ni)
  phiMat <- model.matrix(model$phi, phiDf)
  phiK <- ncol(phiMat)
  pDf <- selectCovars(model$p, dataList, ni)
  pMat <- model.matrix(model$p, pDf)
  pK <- ncol(pMat)
  K <- phiK + pK
  if(nrow(phiMat) != ni || nrow(pMat) != ni)
    stop("Missing values not allowed in covariates.")

  # Run MLE version to get starting values
  nll <- function(param){
    phiBeta <- param[1:phiK]
    pBeta <- param[(phiK+1):K]
    phiProb <- plogis(phiMat %*% phiBeta)
    pProb <- plogis(pMat %*% pBeta)
    if(any(pProb * phiProb == 1))
      return(.Machine$double.max)
    # Output the negative log(likelihood) value:
    nll <- -sum(mArray * log(qArray(phiProb, pProb)), na.rm=TRUE)
    return(min(nll, .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  res <- nlm(nll, param)
  # if(res$code > 2)   # exit code 1 or 2 is ok.
    # stop("MLE estimation failed.")
  start <- res$estimate

  # Do the model:
  modeltext <- "
    model{
      # priors for beta parameters
      for(i in 1:phiK)  {
        phiBeta[i] ~ dnorm(0, 0.01)
      }
      for(i in 1:pK)  {
        pBeta[i] ~ dnorm(0, 0.01)
      }
      # Calculate p and phi
      for(t in 1:(nocc-1)) {
        logit(phi[t]) <- sum(phiBeta[] * phiMat[t, ])
        logit(p[t]) <- sum(pBeta[] * pMat[t, ])
      }
      # Multinomial likelihood
      for(t in 1:(nocc-1)) {
        marr[t, 1:nocc] ~ dmulti(pr[t, ], rel[t])
      }
      # Cell probs of the m-array
      for(t in 1:(nocc-1)) {
        q[t] <- 1 - p[t]  # prob of non-recapture
        # main diagonal
        pr[t, t] <- phi[t] * p[t]
        # above main diagonal
        for(j in (t+1):(nocc-1)) {
          pr[t, j] <- prod(phi[t:j]) * prod(q[t:(j-1)]) * p[j]
        }
        # below main diagonal
        for(j in 1:(t-1)) {
          pr[t, j] <- 0
        }
      }
      # last column, prob of never recaptured
      for(t in 1:(nocc-1)) {
        pr[t, nocc] <- 1 - sum(pr[t, 1:(nocc-1)])
      }
    }
  "
    
   # organise the data:
  jagsData <- list(nocc = ncol(mArray), rel=rowSums(mArray), 
                      pK = pK, phiK = phiK, pMat=pMat, phiMat=phiMat)
  if(!priorOnly)
      jagsData$marr <- mArray
      
  inits <- function() {
    start1 <- start * runif(K, 0.9, 1.1)
    list(phiBeta = start1[1:phiK], pBeta = start1[(phiK+1):K])
  }
  wanted <- c("phi", "p")
  
  # Run the model:
  resB <- autorun.jags(modeltext, wanted, jagsData, n.chains=3, inits, ...)
    
  return(as.Bwiqid(resB, 
      header = "Model fitted in JAGS with runjags::autorun.jags",
      defaultPlot = "phi1"))
}

