# Bayesian Single season occupancy with site and survey covariates.

# This doesn'twork properly with the weta data with site and survey covars.

# DO NOT EXPORT

BoccSScov <- function(DH, model=list(psi~1, p~1), data=NULL, 
    numSavedSteps=1e4, thinSteps=1, burnInSteps = 1e3) {
  # single-season occupancy models with site and survey covatiates
  # ** DH is detection data in a 1/0/NA matrix or data frame, sites in rows, 
  #    detection occasions in columns..
  # ** model is a list of 2-sided formulae for psi and p; can also be a single
  #   2-sided formula, eg, model = psi ~ habitat.
  # ** data is a DATA FRAME with single columns for site covariates and a column for each survey occasion for each survey covariate.

  # Standardise the model:
  if(inherits(model, "formula"))
    model <- list(model)
  model <- stdform (model)
  model0 <- list(psi=~1, p=~1)
  model <- replace (model0, names(model), model)

  # Summarize detection history
  site.names <- rownames(DH)
  DH <- as.matrix(DH)
  nSites <- nrow(DH)
  nSurv <- ncol(DH)
  if(is.null(site.names))
    site.names <- 1:nSites

  # Convert the covariate data frame into a list
  if(is.data.frame(data))
    data <- stddata(data, nSurv)
  
  # Sanity checks:
  if (nSurv < 2)
    stop("More than one survey occasion is needed")
  survey.done <- !is.na(as.vector(DH))
  DHvec <- as.vector(DH)[survey.done]
  # siteID <- as.factor(row(DH))[survey.done]
  siteID <- row(DH)[survey.done]
  survID <- as.factor(col(DH))[survey.done]
  if(!is.null(data))  {
    covLen <- lapply(data, length)
    # Covars for occupancy:
    psiList <- data[covLen == nSites]
    psiList <- lapply(psiList, function(x) if(is.numeric(x)) scale(x) else x)
    psiDf <- as.data.frame(psiList)
    if(any(is.na(psiDf)))
      stop("Missing site covariates are not allowed.")
    # Covars for probability of detection:
    pList <- data[covLen == nSites | covLen == nSites*nSurv]
    pList$.Time <- col(DH) 
    pList$.Time2 <- col(DH)^2
    pList <- lapply(pList, as.vector)
    pList <- lapply(pList, function(x) if(is.numeric(x)) scale(x) else x)
        # .Time and .Time2 standardized separately ??
    pList$.time <- as.factor(col(DH))
    pDf <- as.data.frame(pList)[survey.done, ]
    if(any(is.na(pDf)))
      stop("Missing survey covariates are not allowed when a survey was done.")
  } else {
    psiDf <- data.frame(.dummy = rep(NA, nSites))
    # model.matrix needs a data frame, NULL won't do.
    pDf <- data.frame(.time = as.factor(survID))
    pDf$.Time <- scale(as.vector(col(DH))[survey.done])
    pDf$.Time2 <- scale(as.vector(col(DH))[survey.done]^2)
  }

  psiModMat <- model.matrix(model$psi, psiDf)
  psiK <- ncol(psiModMat)
  pModMat <- model.matrix(model$p, pDf)
  pK <- ncol(pModMat)
  K <- psiK + pK

  # Do MLE estimate to get starting values
  # Negative log likelihood function
  nll <- function(param){
    psiBeta <- param[1:psiK]
    pBeta <- param[(psiK+1):K]
    psiProb <- as.vector(plogis(psiModMat %*% psiBeta))
    pProb <- plogis(pModMat %*% pBeta)
    Lik1 <- DHvec*pProb + (1-DHvec) * (1-pProb)
    Lik2 <- tapply(Lik1, as.factor(siteID), prod)
    llh <- sum(log(psiProb * Lik2 + 
          (1 - psiProb) * (rowSums(DH, na.rm=TRUE) == 0)))
    return(min(-llh, .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  res <- nlm(nll, param)
  if(res$code > 3)   # exit code 1 or 2 is ok, 3 dodgy but...
    stop("MLE estimation failed")
  start <- res$estimate

  # JAGS model
  modelFile <- tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".txt")
  # modelFile <- "model.txt"
  modeltext <- "
    model{
      # priors for beta parameters
      for(i in 1:psiK)  {
        # psiBeta[i] ~ dnorm(0, 0.01)
        psiBeta[i] ~ dunif(-5, 5)
      }
      for(i in 1:pK)  {
        pBeta[i] ~ dnorm(0, 0.01)
      }
      # Calculate psi, do z
      for(i in 1:nSites) {
        logit(psi[i]) <- sum(psiBeta[] * psiMat[i, ])
        # lp[i] <- sum(psiBeta[] * psiMat[i, ])
        # psi[i] <- 1/(1+exp(-lp[i]))
        z[i] ~ dbern(psi[i])
      }
      for(j in 1:nObs) {
        logit(p[j]) <- sum(pBeta[] * pMat[j, ])
        DHvec[j] ~ dbern(p[j] * z[siteID[j]])
      }
    }
  "
  writeLines(modeltext, con=modelFile)
  jagsData <- list(DHvec=DHvec, nSites = nSites, siteID=siteID,
                nObs=length(DHvec), 
                pK = pK, psiK = psiK, pMat=pModMat, psiMat=psiModMat)
  inits <- function() {list(psiBeta = start[1:psiK],
                        pBeta = start[(psiK+1):K],
                        # z=rep(1, nSites))}
                        z=(rowSums(DH, na.rm=TRUE) > 0)*1)}
  wanted <- c("psi", "psiBeta", "pBeta")
  # Create the model and run:
  jm <-jags.model(modelFile, jagsData, inits,
            n.chains = 3, n.adapt=1000, quiet=FALSE)
  update(jm, n.iter=burnInSteps)
  codaSamples <- coda.samples(jm, variable.names=wanted, 
                        n.iter= numSavedSteps * thinSteps, thin=thinSteps)

  out <- as.data.frame(as.matrix(codaSamples))
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  attr(out, "n.eff") <- effectiveSize(codaSamples)
  attr(out, "data") <- DH
  attr(out, "defaultPlot") <- "psi"
  return(out)

}


