

# Bayesian version of CJS models
# This uses runjags

# THIS DOESN'T WORK PROPERLY

# Not exported

BsurvCJS_group <- function(DH, model=list(phi~1, p~1), data=NULL, freq=1, group,
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
  if(length(freq) == 1)
    freq <- rep(freq, nrow(DH))
  if (length(freq) != nrow(DH))
    stop("freq must have a value for each row of the detection history matrix.")

  if(!missing(group))  {
    if(length(group) != nrow(DH))
      stop("Group must have a value for each row of the detection history matrix.")
    group <- as.factor(group)
    nGroup <- nlevels(group)
    # groupNames <- levels(group)
    # data <- as.data.frame(cbind(data, group=rep(groupNames, each=ni)))
  } else {
    group <- factor(rep("all", ni))
    nGroup <- 1
  }
  groupNames <- levels(group)
  data <- as.data.frame(cbind(data, group=rep(groupNames, each=ni)))

  # Convert detection history to 3d array of m-arrays to facilitate use of multinomial likelihood
  mARRAY <- array(0, c(ni, ni+1, nGroup))
  if(nGroup == 1) {
    mARRAY[, , 1] <- ch2mArray(CH=DH, freq=freq)
  } else {
    for(i in 1:nGroup) {
      DHgrp <- subset(DH, group==groupNames[i])
      freqgrp <- subset(freq, group==groupNames[i])
      mARRAY[, , i] <- ch2mArray(CH=DHgrp, freq=freqgrp)
    }
  }

  # Standardise the model:
  model <- stdModel(model, defaultModel=list(phi=~1, p=~1))

  # Standardize the data
  dataList <- stddata(data, NULL)
  dataList$.Time <- as.vector(scale(1:ni)) /2
  dataList$.time <- as.factor(1:ni)

  # Set up model matrices
  phiDf <- selectCovars(model$phi, dataList, ni*nGroup)
  phiMat <- model.matrix(model$phi, phiDf)
  phiK <- ncol(phiMat)
  pDf <- selectCovars(model$p, dataList, ni*nGroup)
  pMat <- model.matrix(model$p, pDf)
  pK <- ncol(pMat)
  K <- phiK + pK
  if(nrow(phiMat) != ni*nGroup || nrow(pMat) != ni*nGroup)
    stop("Missing values not allowed in covariates.")
  # Turn model matrices into arrays, one page per group:
  phiArray <- array(unlist(split(phiMat, data$group)), c(ni, ncol(phiMat), nGroup))
  pArray <- array(unlist(split(pMat, data$group)), c(ni, ncol(pMat), nGroup))
  # Get names for coefficiants
  betaNames <- c(paste0("p:", colnames(pMat)),
                  paste0("phi:", colnames(phiMat)))

  # Run MLE version to get starting values
  nll <- function(param){
    phiBeta <- param[1:phiK]
    pBeta <- param[(phiK+1):K]
    phiProb <- plogis(phiMat %*% phiBeta)
    pProb <- plogis(pMat %*% pBeta)
    if(any(pProb * phiProb == 1))
      return(.Machine$double.max)
    if(nGroup == 1) {
      nll <-  -sum(mARRAY[, , 1] * log(qArray(phiProb, pProb)), na.rm=TRUE)
    } else {
      nll <- numeric(nGroup)
      for(i in 1:nGroup) {
        phiProb0 <- phiProb[data$group == groupNames[i]]
        pProb0 <- pProb[data$group == groupNames[i]]
        nll[i] <- -sum(mARRAY[, , i] * log(qArray(phiProb0, pProb0)), na.rm=TRUE)
      }
    }
    return(min(sum(nll), .Machine$double.xmax))
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
        betaPhi[i] ~ dnorm(0, 0.01)
      }
      for(i in 1:pK)  {
        betaP[i] ~ dnorm(0, 0.01)
      }
      for(grp in 1:nGroup) {
        for(t in 1:(nocc-1)) {
          logit(phi[t, grp]) <- sum(betaPhi[] * phiMat[t, , grp])
          logit(p[t, grp]) <- sum(betaP[] * pMat[t, , grp])
          # Multinomial likelihood
          marr[t, 1:nocc, grp] ~ dmulti(pr[t, , grp], rel[t, grp])
          # Cell probs of the m-array
          q[t, grp] <- 1 - p[t, grp]  # prob of non-recapture
          # main diagonal
          pr[t, t, grp] <- phi[t, grp] * p[t, grp]
          # above main diagonal
          for(j in (t+1):(nocc-1)) {
            pr[t, j, grp] <- prod(phi[t:j, grp]) * prod(q[t:(j-1), grp]) * p[j, grp]
          }
          # below main diagonal
          for(j in 1:(t-1)) {
            pr[t, j, grp] <- 0
          }
          # last column, prob of never recaptured
          pr[t, nocc, grp] <- 1 - sum(pr[t, 1:(nocc-1), grp])
        }
      }
    }
  "

   # organise the data:
  jagsData <- list(nocc = ni+1, rel=apply(mARRAY, c(1, 3), "sum"), nGroup=nGroup,
                      pK = pK, phiK = phiK, pMat=pArray, phiMat=phiArray)
  if(!priorOnly)
      jagsData$marr <- mARRAY

  inits <- function() {
    start1 <- start * runif(K, 0.9, 1.1)
    list(betaPhi = start1[1:phiK], betaP = start1[(phiK+1):K])
  }
  wanted <- c("phi", "p", "betaPhi", "betaP")

  # Run the model:
  resB <- autorun.jags(modeltext, wanted, jagsData, n.chains=3, inits)#, ...)

  out <- as.Bwiqid(resB,
      header = "Model fitted in JAGS with runjags::autorun.jags")
  attr(out, "defaultPlot") <- names(out)[K+1]
  names(out)[1:K] <- betaNames
  return(out)
}

