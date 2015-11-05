# Estimation of apparent survival - CJS models

# These functions allow for Adult and Juvenile survival to differ,
#  eg. when birds are ringed as nestlings.

# Data should be organised as 2 detection history matrices, one for juveniles and one for
#   adults (when juveniles are recaptured they are already adults).

qArrayAJ <- function(phi, p, phiJ=phi) {
  # Calculates the matrix of multinomial cell probabilities
  #   corresponding to an m-array.
  # phi = vector of apparent survival probabilities
  # phiJ = vector of apparent survival probabilities for juveniles
  # p = vector of recapture probabilities
  # NO SANITY CHECKS, calling function must take care of this.

  n <- length(phi)

  q <- diag(as.vector(p * phiJ), n, n+1)  # Create n x n+1 matrix and fill diagonal
  for (i in 1:(n-1)){ # Fill the upper triangle
    for (j in (i+1):n) {
      q[i,j] <- phiJ[i]*prod(phi[(i+1):j])*prod(1-p[i:(j-1)])*p[j]
    }
  }
  q[, n+1] <- 1 - rowSums(q, na.rm=TRUE)  # Add the last column and return
  return(q)
}
# ..........................................................................

survCJSaj <- function(DHj, DHa=NULL, model=list(phiJ~1, phiA~1, p~1), data=NULL,
    freqj=1, freqa=1, ci = 0.95, link=c("logit", "probit")) {
  # phi(t) p(t) model or models with time covariates for Cormack-Joly-Seber
  # estimation of apparent survival.
  # ** DHj is detection history matrix/data frame, animals x occasions, for animals marked as juveniles; DHa (optional) has detection histories for animals marked as adults.
  # ** freqj and freqa are vectors of frequencies for each detection history
  # ** model is a list of 2-sided formulae for psiJ, psiA and p; can also be a single
  #   2-sided formula, eg, model = psiJ ~ habitat.
  # ** data a data frame with the covariates.
  # ** ci is required confidence interval.

    if(match.arg(link) == "logit") {
    plink <- plogis
  } else {
    plink <- pnorm
  }

  # Sanity checks:
  # Check DHj and DHa have same no. of columns ...
  nocc <- ncol(DHj)
  ni <- nocc - 1  # number of survival intervals and REcapture occasions
  stopifnot(is.null(DHa) || ncol(DHa) == nocc)
  if (length(freqj) == 1)
    freqj <- rep(freqj, nrow(DHj))
  # if (length(freqa) == 1)
    # freqa <- rep(freqa, nrow(DHa))  # Not needed

  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  # Deal with grownup juveniles, do m-array for these:
  grown <- DHj
  # Remove first capture
  getFirst <- function(x) min(which(x == 1))
  first <- apply(DHj, 1, getFirst)
  for(i in 1:nrow(grown))
    grown[i, first[i]] <- 0
  marrayA <- ch2mArray(grown, freqj)

  # Do m-array for juvenile juveniles
  ma <- matrix(0, nocc, nocc+1)
  for(i in 1:nrow(DHj)) {
    cht <- which(DHj[i, ] != 0) # When was animal caught?
    # Fill in release/recapture data
    # we are only interested in the first recapture
    if(length(cht) > 1)
      ma[cht[1], cht[2]] <- ma[cht[1], cht[2]] + freqj[i]
  }
  # Juveniles never seen again:
  ringed <- tapply(freqj, first, sum)
  ma[, nocc+1] <- c(ringed, 0) - rowSums(ma)
  marrayJ <- ma[-nocc, -1]

  # Add data for adults
  if(!is.null(DHa))
    marrayA <- marrayA + ch2mArray(DHa, freqa)

  # Standardise the model:
  model <- stdModel(model, defaultModel=list(phiJ=~1, phiA=~1, p=~1))

  # Standardize the data
  dataList <- stddata(data, NULL)
  dataList$.Time <- as.vector(scale(1:ni)) /2
  dataList$.time <- as.factor(1:ni)

  # Set up model matrices
  phiADf <- selectCovars(model$phiA, dataList, ni)
  phiAMat <- model.matrix(model$phiA, phiADf)
  phiAK <- ncol(phiAMat)
  phiJDf <- selectCovars(model$phiJ, dataList, ni)
  phiJMat <- model.matrix(model$phiJ, phiJDf)
  phiJK <- ncol(phiJMat)
  pDf <- selectCovars(model$p, dataList, ni)
  pMat <- model.matrix(model$p, pDf)
  pK <- ncol(pMat)
  K <- phiAK + phiJK + pK
  parID <- rep(1:3, c(phiAK, phiJK, pK))
  if(nrow(phiAMat) != ni || nrow(phiJMat) != ni || nrow(pMat) != ni)
    stop("Missing values not allowed in covariates.")

  # Objects to hold results
  beta.mat <- matrix(NA_real_, K, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c(
    paste("phiA:", colnames(phiAMat)),
    paste("phiJ:", colnames(phiJMat)),
    paste("p:", colnames(pMat)))
  lp.mat <- matrix(NA_real_, ni*3, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c(
    paste("phiA", 1:ni, sep=""),
    paste("phiJ", 1:ni, sep=""),
    paste("p", 1:ni, sep=""))
  logLik <- NA_real_
  varcov <- NULL

  nll <- function(param){
    phiABeta <- param[parID==1]
    phiJBeta <- param[parID==2]
    pBeta <- param[parID==3]
    phiAProb <- plink(phiAMat %*% phiABeta)
    phiJProb <- plink(phiJMat %*% phiJBeta)
    pProb <- plink(pMat %*% pBeta)
    if(any(pProb * phiAProb == 1) || any(pProb * phiJProb == 1) )
      return(.Machine$double.max)
    # Calculate the negative log(likelihood) value:
    return(min(-sum(marrayA  * log(qArrayAJ(phiAProb, pProb, phiAProb)),   # adults
            marrayJ * log(qArrayAJ(phiAProb, pProb, phiJProb)), na.rm=TRUE),   # juveniles
          .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  res <- nlm(nll, param, hessian=TRUE, stepmax=10) # 2015-03-01
  if(res$code > 2)   # exit code 1 or 2 is ok.
    warning(paste("Convergence may not have been reached (nlm code", res$code, ")"))

  # Organise the output
  beta.mat[,1] <- res$estimate
  lp.mat[, 1] <- c(phiAMat %*% beta.mat[parID==1, 1],
                   phiJMat %*% beta.mat[parID==2, 1],
                   pMat %*% beta.mat[parID==3, 1])
  # varc <- try(solve(res$hessian), silent=TRUE)
  varcov0 <- try(chol2inv(chol(res$hessian)), silent=TRUE)
  # if (!inherits(varc, "try-error") && all(diag(varc) > 0)) {
  if (!inherits(varcov0, "try-error")) {
    varcov <- varcov0
    SE <- suppressWarnings(sqrt(diag(varcov)))
    beta.mat[, 2] <- SE
    beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
    temp <- c(sqrt(diag(phiAMat %*% varcov[parID==1, parID==1] %*% t(phiAMat))),
             sqrt(diag(phiJMat %*% varcov[parID==2, parID==2] %*% t(phiJMat))),
             sqrt(diag(pMat %*% varcov[parID==3, parID==3] %*% t(pMat))))
    if(all(temp >= 0))  {
      SElp <- sqrt(temp)
      lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
      logLik <- -res$minimum
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plink(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=sum(marrayJ, marrayA)))
  class(out) <- c("wiqid", "list")
  return(out)
}


