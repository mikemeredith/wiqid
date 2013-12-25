# Calculation of likelihood for m-array of recaptures 

# These functions allow for Adult and Juvenile survival to differ,
#  eg. when birds are ringed as nestlings.

# Data should be organised as 2 m-arrays, one for juveniles and one for
#   adults, noting that when juveniles are recaptured they are already adults.

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

survCJSaj <- function(mArray, mArrayJ, model=list(phiA~1, phiJ~1, p~1), data=NULL, ci = 0.95) {
  # phi(t) p(t) model or models with time covariates for Cormack-Joly-Seber
  # estimation of apparent survival.
  # ** mArray is capture-recapture data in m-array format, with zeros in the lower
  #   triangle.
  # ** model is a list of 2-sided formulae for psi and p; can also be a single
  #   2-sided formula, eg, model = psi ~ habitat.
  # ** data a data frame with the covariates.
  # ** ci is required confidence interval.

  # Sanity checks:
  ni <- nrow(mArray)
  if (ni < 2)
    stop("More than one recapture occasion is needed")
  if (ncol(mArray) != ni + 1 || any(mArray[lower.tri(mArray)] != 0))
    stop("mArray is not a valid m-array format")
  stopifnot(nrow(mArrayJ) == ni)
  stopifnot(ncol(mArrayJ) == ni+1)
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  # Standardise the model:
  if(inherits(model, "formula"))
    model <- list(model)
  model <- stdform (model)
  model0 <- list(phiA=~1, phiJ=~1, p=~1)
  model <- replace (model0, names(model), model)

  data$time <- as.factor(1:ni)
  ddf <- as.data.frame(data)
  phiAMat <- model.matrix(model$phiA, ddf)
  phiAK <- ncol(phiAMat)
  phiJMat <- model.matrix(model$phiJ, ddf)
  phiJK <- ncol(phiJMat)
  pMat <- model.matrix(model$p, ddf)
  pK <- ncol(pMat)
  K <- phiAK + phiJK + pK
  parID <- rep(1:3, c(phiAK, phiJK, pK))

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
    phiAProb <- plogis(phiAMat %*% phiABeta)
    phiJProb <- plogis(phiJMat %*% phiJBeta)
    pProb <- plogis(pMat %*% pBeta)
    if(any(pProb * phiAProb == 1) || any(pProb * phiJProb == 1) )
      return(.Machine$double.max)
    # Calculate the negative log(likelihood) value:
    return(min(-sum(mArray  * log(qArrayAJ(phiAProb, pProb, phiAProb)),   # adults
            mArrayJ * log(qArrayAJ(phiAProb, pProb, phiJProb)), na.rm=TRUE),   # juveniles
          .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  res <- nlm(nll, param, hessian=TRUE)
  if(res$code < 3)  {  # exit code 1 or 2 is ok.
    beta.mat[,1] <- res$estimate
    lp.mat[, 1] <- c(phiAMat %*% beta.mat[parID==1, 1],
                     phiJMat %*% beta.mat[parID==2, 1],
                     pMat %*% beta.mat[parID==3, 1])
    varc <- try(solve(res$hessian), silent=TRUE)
    if (!inherits(varc, "try-error") && all(diag(varc) > 0)) {
      varcov <- varc
      SE <- sqrt(diag(varcov))
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
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=sum(mArray, mArrayJ)))
  class(out) <- c("wiqid", "list")
  return(out)
}


