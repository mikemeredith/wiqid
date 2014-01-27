
# Maximum likelihood estimation for Cormack-Jolly-Seber apparent-survival models
#    of type phi(t) p(t) or with time covariates
# Uses detection history matrix as input, not m-array.
# Uses multinomial likelihood.

qArray <- function(phi, p) {
  # Calculates the matrix of multinomial cell probabilities
  #   corresponding to an m-array.
  # phi = vector of apparent survival probabilities
  # p = vector of recapture probabilities
  # NO SANITY CHECKS, calling function must take care of this.

  n <- length(phi)

  q <- diag(as.vector(p * phi), n, n+1)  # Create n x n+1 matrix and fill diagonal
  for (i in 1:(n-1)){ # Fill the upper triangle
    for (j in (i+1):n) {
      q[i,j] <- prod(phi[i:j])*prod(1-p[i:(j-1)])*p[j]
    }
  }
  q[, n+1] <- 1 - rowSums(q, na.rm=TRUE)  # Add the last column and return
  return(q) 
}
# ..........................................................................

survCJS <- function(DH, model=list(phi~1, p~1), data=NULL, freq=1, ci = 0.95) {
  # ** DH is detection history matrix/data frame, animals x occasions.
  # ** freq is vector of frequencies for each detection history
  # ** model is a list of 2-sided formulae for psi and p; can also be a single
  #   2-sided formula, eg, model = psi ~ habitat.
  # ** data a data frame with the covariates.
  # ** ci is required confidence interval.

  # Sanity checks:  for DH??
  ni <- ncol(DH) - 1  # number of survival intervals and REcapture occasions
  if(!is.null(data) && nrow(data) != ni)
    stop("The 'data' argument is not a valid data frame.")
  
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  # Convert detection history to m-array to facilitate use of multinomial likelihood
  mArray <- ch2mArray(CH=DH, freq=freq)
  
  # Standardise the model:
  model <- stdModel(model, defaultModel=list(phi=~1, p=~1))

  data$.time <- as.factor(1:ni)  # add .Time trend
  ddf <- as.data.frame(data)  # Should standardise
  phiMat <- model.matrix(model$phi, ddf)
  phiK <- ncol(phiMat)
  pMat <- model.matrix(model$p, ddf)
  pK <- ncol(pMat)
  K <- phiK + pK

  beta.mat <- matrix(NA_real_, K, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c(
    paste("phi:", colnames(phiMat)),
    paste("p:", colnames(pMat)))
  lp.mat <- matrix(NA_real_, ni*2, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c(
    paste("phi", 1:ni, sep=""),
    paste("p", 1:ni, sep=""))
  logLik <- NA_real_

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
  res <- nlm(nll, param, hessian=TRUE)
  if(res$code < 3)  {  # exit code 1 or 2 is ok.
    beta.mat[,1] <- res$estimate
    lp.mat[, 1] <- c(phiMat %*% beta.mat[1:phiK, 1],
                     pMat %*% beta.mat[(phiK+1):K, 1])
    varcov <- try(solve(res$hessian), silent=TRUE)
    if (!inherits(varcov, "try-error") && all(diag(varcov) > 0)) {
      SE <- sqrt(diag(varcov))
      beta.mat[, 2] <- SE
      beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      SElp <- c(sqrt(diag(phiMat %*% varcov[1:phiK, 1:phiK] %*% t(phiMat))),
                sqrt(diag(pMat %*% varcov[(phiK+1):K, (phiK+1):K] %*% t(pMat))))
      lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
      logLik <- -res$minimum
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=sum(mArray)))
  class(out) <- c("wiqid", "list")
  return(out)
}


