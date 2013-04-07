# Calculation of likelihood for m-array of recaptures 


qArray <- function(phi, p) {
  # Calculates the matrix of multinomial cell probabilities
  #   corresponding to an m-array.
  # phi = vector of apparent survival probabilities
  # p = vector of recapture probabilities

  n <- length(phi)

  # NO SANITY CHECKS, calling function must take care of this.
  # if(length(p) != n)
   # stop("p and phi must have the same length.")
  # if(n < 2)
   # stop("More than one recapture occasion is needed.")

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

cjs0 <- function(mArray, ci = 0.95) {
  # phi(.) p(.) model for Cormack-Joly-Seber estimation of apparent survival.
  # nArray is capture-recapture data in m-array format, with zeros in the lower
  #   triangle.
  # ci is the required confidence interval
  # Original code from Ruth King, but cf. Kery & Schaub
  # Largely altered by MM

  # Sanity checks:
  ni <- nrow(mArray)
  if (ni < 2)
    stop("More than one recapture occasion is needed")
  if (ncol(mArray) != ni + 1 || any(mArray[lower.tri(mArray)] != 0))
    stop("mArray is not a valid m-array format")
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  beta.mat <- matrix(NA_real_, 2, 3)
  AIC <- NA_real_

  nll <- function(param){
    # This is for constant p and phi, modify for other models.
    p <- rep(plogis(param[1]), ni)
    phi <- rep(plogis(param[2]), ni)
    if(any(p * phi == 1))
      return(.Machine$double.max)
    # Output the negative log(likelihood) value:
    return(-sum(mArray * log(qArray(p, phi)), na.rm=TRUE))
  }

  # Run mle estimation with nlm:
  param <- c(0, 0)
  res <- nlm(nll, param, hessian=TRUE)
  if(res$code < 3)  {  # exit code 1 or 2 is ok.
    beta.mat[,1] <- res$estimate
    AIC <- 2*res$minimum + 4
    if (det(res$hessian) > 0) {
      SE <- sqrt(diag(solve(res$hessian)))
      beta.mat[, 2:3] <- sweep(outer(SE, crit), 1, res$estimate, "+")
    }
  }
  out.mat <- plogis(beta.mat)
  colnames(out.mat) <- c("est", "lowCI", "uppCI")
  rownames(out.mat) <- c("phiHat", "pHat")
  attr(out.mat, "AIC") <- AIC
  return(out.mat)
}
# .................................................................

cjsTime <- function(mArray, phi=~1, p=~1, data=NULL, ci = 0.95) {
  # phi(t) p(t) model or models with time covariates for Cormack-Joly-Seber
  # estimation of apparent survival.
  # ** mArray is capture-recapture data in m-array format, with zeros in the lower
  #   triangle.
  # ** phi and p are one-sided formulae describing the model
  # ** data a data frame with the covariates.
  # ** ci is required confidence interval.

  # Sanity checks:
  ni <- nrow(mArray)
  if (ni < 2)
    stop("More than one recapture occasion is needed")
  if (ncol(mArray) != ni + 1 || any(mArray[lower.tri(mArray)] != 0))
    stop("mArray is not a valid m-array format")
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  data$time <- as.factor(1:ni)
  ddf <- as.data.frame(data)
  phiMat <- model.matrix(phi, ddf)
  phiK <- ncol(phiMat)
  pMat <- model.matrix(p, ddf)
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
  AIC <- NA_real_

  nll <- function(param){
    # This is for constant p and phi, modify for other models.
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
    AIC <- 2*res$minimum + 2 * K
    if (det(res$hessian) > 0) {
      varcov <- solve(res$hessian)
      SE <- sqrt(diag(varcov))
      beta.mat[, 2] <- SE
      beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      SElp <- c(sqrt(diag(phiMat %*% varcov[1:phiK, 1:phiK] %*% t(phiMat))),
                sqrt(diag(pMat %*% varcov[(phiK+1):K, (phiK+1):K] %*% t(pMat))))
      lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
    }
  }
  out <- list(beta = beta.mat, real = plogis(lp.mat), AIC = AIC)
  return(out)
}


