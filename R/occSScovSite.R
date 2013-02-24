# Single season occupancy with site covariates (not survey covariates)

occSScovSite <- function(DH, psi=~1, p=~1, data=NULL) {
  # single-season occupancy models with site-specific covatiates
  # ** DH is detection data in a 1/0/NA matrix, sites in rows, detection occasions
  #    in columns..
  # ** psi and p are one-sided formulae describing the model
  # ** data a data frame with the site covariates.

  nSites <- nrow(DH)
  nSurv <- ncol(DH)
  if(is.null(rownames(DH)))
    rownames(DH) <- 1:nSites
  # Sanity checks:
  if (nSurv < 2)
    stop("More than one survey occasion is needed")
  if(!is.null(data))  {
    ddf <- as.data.frame(data)
    if(any(is.na(ddf)))
      stop("Missing site covariates are not allowed.")
    if(nrow(ddf) != nSites)
      stop("'data' must have a row for each site.")
  } else {
    ddf <- data.frame(.dummy = rep(NA, nSites))
    # model.matrix needs a data frame, NULL won't do.
  }

  psiMat <- model.matrix(psi, ddf)
  psiK <- ncol(psiMat)
  pMat <- model.matrix(p, ddf)
  pK <- ncol(pMat)
  K <- psiK + pK

  beta.mat <- matrix(NA_real_, K, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c(
    paste("psi:", colnames(psiMat)),
    paste("p:", colnames(pMat)))
  lp.mat <- matrix(NA_real_, nSites*2, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c(
    paste("psi:", rownames(DH), sep=""),
    paste("p:", rownames(DH), sep=""))
  AIC <- NA_real_

  nll <- function(param){
    # This is for constant p and psi, modify for other models.
    psiBeta <- param[1:psiK]
    pBeta <- param[(psiK+1):K]
    psiProb <- plogis(psiMat %*% psiBeta)
    pProb <- plogis(pMat %*% pBeta)
    p.DH <- sweep(DH, 1, pProb, "*") + sweep((1-DH), 1, (1-pProb), "*")
    llh <- sum(log(psiProb * apply(p.DH, 1, prod, na.rm=TRUE) + 
          (1 - psiProb) * (rowSums(DH, na.rm=TRUE) == 0)))
    return(min(-llh, .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  res <- nlm(nll, param, hessian=TRUE)
  if(res$code < 3)  {  # exit code 1 or 2 is ok.
    beta.mat[,1] <- res$estimate
    lp.mat[, 1] <- c(psiMat %*% beta.mat[1:psiK, 1],
                     pMat %*% beta.mat[(psiK+1):K, 1])
    AIC <- 2*res$minimum + 2 * K
    if (det(res$hessian) > 0) {
      varcov <- solve(res$hessian)
      SE <- sqrt(diag(varcov))
      beta.mat[, 2] <- SE
      crit <- qnorm(c(0.025, 0.975))
      beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      SElp <- c(sqrt(diag(psiMat %*% varcov[1:psiK, 1:psiK] %*% t(psiMat))),
                sqrt(diag(pMat %*% varcov[(psiK+1):K, (psiK+1):K] %*% t(pMat))))
      lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
    }
  }
  out <- list(beta = beta.mat, real = plogis(lp.mat), AIC = AIC)
  return(out)
}


