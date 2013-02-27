# Single season occupancy with site covariates (not survey covariates)

occSScovSite <- function(y, n, psi=~1, p=~1, data=NULL) {
  # single-season occupancy models with site-specific covatiates
  # new version with y/n input; much faster!
  # n is a vector with the number of occasions at each site.
  # y is a vector with the number of detections at each site.
  # ** psi and p are one-sided formulae describing the model
  # ** data a data frame with the site covariates.

  if(length(n) == 1)
    n <- rep(n, length(y))
  if(length(y) != length(n))
    stop("y and n must have the same length")
  nSites <- length(y)
  # nSurv <- ncol(DH)
  # if(is.null(rownames(DH)))
    # rownames(DH) <- 1:nSites
  # Sanity checks:
  # if (nSurv < 2)
    # stop("More than one survey occasion is needed")
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

  psiModMat <- model.matrix(psi, ddf)
  psiK <- ncol(psiModMat)
  pModMat <- model.matrix(p, ddf)
  pK <- ncol(pModMat)
  K <- psiK + pK

  beta.mat <- matrix(NA_real_, K, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c(
    paste("psi:", colnames(psiModMat)),
    paste("p:", colnames(pModMat)))
  lp.mat <- matrix(NA_real_, nSites*2, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c(
    paste("psi:", 1:nSites, sep=""),
    paste("p:", 1:nSites, sep=""))
  AIC <- NA_real_

  nll <- function(param){
    psiBeta <- param[1:psiK]
    pBeta <- param[(psiK+1):K]
    psiProb <- plogis(psiModMat %*% psiBeta)
    pProb <- plogis(pModMat %*% pBeta)
    prob <- psiProb * pProb^y * (1-pProb)^(n - y) + (1 - psiProb) * (y==0)
    return(min(-sum(log(prob)), .Machine$double.xmax))
    # p.DH <- sweep(DH, 1, pProb, "*") + sweep((1-DH), 1, (1-pProb), "*")
    # llh <- sum(log(psiProb * apply(p.DH, 1, prod, na.rm=TRUE) + 
          # (1 - psiProb) * (rowSums(DH, na.rm=TRUE) == 0)))
    # return(min(-llh, .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  res <- nlm(nll, param, hessian=TRUE)
  if(res$code < 3)  {  # exit code 1 or 2 is ok.
    beta.mat[,1] <- res$estimate
    lp.mat[, 1] <- c(psiModMat %*% beta.mat[1:psiK, 1],
                     pModMat %*% beta.mat[(psiK+1):K, 1])
    AIC <- 2*res$minimum + 2 * K
    if (det(res$hessian) > 0) {
      varcov <- solve(res$hessian)
      SE <- sqrt(diag(varcov))
      beta.mat[, 2] <- SE
      crit <- qnorm(c(0.025, 0.975))
      beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      SElp <- c(sqrt(diag(psiModMat %*% varcov[1:psiK, 1:psiK] %*% t(psiModMat))),
                sqrt(diag(pModMat %*% varcov[(psiK+1):K, (psiK+1):K] %*% t(pModMat))))
      lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              real = plogis(lp.mat),
              AIC = AIC)
  class(out) <- c("occupancy", "list")
  return(out)
}


