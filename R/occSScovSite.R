# Single season occupancy with site covariates (not survey covariates)

occSScovSite <- function(y, n, psi=~1, p=~1, data=NULL) {
  # single-season occupancy models with site-specific covatiates
  # new version with y/n input; much faster!
  # y is a vector with the number of detections at each site.
  # n is a vector with the number of occasions at each site.
  # ** psi and p are one-sided formulae describing the model
  # ** data a data frame or list with the site covariates.

  if(length(n) == 1)
    n <- rep(n, length(y))
  if(length(y) != length(n))
    stop("y and n must have the same length")
  if(any(y > n))
    stop("y cannot be greater than n")
  nSites <- length(y)
  if(!is.null(data))  {
    data <- lapply(data, function(x) if(is.numeric(x)) scale(x) else x)
    ddf <- as.data.frame(data)
    if(any(is.na(ddf)))
      stop("Missing site covariates are not allowed.")
    if(nrow(ddf) != nSites)
      stop("'data' must have a row for each site.")
  } else {
    ddf <- data.frame(.dummy = rep(NA, nSites))
    # model.matrix needs a data frame, NULL won't do.
  }

  psiModMat <- model.matrix(as.formula(psi), ddf)
  psiK <- ncol(psiModMat)
  pModMat <- model.matrix(as.formula(p), ddf)
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
  logLik <- NA_real_

  nll <- function(param){
    psiBeta <- param[1:psiK]
    pBeta <- param[(psiK+1):K]
    psiProb <- as.vector(plogis(psiModMat %*% psiBeta))
    pProb <- as.vector(plogis(pModMat %*% pBeta))
    prob <- psiProb * pProb^y * (1-pProb)^(n - y) + (1 - psiProb) * (y==0)
    return(min(-sum(log(prob)), .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  res <- nlm(nll, param, hessian=TRUE)
  if(res$code < 4)  {  # exit code 1 or 2 is ok, 3 doubtful.
    beta.mat[,1] <- res$estimate
    lp.mat[, 1] <- c(psiModMat %*% beta.mat[1:psiK, 1],
                     pModMat %*% beta.mat[(psiK+1):K, 1])
    varcov <- try(solve(res$hessian), silent=TRUE)
    if (!inherits(varcov, "try-error") && 
        all(diag(varcov) > 0)) {
      SE <- sqrt(diag(varcov))
      beta.mat[, 2] <- SE
      crit <- qnorm(c(0.025, 0.975))
      beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      temp <- c(diag(psiModMat %*% varcov[1:psiK, 1:psiK] %*% t(psiModMat)),
                diag(pModMat %*% varcov[(psiK+1):K, (psiK+1):K] %*% t(pModMat)))
      if(all(temp >= 0))  {
        SElp <- sqrt(temp)
        lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
        logLik <- -res$minimum
      }
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=length(y)))
  class(out) <- c("occupancy", "list")
  return(out)
}


