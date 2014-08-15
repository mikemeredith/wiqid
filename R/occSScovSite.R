# Single season occupancy with site covariates (not survey covariates)

# 'model' argument added 2013-12-02

occSScovSite <- function(y, n, model=NULL, data=NULL, ci=0.95) {
  # single-season occupancy models with site-specific covatiates
  # new version with y/n input; much faster!
  # y is a vector with the number of detections at each site.
  # n is a vector with the number of occasions at each site.
  # model is a list of 2-sided formulae for psi and p; can also be a single
  #   2-sided formula, eg, model = psi ~ habitat.
  # ci is the required confidence interval.
  if(length(n) == 1)
    n <- rep(n, length(y))
  if(length(y) != length(n))
    stop("y and n must have the same length")
  if(any(y > n))
    stop("y cannot be greater than n")

  crit <- fixCI(ci)
  
  # Standardise the model:
  model <- stdModel(model, list(psi=~1, p=~1))

  # Convert the covariate data frame into a list
  nSites <- length(y)
  dataList <- stddata(data, nocc=NULL)

  psiDf <- selectCovars(model$psi, dataList, nSites)
  if (nrow(psiDf) != nSites)
    stop("Number of site covars doesn't match sites.")
  psiModMat <- model.matrix(model$psi, psiDf)
  psiK <- ncol(psiModMat)
  pDf <- selectCovars(model$p, dataList, nSites)
  if (nrow(pDf) != nSites)
    stop("Number of site covars doesn't match sites.")
  pModMat <- model.matrix(model$p, pDf)
  pK <- ncol(pModMat)
  K <- psiK + pK
  # model.matrix removes rows with NAs:
  if(nrow(psiModMat) != nSites || nrow(pModMat) != nSites)
    stop("Missing site covariates are not allowed.")

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
  if(res$code > 2)   # exit code 1 or 2 is ok.
    warning(paste("Convergence may not have been reached (code", res$code, ")"))

  # Process output
  beta.mat[,1] <- res$estimate
  lp.mat[, 1] <- c(psiModMat %*% beta.mat[1:psiK, 1],
                   pModMat %*% beta.mat[(psiK+1):K, 1])
  varcov <- try(solve(res$hessian), silent=TRUE)
  if (!inherits(varcov, "try-error") && 
      all(diag(varcov) > 0)) {
    SE <- sqrt(diag(varcov))
    beta.mat[, 2] <- SE
    beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
    temp <- c(diag(psiModMat %*% varcov[1:psiK, 1:psiK] %*% t(psiModMat)),
              diag(pModMat %*% varcov[(psiK+1):K, (psiK+1):K] %*% t(pModMat)))
    if(all(temp >= 0))  {
      SElp <- sqrt(temp)
      lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
      logLik <- -res$minimum
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=length(y)))
  class(out) <- c("wiqid", "list")
  return(out)
}


