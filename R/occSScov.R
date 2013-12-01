# Single season occupancy with site and survey covariates.

occSScov <- function(DH, psi=~1, p=~1, data=NULL, ci=0.95) {
  # single-season occupancy models with site and survey covatiates
  # ** DH is detection data in a 1/0/NA matrix or data frame, sites in rows, 
  #    detection occasions in columns..
  # ** psi and p are one-sided formulae describing the model
  # ** data is a LIST with vectors for site covariates and site x survey matrices
  #     for survey covariates.
  # ci is the required confidence interval.
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  site.names <- rownames(DH)
  DH <- as.matrix(DH)
  nSites <- nrow(DH)
  nSurv <- ncol(DH)
  if(is.null(site.names))
    site.names <- 1:nSites
  # Sanity checks:
  if (nSurv < 2)
    stop("More than one survey occasion is needed")
  survey.done <- !is.na(as.vector(DH))
  DHvec <- as.vector(DH)[survey.done]
  siteID <- as.factor(row(DH))[survey.done]
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

  psiModMat <- model.matrix(as.formula(psi), psiDf)
  psiK <- ncol(psiModMat)
  pModMat <- model.matrix(as.formula(p), pDf)
  pK <- ncol(pModMat)
  K <- psiK + pK

  beta.mat <- matrix(NA_real_, K, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c(
    paste("psi:", colnames(psiModMat)),
    paste("p:", colnames(pModMat)))
  lp.mat <- matrix(NA_real_, nSites + sum(survey.done), 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c(
    paste("psi:", site.names, sep=""),
    paste("p:", siteID, ",", survID, sep=""))
  logLik <- NA_real_

  nll <- function(param){
    psiBeta <- param[1:psiK]
    pBeta <- param[(psiK+1):K]
    psiProb <- as.vector(plogis(psiModMat %*% psiBeta))
    pProb <- plogis(pModMat %*% pBeta)
    Lik1 <- DHvec*pProb + (1-DHvec) * (1-pProb)
    Lik2 <- tapply(Lik1, siteID, prod)
    llh <- sum(log(psiProb * Lik2 + 
          (1 - psiProb) * (rowSums(DH, na.rm=TRUE) == 0)))
    return(min(-llh, .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  res <- nlm(nll, param, hessian=TRUE)
  if(res$code < 4)  {  # exit code 1 or 2 is ok, 3 dodgy but...
    beta.mat[,1] <- res$estimate
    lp.mat[, 1] <- c(psiModMat %*% beta.mat[1:psiK, 1],
                     pModMat %*% beta.mat[(psiK+1):K, 1])
    if(res$code < 3) # Keep NA if in doubt
      logLik <- -res$minimum
    varcov <- try(solve(res$hessian), silent=TRUE)
    if (!inherits(varcov, "try-error") &&
        all(diag(varcov) > 0)) {
      SE <- sqrt(diag(varcov))
      beta.mat[, 2] <- SE
      beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      SElp <- c(sqrt(diag(psiModMat %*% varcov[1:psiK, 1:psiK] %*% t(psiModMat))),
                sqrt(diag(pModMat %*% varcov[(psiK+1):K, (psiK+1):K] %*% t(pModMat))))
      lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=nrow(DH)))
  class(out) <- c("occupancy", "list")
  return(out)
}


