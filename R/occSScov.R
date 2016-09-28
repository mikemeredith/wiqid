# Single season occupancy with site and survey covariates.

# 'link' argument added 2015-02-20

occSS <- function(DH, model=NULL, data=NULL, ci=0.95, link=c("logit", "probit"),
  verify=TRUE) {
  # single-season occupancy models with site and survey covatiates
  # ** DH is detection data in a 1/0/NA matrix or data frame, sites in rows,
  #    detection occasions in columns..
  # ** model is a list of 2-sided formulae for psi and p; can also be a single
  #   2-sided formula, eg, model = psi ~ habitat.
  # ** data is a DATA FRAME with single columns for site covariates and a column for each survey occasion for each survey covariate.
  # ci is the required confidence interval.
  
  if(verify)
    DH <- verifyDH(DH, allowNA=TRUE)

  if(is.null(model)) {
    y <- rowSums(DH, na.rm=TRUE)
    n <- rowSums(!is.na(DH))
    return(occSS0(y, n, ci=ci, link=link))
  }
  crit <- fixCI(ci)

  if(match.arg(link) == "logit") {
    plink <- plogis
  } else {
    plink <- pnorm
  }

  # Standardise the model:
  model <- stdModel(model, list(psi=~1, p=~1))

  # Summarize detection history
  site.names <- rownames(DH)
  DH <- as.matrix(DH)
  nSites <- nrow(DH)
  nSurv <- ncol(DH)
  if (nSurv < 2)
    stop("More than one survey occasion is needed")
  if(is.null(site.names))
    site.names <- 1:nSites

  # Convert the covariate data frame into a list
  dataList <- stddata(data, nSurv)
  time <- rep(1:nSurv, each=nSites)
  dataList$.Time <- as.vector(scale(time)) /2
  dataList$.time <- as.factor(time)
  before <- cbind(0L, DH[, 1:(nSurv - 1)]) # 1 if animal seen on previous occasion
  dataList$.b <- as.vector(before)

  survey.done <- !is.na(as.vector(DH))
  DHvec <- as.vector(DH)[survey.done]
  siteID <- as.factor(row(DH))[survey.done]
  survID <- as.factor(col(DH))[survey.done]

  psiDf <- selectCovars(model$psi, dataList, nSites)
  if (nrow(psiDf) != nSites)
    stop("Number of site covars doesn't match sites.\nAre you using survey covars?")
  psiModMat <- model.matrix(model$psi, psiDf)
  if(nrow(psiModMat) != nrow(psiDf))
      stop("Missing site covariates are not allowed.")
  psiK <- ncol(psiModMat)
  pDf0 <- selectCovars(model$p, dataList, nSites*nSurv)
  pDf <- pDf0[survey.done, , drop=FALSE]
  pModMat <- model.matrix(model$p, pDf)
  if(nrow(pModMat) != nrow(pDf))
      stop("Missing survey covariates are not allowed when a survey was done.")
  pK <- ncol(pModMat)
  K <- psiK + pK

  # objects to hold output
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
  varcov <- NULL

  # Negative log likelihood function
  nll <- function(param){
    psiBeta <- param[1:psiK]
    pBeta <- param[(psiK+1):K]
    psiProb <- as.vector(plink(psiModMat %*% psiBeta))
    pProb <- plink(pModMat %*% pBeta)
    Lik1 <- DHvec*pProb + (1-DHvec) * (1-pProb)
    Lik2 <- tapply(Lik1, siteID, prod)
    llh <- sum(log(psiProb * Lik2 +
          (1 - psiProb) * (rowSums(DH, na.rm=TRUE) == 0)))
    return(min(-llh, .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  res <- nlm(nll, param, hessian=TRUE)
  if(res$code > 2)   # exit code 1 or 2 is ok.
    warning(paste("Convergence may not have been reached (code", res$code, ")"))
  beta.mat[,1] <- res$estimate
  lp.mat[, 1] <- c(psiModMat %*% beta.mat[1:psiK, 1],
                   pModMat %*% beta.mat[(psiK+1):K, 1])
  if(res$code < 3) # Keep NA if in doubt
    logLik <- -res$minimum
  # varcov0 <- try(solve(res$hessian), silent=TRUE)
  varcov0 <- try(chol2inv(chol(res$hessian)), silent=TRUE)
  # if (!inherits(varcov0, "try-error") && all(diag(varcov0) > 0)) {
  if (!inherits(varcov0, "try-error")) {
    varcov <- varcov0
    SE <- suppressWarnings(sqrt(diag(varcov)))
    beta.mat[, 2] <- SE
    beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
    SElp <- c(sqrt(diag(psiModMat %*% varcov[1:psiK, 1:psiK] %*% t(psiModMat))),
              sqrt(diag(pModMat %*% varcov[(psiK+1):K, (psiK+1):K] %*% t(pModMat))))
    lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
  }
  out <- list(call = match.call(),
              link = match.arg(link),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plink(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=nrow(DH)))
  class(out) <- c("wiqid", "list")
  return(out)
}


