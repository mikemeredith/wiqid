
# Royle-Nichols occupancy models, with abundance-induced heterogeneity in detection probability.

## WORK IN PROGRESS ##

occSSrn <- function(DH, model=NULL, data=NULL, ci=0.95) {
  # single-season Royle-Nichols model with site and survey covariates
  # ** DH is detection data in a 1/0/NA matrix or data frame, sites in rows, 
  #    detection occasions in columns.
  # ** model is a list of 2-sided formulae for lambda and r; can also be a single
  #   2-sided formula, eg, model = lambda ~ habitat.
  # ** data is a DATA FRAME with single columns for site covariates and a column for each survey occasion for each survey covariate.
  # ci is the required confidence interval.
  
  if (!is.matrix(DH) && !is.data.frame(DH))
    stop("DH should be a detection history matrix (or data frame)")
  
  if(TRUE) {  # if(is.null(model)) {
    y <- rowSums(DH, na.rm=TRUE)
    n <- rowSums(!is.na(DH))
    if(is.null(model)) {
      return(occSSrn0(y, n, ci=ci))
    } else {
      return(occSSrnSite(y, n, model=model, data=data, ci=ci))
    }
  }
  
  crit <- fixCI(ci)

  # Standardise the model:
  model <- stdModel(model, list(lambda=~1, r=~1))

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
  siteDex <- row(DH)[survey.done]
  siteID <- as.factor(row(DH))[survey.done]
  survID <- as.factor(col(DH))[survey.done]

  lamDf <- selectCovars(model$lambda, dataList, nSites)
  if (nrow(lamDf) != nSites)
    stop("Number of site covars doesn't match sites.\nAre you using survey covars?")
  lamModMat <- model.matrix(model$lambda, lamDf)
  if(nrow(lamModMat) != nrow(lamDf))
      stop("Missing site covariates are not allowed.")
  lamK <- ncol(lamModMat)
  rDf0 <- selectCovars(model$r, dataList, nSites*nSurv)
  rDf <- rDf0[survey.done, , drop=FALSE]
  rModMat <- model.matrix(model$r, rDf)
  if(nrow(rModMat) != nrow(rDf))
      stop("Missing survey covariates are not allowed when a survey was done.")
  rK <- ncol(rModMat)
  K <- lamK + rK

  # objects to hold output
  beta.mat <- matrix(NA_real_, K, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c(
    paste("lambda:", colnames(lamModMat)),
    paste("r:", colnames(rModMat)))
  lp.mat <- matrix(NA_real_, nSites + sum(survey.done), 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c(
    paste("lambda:", site.names, sep=""),
    paste("r:", siteID, ",", survID, sep=""))
  logLik <- NA_real_
  varcov <- NULL

  # Negative log likelihood function
  nll <- function(param){
    lamBeta <- param[1:lamK]
    rBeta <- param[(lamK+1):K]
    lambda <- as.vector(exp(lamModMat %*% lamBeta))
    s <- 1 - as.vector(plogis(rModMat %*% rBeta))    # s = 1 - r
    llh <- numeric(length(DHvec))
    for(i in seq_along(llh)) {
      qN <- s[i]^(0:Nmax) # prob of detecting no one for N=0, 1, 2, 3...
      Npart <- dpois(0:Nmax, lambda[siteDex[i]]) # prob N = 0, 1, 2, 3...
      llh[i] <- log(sum((if(DHvec[i]==1) 1-qN else qN) * Npart))
    }
    return(min(-sum(llh), .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  Nmax <- 100 # See later if this is sensible

  res <- nlm(nll, param, hessian=TRUE)
  if(res$code > 2)   # exit code 1 or 2 is ok.
    warning(paste("Convergence may not have been reached (code", res$code, ")"))
  beta.mat[,1] <- res$estimate
  lp.mat[, 1] <- c(psiModMat %*% beta.mat[1:psiK, 1],
                   pModMat %*% beta.mat[(psiK+1):K, 1])
  if(res$code < 3) # Keep NA if in doubt
    logLik <- -res$minimum
  varcov0 <- try(solve(res$hessian), silent=TRUE)
  if (!inherits(varcov0, "try-error") && all(diag(varcov0) > 0)) {
    varcov <- varcov0
    SE <- sqrt(diag(varcov))
    beta.mat[, 2] <- SE
    beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
    SElp <- c(sqrt(diag(psiModMat %*% varcov[1:psiK, 1:psiK] %*% t(psiModMat))),
              sqrt(diag(pModMat %*% varcov[(psiK+1):K, (psiK+1):K] %*% t(pModMat))))
    lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=nrow(DH)))
  class(out) <- c("wiqid", "list")
  return(out)
}


# ------------------------------------------------------------------

occSSrn0 <-
function(y, n, ci=0.95)  {
  # Fast version without covariates.
  # y is a vector with the number of detections at each site.
  # n is a vector with the number of occasions at each site.
  # ci is the required confidence interval.
  
   crit <- fixCI(ci)

	# Starting values:
  beta.mat <- matrix(NA_real_, 2, 4) 
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("lambda", "r")
  logLik <- NA_real_
  varcov <- NULL
  
  if(sum(n) > 0 && sum(y) > 0 && any(y < n)) {    # If all n's are 0, no data available.
    params <- c(0, 0)
    Nmax <- 100 # See later if this is sensible
    # Negative log-likelihood function:
    nll <- function(params) {
      lambda <- exp(params[1])
      r <- plogis(params[2])
      rpart <- (1-r)^(0:Nmax)
      Npart <- dpois(0:Nmax, lambda)
      llh <- 0
      for(i in seq_along(n)) {
        llh <- llh + log(sum((1-rpart)^y[i] * rpart^(n[i]-y[i]) * Npart))
      }
      return(min(-llh, .Machine$double.xmax)) # min(..) stops Inf being returned
    }
    res <- nlm(nll, params, hessian=TRUE)
    if(res$code > 2)   # exit code 1 or 2 is ok.
      warning(paste("Convergence may not have been reached (code", res$code, ")"))
    beta.mat[,1] <- res$estimate
    varcov0 <- try(solve(res$hessian), silent=TRUE)
    if (!inherits(varcov0, "try-error") && all(diag(varcov0) > 0)) {
      varcov <- varcov0
      SE <- sqrt(diag(varcov))
      beta.mat[, 2] <- SE
      beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      logLik <- -res$minimum
    } 
  }
	lambda <- exp(beta.mat[1, -2])
	real <- rbind(1-dpois(0, lambda), lambda, plogis(beta.mat[2, -2]))
  colnames(real) <- c("est", "lowCI", "uppCI")
  rownames(real) <- c("psiHat", "lambdaHat", "rHat")
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = real,
              logLik = c(logLik=logLik, df=2, nobs=length(y)))
  class(out) <- c("wiqid", "list")
  return(out)
}
