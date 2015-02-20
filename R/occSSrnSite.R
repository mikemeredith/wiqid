# Single-season Royle-Nichols occupancy models, with abundance-induced
#   heterogeneity in detection probability, with site covariates
#   (not survey covariates)

# 'link' argument added 2015-02-20

occSSrnSite <- function(y, n, model=NULL, data=NULL,
    ci=0.95, link=c("logit", "probit")) {
  # single-season occupancy models with site-specific covatiates
  # y is a vector with the number of detections at each site.
  # n is a vector with the number of occasions at each site.
  # model is a list of 2-sided formulae for lambda and r; can also be a single
  #   2-sided formula, eg, model = lambda ~ habitat.
  # ci is the required confidence interval.
  if(length(n) == 1)
    n <- rep(n, length(y))
  if(length(y) != length(n))
    stop("y and n must have the same length")
  if(any(y > n))
    stop("y cannot be greater than n")
  crit <- fixCI(ci)

  if(match.arg(link) == "logit") {
    plink <- plogis
  } else {
    plink <- pnorm
  }

  # Standardise the model:
  model <- stdModel(model, list(lambda=~1, r=~1))

  # Convert the covariate data frame into a list
  nSites <- length(y)
  dataList <- stddata(data, nocc=NULL)

  lamDf <- selectCovars(model$lambda, dataList, nSites)
  if (nrow(lamDf) != nSites)
    stop("Number of site covars doesn't match sites.")
  lamModMat <- model.matrix(model$lambda, lamDf)
  lamK <- ncol(lamModMat)
  rDf <- selectCovars(model$r, dataList, nSites)
  if (nrow(rDf) != nSites)
    stop("Number of site covars doesn't match sites.")
  rModMat <- model.matrix(model$r, rDf)
  rK <- ncol(rModMat)
  K <- lamK + rK
  # model.matrix removes rows with NAs:
  if(nrow(lamModMat) != nSites || nrow(rModMat) != nSites)
    stop("Missing site covariates are not allowed.")

  beta.mat <- matrix(NA_real_, K, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c(
    paste("lambda:", colnames(lamModMat)),
    paste("r:", colnames(rModMat)))
  lp.mat <- matrix(NA_real_, nSites*2, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c(
    paste("lambda:", 1:nSites, sep=""),
    paste("r:", 1:nSites, sep=""))
  logLik <- NA_real_
  varcov <- NULL

  nll <- function(param){
    lamBeta <- param[1:lamK]
    rBeta <- param[(lamK+1):K]
    lambda <- as.vector(exp(lamModMat %*% lamBeta))
    s <- 1 - as.vector(plink(rModMat %*% rBeta)) # s = 1 - r
    llh <- numeric(nSites)
    for(i in 1:nSites) {
      rpart <- s[i]^(0:Nmax)
      Npart <- dpois(0:Nmax, lambda[i])
      llh[i] <- log(sum((1-rpart)^y[i] * rpart^(n[i]-y[i]) * Npart))
    }
    return(min(-sum(llh), .Machine$double.xmax))
  }

  # Run mle estimation with nlm:
  param <- rep(0, K)
  Nmax <- 100
  res <- nlm(nll, param, hessian=TRUE)
  if(res$code > 2)   # exit code 1 or 2 is ok.
    warning(paste("Convergence may not have been reached (code", res$code, ")"))

  # Process output
  beta.mat[,1] <- res$estimate
  lp.mat[, 1] <- c(lamModMat %*% beta.mat[1:lamK, 1],
                   rModMat %*% beta.mat[(lamK+1):K, 1])
  varcov0 <- try(solve(res$hessian), silent=TRUE)
  if (!inherits(varcov0, "try-error") && all(diag(varcov0) > 0)) {
    varcov <- varcov0
    SE <- sqrt(diag(varcov))
    beta.mat[, 2] <- SE
    beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
    temp <- c(diag(lamModMat %*% varcov[1:lamK, 1:lamK] %*% t(lamModMat)),
              diag(rModMat %*% varcov[(lamK+1):K, (lamK+1):K] %*% t(rModMat)))
    if(all(temp >= 0))  {
      SElp <- sqrt(temp)
      lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
      logLik <- -res$minimum
    }
  }
  realLam <- exp(lp.mat[1:nSites, ])
  realR <- plink(lp.mat[(nSites+1):(nSites*2), ])
  realPsi <- 1-dpois(0, realLam)
  rownames(realPsi) <- paste("psi:", 1:nSites, sep="")

  out <- list(call = match.call(),
              link = match.arg(link),
              beta = beta.mat,
              beta.vcv = varcov,
              real = rbind(realPsi, realLam, realR),
              logLik = c(logLik=logLik, df=K, nobs=length(y)))
  class(out) <- c("wiqid", "list")
  return(out)
}


