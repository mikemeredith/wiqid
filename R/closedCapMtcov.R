# Mtcov model, capture probability a function of time dependent covariates

closedCapMtcov <-
function(CH, model=list(p~1), data=NULL, ci = 0.95, ciType=c("normal", "MARK")) {
  # CH is a 1/0 capture history matrix, animals x occasions
  # ci is the required confidence interval
  
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))
  ciType <- match.arg(ciType)

  # Standardise the model:
  model <- stdModel(model, list(p=~1))

  CH <- round(as.matrix(CH))
  nocc <- ncol(CH)    # number of capture occasions
  N.cap <- nrow(CH)   # total number of individual animals captured
  n <- colSums(CH)    # vector of number of captures on each occasion

  # Convert the covariate data frame into a list
  dataList <- stddata(data, nocc)
  dataList$.Time <- scale(1:nocc)
  dataList$.time <- as.factor(1:nocc)
  ddf <- as.data.frame(dataList)

  pModMat <- model.matrix(model$p, ddf)
  K <- ncol(pModMat)

  beta.mat <- matrix(NA_real_, K+1, 4) # objects to hold output
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("Nhat", colnames(pModMat))
  lp.mat <- matrix(NA_real_, nocc, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- paste0("p", 1:nocc)
  logLik <- NA_real_
  varcov <- NULL
  
  if(N.cap > 0)  {
    nll <- function(params) {
      N <- min(exp(params[1]) + N.cap, 1e+300, .Machine$double.xmax)
      pBeta <- params[-1]
      p <- as.vector(plogis(pModMat %*% pBeta))
      tmp <- lgamma(N + 1) - lgamma(N - N.cap + 1) + 
        sum(n * log(p) + (N - n) * log(1-p))
      return(min(-tmp, .Machine$double.xmax))
    }
    params <- c(log(5), rep(0, K))
    res <- nlm(nll, params, hessian=TRUE, iterlim=1000)
    if(res$code > 2)   # exit code 1 or 2 is ok.
      warning(paste("Convergence may not have been reached (code", res$code, ")"))
    beta.mat[,1] <- res$estimate
    lp.mat[, 1] <- pModMat %*% beta.mat[-1, 1]
    # varcov0 <- try(solve(res$hessian), silent=TRUE)
    varcov0 <- try(chol2inv(chol(res$hessian)), silent=TRUE)
    # if (!inherits(varcov0, "try-error") && all(diag(varcov0) > 0)) {
    if (!inherits(varcov0, "try-error")) {
      varcov <- varcov0
      beta.mat[, 2] <- suppressWarnings(sqrt(diag(varcov)))
      beta.mat[, 3:4] <- sweep(outer(beta.mat[, 2], crit), 1, res$estimate, "+")
      temp <- diag(pModMat %*% varcov[-1, -1] %*% t(pModMat))
      if(all(temp >= 0))  {
        SElp <- sqrt(temp)
        lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
        logLik <- -res$minimum
      }
    }
  }
  if(ciType == "normal") {
    Nhat <- exp(beta.mat[1, -2]) + N.cap
  } else {
    Nhat <- getMARKci(beta.mat[1, 1], beta.mat[1, 2], ci) + N.cap
  }
  out <- list(call = match.call(),
          beta = beta.mat,
          beta.vcv = varcov,
          real = rbind(Nhat, plogis(lp.mat)),
          logLik = c(logLik=logLik, df=K+1, nobs=length(CH)))
  class(out) <- c("wiqid", "list")
  return(out)
}
