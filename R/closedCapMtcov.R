# Mtcov model, capture probability a function of time dependent covariates

closedCapMtcov <-
function(CH, p=~1, data=NULL, ci = 0.95, ciType=c("normal", "MARK")) {
  # CH is a 1/0 capture history matrix, animals x occasions
  # ci is the required confidence interval
  
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))
  ciType <- match.arg(ciType)

  CH <- round(as.matrix(CH))
  nocc <- ncol(CH)    # number of capture occasions
  N.cap <- nrow(CH)   # total number of individual animals captured
  n <- colSums(CH)    # vector of number of captures on each occasion

  if(!is.null(data)) {
    stopifnot(nrow(data) == nocc)
    stopifnot(sum(is.na(data)) == 0)
    data <- lapply(data, function(x) if(is.numeric(x)) scale(x) else x)
    ddf <- as.data.frame(data)
  } else {
    ddf <- data.frame(.dummy = rep(NA, nocc))
    # model.matrix needs a data frame, NULL won't do.
  }
  pModMat <- model.matrix(as.formula(p), ddf)
  K <- ncol(pModMat)

  beta.mat <- matrix(NA_real_, K+1, 4) # objects to hold output
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("Nhat", colnames(pModMat))
  lp.mat <- matrix(NA_real_, nocc, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- paste0("p", 1:nocc)
  logLik <- NA_real_
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
    if(res$code < 3)  {  # exit code 1 or 2 is ok.
      beta.mat[,1] <- res$estimate
      lp.mat[, 1] <- pModMat %*% beta.mat[-1, 1]
      varcov <- try(solve(res$hessian), silent=TRUE)
      if (!inherits(varcov, "try-error") && all(diag(varcov) > 0)) {
        beta.mat[, 2] <- sqrt(diag(varcov))
        beta.mat[, 3:4] <- sweep(outer(beta.mat[, 2], crit), 1, res$estimate, "+")
        temp <- diag(pModMat %*% varcov[-1, -1] %*% t(pModMat))
        if(all(temp >= 0))  {
          SElp <- sqrt(temp)
          lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
          logLik <- -res$minimum
        }
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
          real = rbind(Nhat, plogis(lp.mat)),
          logLik = c(logLik=logLik, df=K+1, nobs=length(CH)))
  class(out) <- c("closedCap", "list")
  return(out)
}
