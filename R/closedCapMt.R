# Mt model, capture probability time dependent

closedCapMt <-
function(CH, ci = 0.95, ciType=c("normal", "MARK")) {
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
  n <- colSums(CH)    # vector of number captures on each occasion

  beta.mat <- matrix(NA_real_, nocc+1, 4) # objects to hold output
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("Nhat", paste0("p", 1:nocc))
  logLik <- NA_real_
  
  if(N.cap > 0)  {
    nll <- function(params) {
      N <- min(exp(params[1]) + N.cap, 1e+300, .Machine$double.xmax)
      p <- plogis(params[-1])
      tmp <- lgamma(N + 1) - lgamma(N - N.cap + 1) + 
        sum(n * log(p) + (N - n) * log(1-p))
      return(min(-tmp, .Machine$double.xmax))
    }
    params <- c(log(5), rep(0, nocc))
    res <- nlm(nll, params, hessian=TRUE, iterlim=1000)
    if(res$code < 3)  {  # exit code 1 or 2 is ok.
      beta.mat[,1] <- res$estimate
      varcov <- try(solve(res$hessian), silent=TRUE)
      if (!inherits(varcov, "try-error") && all(diag(varcov) > 0)) {
        beta.mat[, 2] <- sqrt(diag(varcov))
        beta.mat[, 3:4] <- sweep(outer(beta.mat[, 2], crit), 1, res$estimate, "+")
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
          real = rbind(Nhat, plogis(beta.mat[-1, -2])),
          logLik = c(logLik=logLik, df=nocc+1, nobs=length(CH)))
  class(out) <- c("wiqid", "list")
  return(out)
}

