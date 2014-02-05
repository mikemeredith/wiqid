closedCapMh2 <-
function(CH, ci = 0.95, ciType=c("normal", "MARK")) {
  # CH is a 1/0 capture history matrix, animals x occasions, OR
  #  a vector of capture frequencies of length equal to the number
  #  of occasions - trailing zeros are required.
  # ci is the required confidence interval
  
  if (is.matrix(CH) || is.data.frame(CH)) {
    n.occ <- ncol(CH)
    freq <- tabulate(rowSums(CH), nbins=n.occ)
  } else {
    freq <- round(CH)
    n.occ <- length(freq)
  }

  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))
  ciType <- match.arg(ciType)

  N.cap <- sum(freq)  # Number of individual animals captured
  beta.mat <- matrix(NA_real_, 4, 4) # objects to hold output
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("Nhat", "piHat", "p1hat", "p2hat")
  logLik <- NA_real_
  varcov <- NULL
  
  if(sum(freq[-1]) > 1)  {  # Do checks here
    # nll1 ensures p2 <= p1, but useless for Hessian
    # nll2 starts with output from nll1 (so p2/p1 not an issue)
    nll1 <- function(params) {
      f0 <- min(exp(params[1]), 1e+300, .Machine$double.xmax)
      f.vect <- c(f0, freq)
      pi <- plogis(params[2])
      p1 <- plogis(params[3])
      p2 <- p1 * plogis(params[4]) # ensures that p2 <= p1
      bin.co <- lgamma(N.cap + f0 + 1) - lgamma(f0 + 1)
      tmp <- numeric(length(f.vect))
      for(i in seq_along(f.vect))
         tmp[i] <- pi * p1^(i-1) * (1-p1)^(n.occ-i+1) +
                  (1-pi) * p2^(i-1) * (1-p2)^(n.occ-i+1)
      llh <- bin.co + sum(f.vect * log(tmp))
      return(min(-llh, .Machine$double.xmax))
    }

    nll2 <- function(params) {
      f0 <- min(exp(params[1]), 1e+300, .Machine$double.xmax)
      f.vect <- c(f0, freq)
      pi <- plogis(params[2])
      p1 <- plogis(params[3])
      p2 <- plogis(params[4])
      bin.co <- lgamma(N.cap + f0 + 1) - lgamma(f0 + 1)
      tmp <- numeric(length(f.vect))
      for(i in seq_along(f.vect))
         tmp[i] <- pi * p1^(i-1) * (1-p1)^(n.occ-i+1) +
                  (1-pi) * p2^(i-1) * (1-p2)^(n.occ-i+1)
      llh <- bin.co + sum(f.vect * log(tmp))
      return(min(-llh, .Machine$double.xmax))
    }

    res1 <- suppressWarnings(nlm(nll1, c(log(5), 0,0,0)))
    params <- res1$estimate
    p2 <- plogis(params[3]) * plogis(params[4])
    params[4] <- qlogis(p2)
    res2 <- nlm(nll2, params, hessian=TRUE)
    if(res2$code > 2)   # exit code 1 or 2 is ok.
      warning(paste("Convergence may not have been reached (code", res2$code, ")"))
    beta.mat[, 1] <- res2$estimate
    varcov <- try(solve(res2$hessian), silent=TRUE)
    if (!inherits(varcov, "try-error") && all(diag(varcov) > 0)) {
      beta.mat[, 2] <- sqrt(diag(varcov))
      beta.mat[, 3:4] <- sweep(outer(beta.mat[, 2], crit), 1, res2$estimate, "+")
      logLik <- -res2$minimum
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
          logLik = c(logLik=logLik, df=4, nobs=N.cap * n.occ))
  class(out) <- c("wiqid", "list")
  return(out)
}
