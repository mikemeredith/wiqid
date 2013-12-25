# Plot added 2013-12-02
# occSS.T2 added 2013-12-02

occSS.T <-
function(DH, ci=0.95, plot=TRUE)  {
  # DH is a 1/0 matrix of detection histories, sites x occasions
  # ci is the required confidence interval.
  # This returns estimate & CI of psi as usual, plus the intercept
  #   and slope of the p vs time relationship on the logistic scale.

  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  DH <- as.matrix(DH)  # in case it's a data frame
	nocc <- ncol(DH)
  Time <- scale(0:(nocc-1))
	n.par <- 3
  beta.mat <- matrix(NA_real_, n.par, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("psi", "p:(Intercept)", "p:Time")
  lp.mat <- matrix(NA_real_, nocc+1, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c("psi", paste("p", 1:nocc, sep=""))
  logLik <- NA_real_
  varcov <- NULL

  if(ncol(DH) > 1 && sum(DH, na.rm=TRUE) > 0)  {
    # Negative log-likelihood function:
    nll <- function(params) {
      psi <- plogis(params[1])
      p <- plogis(params[2] + params[3] * Time)
      p.DH <- sweep(DH, 2, p, "*") + sweep((1-DH), 2, (1-p), "*")
      llh <- sum(log(psi * apply(p.DH, 1, prod, na.rm=TRUE) + 
          (1 - psi) * (rowSums(DH, na.rm=TRUE) == 0)))
      return(min(-llh, .Machine$double.xmax)) # min(..) needed to stop Inf being returned
    }
    params <- rep(0, n.par)
    res <- nlm(nll, params, hessian=TRUE)
    if(res$code < 3)  {  # exit code 1 or 2 is ok.
      beta.mat[,1] <- res$estimate
      varcov <- try(solve(res$hessian), silent=TRUE)
      if (!inherits(varcov, "try-error") && all(diag(varcov) > 0)) {
        SE <- sqrt(diag(varcov))
        beta.mat[, 'SE'] <- sqrt(diag(varcov))
        beta.mat[, 3:4] <- sweep(outer(beta.mat[, 'SE'], crit), 1, res$estimate, "+")
        logLik <- -res$minimum
      } else {
        varcov <- NULL
      }
    }
    if(mean(DH, na.rm=TRUE) == 1) {
      beta.mat[3:4 ] <- NA
      logLik <- NA_real_
    }
    lp.mat[1, ]  <- beta.mat[1, -2]
    lp.mat[-1, 1]  <- beta.mat[2, 1] + beta.mat[3, 1] * Time
    if(!is.null(varcov)) {
      pModMat <- cbind(1, Time)
      SElp <- sqrt(diag(pModMat %*% varcov[-1, -1] %*% t(pModMat)))
      lp.mat[-1, -1] <- sweep(outer(SElp, crit), 1, lp.mat[-1, 1], "+")
    }
    # Do the plot
    if(plot) {
      real.p <- plogis(lp.mat[-1, ])
      ylim <- range(0, real.p, na.rm=TRUE)
      plot(1:nocc, real.p[, 1], type='l', ylim=ylim,
        xlab="Time", ylab="Probability of detection")
      lines(1:nocc, real.p[, 2], lty=3)
      lines(1:nocc, real.p[, 3], lty=3)
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=n.par, nobs=nrow(DH)))
  class(out) <- c("wiqid", "list")
  return(out)
}

# ============================================================

# This does a quadratic model, p ~ Time + Time2

occSS.T2 <-
function(DH, ci=0.95, plot=TRUE)  {
  # DH is a 1/0 matrix of detection histories, sites x occasions
  # ci is the required confidence interval.
  # This returns estimate & CI of psi as usual, plus the intercept
  #   and slope of the p vs time relationship on the logistic scale.

  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  DH <- as.matrix(DH)  # in case it's a data frame
	nocc <- ncol(DH)
  Time <- scale(0:(nocc-1))
  Time2 <- Time^2
	n.par <- 4
  beta.mat <- matrix(NA_real_, n.par, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("psi", "p:(Intercept)", "p:Time", "p:Time2")
  lp.mat <- matrix(NA_real_, nocc+1, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c("psi", paste("p", 1:nocc, sep=""))
  logLik <- NA_real_
  varcov <- NULL
  
  if(ncol(DH) > 1 && sum(DH, na.rm=TRUE) > 0)  {
    # Negative log-likelihood function:
    nll <- function(params) {
      psi <- plogis(params[1])
      p <- plogis(params[2] + params[3] * Time + params[4] * Time2)
      p.DH <- sweep(DH, 2, p, "*") + sweep((1-DH), 2, (1-p), "*")
      llh <- sum(log(psi * apply(p.DH, 1, prod, na.rm=TRUE) + 
          (1 - psi) * (rowSums(DH, na.rm=TRUE) == 0)))
      return(min(-llh, .Machine$double.xmax)) # min(..) needed to stop Inf being returned
    }
    params <- rep(0, n.par)
    res <- nlm(nll, params, hessian=TRUE)
    if(res$code < 3)  {  # exit code 1 or 2 is ok.
      beta.mat[,1] <- res$estimate
      logLik <- -res$minimum
      varcov <- try(solve(res$hessian), silent=TRUE)
      if (!inherits(varcov, "try-error") && all(diag(varcov) > 0)) {
        SE <- sqrt(diag(varcov))
        beta.mat[, 'SE'] <- sqrt(diag(varcov))
        beta.mat[, 3:4] <- sweep(outer(beta.mat[, 'SE'], crit), 1, res$estimate, "+")
      } else {
        varcov <- NULL
      }
    }
    if(mean(DH, na.rm=TRUE) == 1) {
      beta.mat[3:4 ] <- NA
      logLik <- NA_real_
    }
    lp.mat[1, ]  <- beta.mat[1, -2]
    lp.mat[-1, 1]  <- beta.mat[2, 1] + beta.mat[3, 1] * Time+ beta.mat[4, 1] * Time2
    if(!is.null(varcov)) {
      pModMat <- cbind(1, Time, Time2)
      SElp <- sqrt(diag(pModMat %*% varcov[-1, -1] %*% t(pModMat)))
      lp.mat[-1, -1] <- sweep(outer(SElp, crit), 1, lp.mat[-1, 1], "+")
    }
    # Do the plot
    if(plot) {
      real.p <- plogis(lp.mat[-1, ])
      ylim <- range(0, real.p, na.rm=TRUE)
      plot(1:nocc, real.p[, 1], type='l', ylim=ylim,
        xlab="Time", ylab="Probability of detection")
      lines(1:nocc, real.p[, 2], lty=3)
      lines(1:nocc, real.p[, 3], lty=3)
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=n.par, nobs=nrow(DH)))
  class(out) <- c("wiqid", "list")
  return(out)
}
