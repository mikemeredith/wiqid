occSS.T <-
function(DH, ci=0.95)  {
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
  Time <- 0:(nocc-1)
	n.par <- 3
  beta.mat <- matrix(NA_real_, n.par, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("psi", "p:(Intercept)", "p:Time")
  lp.mat <- matrix(NA_real_, nocc+1, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c("psi", paste("p", 1:nocc, sep=""))
  # AIC <- NA_real_
  logLik <- NA_real_
  varcovar <- NULL
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
      # AIC <- 2*res$minimum + 2 * n.par
      logLik <- -res$minimum
      if (det(res$hessian) > 1e-6) {
        varcovar <- solve(res$hessian)
        beta.mat[, 'SE'] <- sqrt(diag(varcovar))
        beta.mat[, 3:4] <- sweep(outer(beta.mat[, 'SE'], crit), 1, res$estimate, "+")
      }
    }
    if(mean(DH, na.rm=TRUE) == 1) {
      beta.mat[3:4 ] <- NA
      # AIC <- NA_real_
      logLik <- NA_real_
    }
    lp.mat[1, ]  <- beta.mat[1, -2]
    lp.mat[-1, 1]  <- beta.mat[2, 1] + beta.mat[3, 1] * Time
    if(!is.null(varcovar)) {
      pModMat <- cbind(1, Time)
      SElp <- sqrt(diag(pModMat %*% varcovar[-1, -1] %*% t(pModMat)))
      lp.mat[-1, -1] <- sweep(outer(SElp, crit), 1, lp.mat[-1, 1], "+")
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=n.par, nobs=nrow(DH)))
  class(out) <- c("occupancy", "list")
  return(out)
}
