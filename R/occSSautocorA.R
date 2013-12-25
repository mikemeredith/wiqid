# Single season occupancy with autocorrelation among successive surveys
# Model psi(.) p(.) p'(.) in Appendix A of 
# Hines, Nichols, Royle, MacKenzie, Gopalaswamy, Kumar & Karanth (2010)
#   Tigers on trails: occupancy modeling for cluster sampling. 
#   Ecological Applications, 20, 1456-1466.


occSSautocorA <- function(DH, ci=0.95)  {
  # DH is a 1/0 matrix of detection histories, sites x occasions,
  #   with NO MISSING VALUES
  # ci is the required confidence interval.
  # This returns estimate & CI of psi as usual, plus p and p'

  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  DH <- as.matrix(DH)  # in case it's a data frame
  # Extract summary statistics:
	K <- ncol(DH)
  s <- nrow(DH)
  det <- rowSums(DH) > 0
  sD <- sum(det)
  prev <- cbind(0, DH[, -K])
  n11 <- sum(DH + prev == 2)
  n10 <- sum(prev == 1 & DH == 0)
  n01 <- sum(prev == 0 & DH == 1)
  n00 <- sum((prev + DH)[det, ] == 0)
  stopifnot(n11+n10+n01+n00 == sD*K)

	n.par <- 3
  beta.mat <- matrix(NA_real_, n.par, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("psi", "p", "pPrime")
#  lp.mat <- matrix(NA_real_, K+1, 3)
#  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
#  rownames(lp.mat) <- c("psi", paste0("p", 1:K))
  # AIC <- NA_real_
  logLik <- NA_real_
  varcov <- NULL
  if(ncol(DH) > 1 && sum(DH, na.rm=TRUE) > 0)  {
    # Negative log-likelihood function:
    nll <- function(params) {
      psi <- plogis(params[1])
      p <- plogis(params[2])
      pPrime <- plogis(params[3])
      llh <- sD*log(psi) + n11*log(pPrime) + n10*log(1-pPrime) + n01*log(p) +
          n00*log(1-p) + (s-sD)*log(1 - psi + psi*(1-p)^K)
      return(min(-llh, .Machine$double.xmax)) # min(..) needed to stop Inf being returned
    }
    params <- rep(0, n.par)
    res <- nlm(nll, params, hessian=TRUE)
    if(res$code < 3)  {  # exit code 1 or 2 is ok.
      beta.mat[,1] <- res$estimate
      # AIC <- 2*res$minimum + 2 * n.par
      varcov <- try(solve(res$hessian), silent=TRUE)
      if (!inherits(varcov, "try-error") && all(diag(varcov) > 0)) {
        SE <- sqrt(diag(varcov))
        beta.mat[, 'SE'] <- sqrt(diag(varcov))
        beta.mat[, 3:4] <- sweep(outer(beta.mat[, 'SE'], crit), 1, res$estimate, "+")
        logLik <- -res$minimum
      }
    }
    if(mean(DH, na.rm=TRUE) == 1) {
      beta.mat[3:4 ] <- NA
      # AIC <- NA_real_
      logLik <- NA_real_
    }   ##############################################################
    # lp.mat[1, ]  <- beta.mat[1, -2]
    # lp.mat[-1, 1]  <- beta.mat[2, 1] + beta.mat[3, 1] * Time
    # if(!is.null(varcov)) {
      # pModMat <- cbind(1, Time)
      # SElp <- sqrt(diag(pModMat %*% varcov[-1, -1] %*% t(pModMat)))
      # lp.mat[-1, -1] <- sweep(outer(SElp, crit), 1, lp.mat[-1, 1], "+")
    # }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(beta.mat)[, -2],
              logLik = c(logLik=logLik, df=n.par, nobs=nrow(DH)))
  class(out) <- c("wiqid", "list")
  return(out)
}
