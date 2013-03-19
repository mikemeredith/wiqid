
# New verion 2013-02-27 without the time argument (no time=FALSE option)

occSS.t <-
function(DH)  {
  # DH is a 1/0 matrix of detection histories, sites x occasions
  # Check for columns with all NAs:
  DH <- as.matrix(DH)  # in case it's a data frame
	nocc <- ncol(DH)
  NAcol <- colSums(!is.na(DH)) == 0
  if(any(NAcol))
    DH <- DH[, !NAcol]
	n.par <- ncol(DH)+1
  beta.mat <- matrix(NA_real_, n.par, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames.beta.mat <- "psi"
  if(any(!NAcol))
    rownames.beta.mat <- c("psi", paste("p", (1:nocc)[!NAcol], sep=""))
  rownames(beta.mat) <- rownames.beta.mat
  # AIC <- NA_real_
  logLik <- NA_real_
  if(ncol(DH) > 1 && sum(DH, na.rm=TRUE) > 0)  {
    # Negative log-likelihood function:
    nll <- function(params) {
      psi <- plogis(params[1])
      p <- plogis(params[-1])
      p.DH <- sweep(DH, 2, p, "*") + sweep((1-DH), 2, (1-p), "*")
      llh <- sum(log(psi * apply(p.DH, 1, prod, na.rm=TRUE) + 
          (1 - psi) * (rowSums(DH, na.rm=TRUE) == 0)))
      # return(min(-llh, .Machine$double.xmax)) # min(..) stops Inf being returned
      return(-llh) # try it and see
    }
    params <- rep(0, n.par)
    res <- nlm(nll, params, hessian=TRUE)
    if(res$code < 3)  {  # exit code 1 or 2 is ok.
      beta.mat[,1] <- res$estimate
      # AIC <- 2*res$minimum + 2 * n.par
      logLik <- -res$minimum
      if (det(res$hessian) > 1e-6) {
        SE <- sqrt(diag(solve(res$hessian)))
        beta.mat[, 2] <- SE
        crit <- qnorm(c(0.025, 0.975))
        beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      }
    }
  }
  real.mat <- matrix(NA, nocc+1, 3)
  real.mat[c(TRUE, !NAcol), ] <- plogis(beta.mat[, -2])
  rownames(real.mat) <- c("psiHat",paste("p",1:(nocc),"Hat", sep=""))
  colnames(real.mat) <- c("est", "lowCI", "uppCI")
  out <- list(call = match.call(),
              beta = beta.mat,
              real = real.mat,
              logLik=c(logLik=logLik, df=n.par, nobs=nrow(DH)))
  class(out) <- c("occupancy", "list")
  return(out)
}
