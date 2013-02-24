occSS.t <-
function(DH, time=TRUE)  {
  # DH is a 1/0 matrix of detection histories, sites x occasions
  # Check for columns with all NAs:
	nocc <- ncol(DH)
  NAcol <- colSums(!is.na(DH)) == 0
  if(any(NAcol))
    DH <- DH[, !NAcol]
	# Starting values:
	if(time) {
		params <- rep(0, ncol(DH)+1)
	} else {
		params <- rep(0, 2)
	}
	n.par <- length(params)
  beta.mat <- matrix(NA_real_, n.par, 3)
  AIC <- NA_real_
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
    res <- nlm(nll, params, hessian=TRUE)
    if(res$code < 3)  {  # exit code 1 or 2 is ok.
      beta.mat[,1] <- res$estimate
      AIC <- 2*res$minimum + 2 * n.par
      if (det(res$hessian) > 1e-6) {
        SE <- sqrt(diag(solve(res$hessian)))
        crit <- qnorm(c(0.025, 0.975))
        beta.mat[, 2:3] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      }
    }
  }
  if(time) { 
    out.mat <- matrix(NA, nocc+1, 3)
    out.mat[c(TRUE, !NAcol), ] <- plogis(beta.mat)
    rownames(out.mat) <- c("psiHat",paste("p",1:(nocc),"Hat", sep=""))
  } else {
    out.mat  <- plogis(beta.mat)
    rownames(out.mat) <- c("psiHat", "pHat")
  }
  colnames(out.mat) <- c("est", "lowCI", "uppCI")
  attr(out.mat, "AIC") <- AIC
  return(out.mat)
}
