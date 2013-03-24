occSSrn <-
function(y, n)  {
  # y is a vector with the number of detections at each site.
  # n is a vector with the number of occasions at each site.
  if(length(n) == 1)
    n <- rep(n, length(y))
  if(length(y) != length(n))
    stop("y and n must have the same length")
  if(any(y > n))
    stop("y cannot be greater than n")
	# Starting values:
  beta.mat <- matrix(NA_real_, 2, 4) 
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("lambda", "r")
  logLik <- NA_real_
  if(sum(n) > 0 && sum(y) > 0) {    # If all n's are 0, no data available.
    params <- c(0, 0)
    Nmax <- 100 # See later if this is sensible
    # Negative log-likelihood function:
    nll <- function(params) {
      lambda <- exp(params[1])
      r <- plogis(params[2])
      rpart <- (1-r)^(0:Nmax)
      Npart <- dpois(0:Nmax, lambda)
      llh <- 0
      for(i in seq_along(n)) {
        llh <- llh + log(sum((1-rpart)^y[i] * rpart^(n[i]-y[i]) * Npart))
      }
      return(min(-llh, .Machine$double.xmax)) # min(..) stops Inf being returned
    }
    res <- nlm(nll, params, hessian=TRUE)
    if(res$code < 3)  {  # exit code 1 or 2 is ok.
      beta.mat[,1] <- res$estimate
      varcov <- try(solve(res$hessian), silent=TRUE)
      if (!inherits(varcov, "try-error") && 
          all(diag(varcov) > 0)) {
        SE <- sqrt(diag(varcov))
        beta.mat[, 2] <- SE
        crit <- qnorm(c(0.025, 0.975))
        beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
        logLik <- -res$minimum
      }
    }
  }
	lambda <- exp(beta.mat[1, -2])
	real <- rbind(1-dpois(0, lambda), lambda, plogis(beta.mat[2, -2]))
  colnames(real) <- c("est", "lowCI", "uppCI")
  rownames(real) <- c("psiHat", "lambdaHat", "rHat")
  out <- list(call = match.call(),
              beta = beta.mat,
              real = real,
              logLik = c(logLik=logLik, df=2, nobs=length(y)))
  class(out) <- c("occupancy", "list")
  return(out)

}
