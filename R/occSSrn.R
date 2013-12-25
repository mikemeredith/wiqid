occSSrn <-
function(y, n, ci=0.95)  {
  # y is a vector with the number of detections at each site.
  # n is a vector with the number of occasions at each site.
  # ci is the required confidence interval.
  y <- round(y)
  n <- round(n)
  if(length(n) == 1)
    n <- rep(n, length(y))
  if(length(y) != length(n))
    stop("y and n must have the same length")
  if(any(y > n))
    stop("y cannot be greater than n")
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

	# Starting values:
  beta.mat <- matrix(NA_real_, 2, 4) 
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("lambda", "r")
  logLik <- NA_real_
  varcov <- NULL
  
  if(sum(n) > 0 && sum(y) > 0 && any(y < n)) {    # If all n's are 0, no data available.
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
      varcov0 <- try(solve(res$hessian), silent=TRUE)
      if (!inherits(varcov0, "try-error") && all(diag(varcov0) > 0)) {
        varcov <- varcov0
        SE <- sqrt(diag(varcov))
        beta.mat[, 2] <- SE
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
              beta.vcv = varcov,
              real = real,
              logLik = c(logLik=logLik, df=2, nobs=length(y)))
  class(out) <- c("wiqid", "list")
  return(out)

}
