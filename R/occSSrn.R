occSSrn <-
function(y, n)  {
  # n is a vector with the number of occasions at each site.
  # y is a vector with the number of detections at each site.
  if(length(n) == 1)
    n <- rep(n, length(y))
  if(length(y) != length(n))
    stop("y and n must have the same length")
	# Starting values:
  beta.mat <- matrix(NA_real_, 2, 3)
  AIC <- NA_real_
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
      AIC <- 2*res$minimum + 4
      if (det(res$hessian) > 1e-6) {
        SE <- sqrt(diag(solve(res$hessian)))
        crit <- qnorm(c(0.025, 0.975))
        beta.mat[, 2:3] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      }
    }
  }
	lambda <- exp(beta.mat[1,])
	out.mat <- rbind(1-dpois(0, lambda), lambda, plogis(beta.mat[2,]))
  colnames(out.mat) <- c("est", "lowCI", "uppCI")
  rownames(out.mat) <- c("psiHat", "lambdaHat", "rHat")
  attr(out.mat, "AIC") <- AIC
  return(out.mat)
}
