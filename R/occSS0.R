occSS0 <-
function(y, n, ci=0.95) {
  # n is a vector with the number of occasions at each site.
  # y is a vector with the number of detections at each site.
  # ci is the required confidence interval.
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

  beta.mat <- matrix(NA_real_, 2, 3)
  colnames(beta.mat) <- c("est", "lowCI", "uppCI")
  rownames(beta.mat) <- c("psiHat", "pHat")
  logLik <- NA_real_
  varcov <- NULL
  
  if(sum(n) > 0) {    # If all n's are 0, no data available.
    nll <- function(params) {
       psi <- plogis(params[1])
       p <- plogis(params[2])
       prob <- psi * p^y * (1-p)^(n - y) + (1 - psi) * (y==0)
      return(min(-sum(log(prob)), .Machine$double.xmax))
    }
    params <- rep(0,2)
    res <- nlm(nll, params, hessian=TRUE)
    if(res$code < 3)  {  # exit code 1 or 2 is ok.
      beta.mat[,1] <- res$estimate
      varcov <- try(solve(res$hessian), silent=TRUE)
      if (!inherits(varcov, "try-error") && all(diag(varcov) > 0)) {
        SE <- sqrt(diag(varcov))
        beta.mat[, 2:3] <- sweep(outer(SE, crit), 1, res$estimate, "+")
        logLik <- -res$minimum
      }
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(beta.mat),
              logLik = c(logLik=logLik, df=2, nobs=length(y)))
  class(out) <- c("wiqid", "list")
  return(out)
}
