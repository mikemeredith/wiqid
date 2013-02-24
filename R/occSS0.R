occSS0 <-
function(y, n) {
  # n is a vector with the number of occasions at each site.
  # y is a vector with the number of detections at each site.
  if(length(n) == 1)
    n <- rep(n, length(y))
  if(length(y) != length(n))
    stop("y and n must have the same length")
  beta.mat <- matrix(NA_real_, 2, 3)
  AIC <- NA_real_
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
      AIC <- 2*res$minimum + 4
      if (det(res$hessian) > 0) {
        SE <- sqrt(diag(solve(res$hessian)))
        crit <- qnorm(c(0.025, 0.975))
        beta.mat[, 2:3] <- sweep(outer(SE, crit), 1, res$estimate, "+")
      }
    }
  }
  out.mat <- plogis(beta.mat)
  colnames(out.mat) <- c("est", "lowCI", "uppCI")
  rownames(out.mat) <- c("psiHat", "pHat")
  attr(out.mat, "AIC") <- AIC
  return(out.mat)
}
