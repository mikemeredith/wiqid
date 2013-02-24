closedCapM0 <-
function(freq, n.occ = length(freq)) {
  # freq is a vector of capture frequencies; trailing zeros are not required.
  # n.occ is the total number of capture occasions
  N.cap <- sum(freq)  # Number of individual animals captured
  n.snap <- sum(freq * (1:length(freq))) # Total number of capture events
  beta.mat <- matrix(NA_real_, 2, 3)
  AIC <- NA_real_
  if(sum(freq[-1]) > 1) {  # Need recaptures
    nll <- function(params) {
      N <- min(exp(params[1]) + N.cap, 1e+300, .Machine$double.xmax)
      p <- plogis(params[2])
      tmp <- lgamma(N + 1) - lgamma(N - N.cap + 1) + n.snap*log(p) +
            (N*n.occ - n.snap)*log(1-p)
      return(min(-tmp, .Machine$double.xmax))
    }
    params <- c(log(5), 0)
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
  out.mat <- rbind(exp(beta.mat[1, ]) + N.cap, 
                    plogis(beta.mat[2, ]))
  colnames(out.mat) <- c("est", "lowCI", "uppCI")
  rownames(out.mat) <- c("Nhat", "pHat")
  attr(out.mat, "AIC") <- AIC
  return(out.mat)
}
