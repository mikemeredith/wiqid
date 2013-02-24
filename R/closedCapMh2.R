closedCapMh2 <-
function(freq, n.occ = length(freq)) {
  # freq is a vector of capture frequencies.
  freq <- c(freq, rep(0, n.occ-length(freq)))
  # n.occ is the total number of capture occasions
  N.cap <- sum(freq)  # Number of individual animals captured
  beta.mat <- matrix(NA_real_, 4, 3)
  AIC <- NA_real_

  if(sum(freq[-1]) > 1)  {  # Do checks here
    nll1 <- function(params) {
      f0 <- min(exp(params[1]), 1e+300, .Machine$double.xmax)
      f.vect <- c(f0, freq)
      pi <- plogis(params[2])
      p1 <- plogis(params[3])
      p2 <- p1 * plogis(params[4]) # ensures that p2 <= p1
      bin.co <- lgamma(N.cap + f0 + 1) - lgamma(f0 + 1)
      tmp <- numeric(length(f.vect))
      for(i in seq_along(f.vect))
         tmp[i] <- pi * p1^(i-1) * (1-p1)^(n.occ-i+1) +
                  (1-pi) * p2^(i-1) * (1-p2)^(n.occ-i+1)
      llh <- bin.co + sum(f.vect * log(tmp))
      return(min(-llh, .Machine$double.xmax))
    }

    nll2 <- function(params) {
      f0 <- min(exp(params[1]), 1e+300, .Machine$double.xmax)
      f.vect <- c(f0, freq)
      pi <- plogis(params[2])
      p1 <- plogis(params[3])
      p2 <- plogis(params[4])
      bin.co <- lgamma(N.cap + f0 + 1) - lgamma(f0 + 1)
      tmp <- numeric(length(f.vect))
      for(i in seq_along(f.vect))
         tmp[i] <- pi * p1^(i-1) * (1-p1)^(n.occ-i+1) +
                  (1-pi) * p2^(i-1) * (1-p2)^(n.occ-i+1)
      llh <- bin.co + sum(f.vect * log(tmp))
      return(min(-llh, .Machine$double.xmax))
    }

    res1 <- suppressWarnings(nlm(nll1, c(log(5), 0,0,0)))
    if(res1$code < 3)  {  # exit code 1 or 2 is ok.
      params <- res1$estimate
      p2 <- plogis(params[3]) * plogis(params[4])
      params[4] <- qlogis(p2)
      #res2 <- suppressWarnings(nlm(nll2, params, hessian=TRUE))
      res2 <- nlm(nll2, params, hessian=TRUE)
      if(res2$code < 3)  {  
        beta.mat[,1] <- res2$estimate
        AIC <- 2*res2$minimum + 8
        if (det(res2$hessian) > 1e-6) {
          SE <- sqrt(diag(solve(res2$hessian)))
          crit <- qnorm(c(0.025, 0.975))
          beta.mat[, 2:3] <- sweep(outer(SE, crit), 1, res2$estimate, "+")
        }
      }
    }
  }
  out.mat <- rbind(exp(beta.mat[1, ]) + N.cap, 
                    plogis(beta.mat[-1, ]))
  colnames(out.mat) <- c("est", "lowCI", "uppCI")
  rownames(out.mat) <- c("Nhat", "piHat", "p1hat", "p2hat")
  attr(out.mat, "AIC") <- AIC
  return(out.mat)
}
