
# Single-season 2-species occupancy function, based on Richmond et al 2010.

# Simplified analysis with no covariates.
# Called by occ2sps, not exported
# modPars is a vector of length 8 which maps the coefficients estimated to the
#   model parameters; order is
#   psiA, psiBA, psiBa, pA, pB, rA, rBa, rBA.

occ2sps0 <- function(DHA, DHB, modPars, ci=0.95)  {
  nSites <- nrow(DHA)
  crit <- fixCI(ci)

  # get the occupancy vector
  # psiX is a vector with elements psiA, psiBa, psiBA
  getlogPHI <- function(psiX) {
    log(c(psiX[1] * psiX[3],             # both
      psiX[1] * (1 - psiX[3]),       # A only
      (1 - psiX[1]) * psiX[2],       # B only
      (1 - psiX[1]) * (1 - psiX[2]))) # neither
  }

  # Do the detection vector for one site
  # pX is a vector with elements pA, pB, rA, rBa, rBA
  getlogP <- function(dhA, dhB, pX)  {
    # prob of detecting B if both present conditional on detection of A
    probCapB <- dhA * pX[5] + (1 - dhA) * pX[4]
    c(
      # Both sps present, use the r's
      sum(log(dhA * pX[3] + (1 - dhA) * (1 - pX[3])),      # A
        log(dhB * probCapB + (1 - dhB) * (1 - probCapB)), na.rm=TRUE),  # B
      # Sps A present, B absent, use pA
      if(sum(dhB, na.rm=TRUE) > 0) { -Inf } else {
        sum(log(dhA * pX[1] + (1 - dhA) * (1 - pX[1])), na.rm=TRUE) # A
      },
      # Sps A absent, B present, use pB
      if(sum(dhA, na.rm=TRUE) > 0) { -Inf } else {
        sum(log(dhB * pX[2] + (1 - dhB) * (1 - pX[2])), na.rm=TRUE) # B
      },
      # Neither present
      if(sum(dhA, dhB, na.rm=TRUE) > 0) { -Inf } else { 0 } )
  }

  # objects to hold output
  beta.mat <- matrix(NA_real_, 8, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- names(modPars)
  logLik <- NA_real_
  varcov <- NULL

  # Do the neg log lik function:
  nll <- function(params) {
    real <- plogis(params[modPars])
    logPHI <- getlogPHI(real[1:3])
    loglik <- numeric(nSites)
    for(i in 1:nSites)
      loglik[i] <- logSumExp(logPHI + getlogP(dhA=DHA[i, ], dhB=DHB[i, ], pX=real[4:8]) )
    return(min(-sum(loglik), .Machine$double.xmax))
  }

  # Run mle estimation with optim:
  params <- rep(0, length(unique(modPars)))
  res <- optim(params, nll, method="L-BFGS-B", lower=-10, upper=10, hessian=TRUE)
  if(res$convergence > 0) {
    warning(paste("Convergence may not have been reached.", res$message))
  } else {
    logLik <- -res$value
  }

  beta.mat[,1] <- res$par[modPars]
  # varcov0 <- try(solve(res$hessian), silent=TRUE)
  varcov0 <- try(chol2inv(chol(res$hessian)), silent=TRUE)
  # if (!inherits(varcov0, "try-error") && all(diag(varcov0) > 0)) {
  if (!inherits(varcov0, "try-error")) {
    varcov <- varcov0
    SE <- suppressWarnings(sqrt(diag(varcov))[modPars])
    beta.mat[, 2] <- SE
    beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, beta.mat[, 1], "+")
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(beta.mat[, -2]),
              logLik = c(logLik=logLik, df=length(unique(modPars)), nobs=nSites),
              ci = ci)
  class(out) <- c("wiqid", "list")
  return(out)
}
