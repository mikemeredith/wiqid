
# Royle-Nichols occupancy models, with abundance-induced heterogeneity in detection probability.

# 'link' argument added 2015-02-20

## WORK IN PROGRESS ## can't yet deal with suvey covariates.

occSSrn <- function(DH, model=NULL, data=NULL,
    ci=0.95, link=c("logit", "probit")) {
  # single-season Royle-Nichols model with site and survey covariates
  # ** DH is detection data in a 1/0/NA matrix or data frame, sites in rows,
  #    detection occasions in columns.
  # ** model is a list of 2-sided formulae for lambda and r; can also be a single
  #   2-sided formula, eg, model = lambda ~ habitat.
  #  NOTE: survey covariates not yet implemented!
  # ** data is a DATA FRAME with single columns for site covariates and a column for each survey occasion for each survey covariate.
  # ci is the required confidence interval.

  if (!is.matrix(DH) && !is.data.frame(DH))
    stop("DH should be a detection history matrix (or data frame)")

  if(TRUE) {  # TODO check that survey covars aren't included in the model
    y <- rowSums(DH, na.rm=TRUE)
    n <- rowSums(!is.na(DH))
    if(is.null(model)) {
      return(occSSrn0(y, n, ci=ci, link=link))
    } else {
      return(occSSrnSite(y, n, model=model, data=data, ci=ci, link=link))
    }
  }
}

# ------------------------------------------------------------------

occSSrn0 <-
function(y, n, ci=0.95, link=c("logit", "probit"))  {
  # Fast version without covariates.
  # y is a vector with the number of detections at each site.
  # n is a vector with the number of occasions at each site.
  # ci is the required confidence interval.
  crit <- fixCI(ci)

  if(match.arg(link) == "logit") {
    plink <- plogis
  } else {
    plink <- pnorm
  }

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
      r <- plink(params[2])
      rpart <- (1-r)^(0:Nmax)
      Npart <- dpois(0:Nmax, lambda)
      llh <- 0
      for(i in seq_along(n)) {
        llh <- llh + log(sum((1-rpart)^y[i] * rpart^(n[i]-y[i]) * Npart))
      }
      return(min(-llh, .Machine$double.xmax)) # min(..) stops Inf being returned
    }
    res <- nlm(nll, params, hessian=TRUE)
    if(res$code > 2)   # exit code 1 or 2 is ok.
      warning(paste("Convergence may not have been reached (code", res$code, ")"))
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
	lambda <- exp(beta.mat[1, -2])
	real <- rbind(1-dpois(0, lambda), lambda, plink(beta.mat[2, -2]))
  colnames(real) <- c("est", "lowCI", "uppCI")
  rownames(real) <- c("psiHat", "lambdaHat", "rHat")
  out <- list(call = match.call(),
              link = match.arg(link),
              beta = beta.mat,
              beta.vcv = varcov,
              real = real,
              logLik = c(logLik=logLik, df=2, nobs=length(y)))
  class(out) <- c("wiqid", "list")
  return(out)
}
