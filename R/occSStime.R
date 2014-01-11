
# New verion 2013-02-27 without the time argument (no time=FALSE option)

occSStime <-
function(DH, model=p~1, data=NULL, ci=0.95, plot=TRUE)  {
  # DH is a 1/0 matrix of detection histories, sites x occasions
  # model is a 2-sided formula for probability of detection, eg, model = p ~ habitat.
  # data is a DATA FRAME with a row for each capture occasion and columns for time covariates.
  # ci is the required confidence interval.

  # Sanity checks and such:
  DH <- as.matrix(DH)  # in case it's a data frame
	nocc <- ncol(DH)
  if (nocc < 2)
    stop("More than one survey occasion is needed")
  stopifnot(is.null(data) || nrow(data) == nocc)
  
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  # Standardise the model:
  model <- stdModel(model, defaultModel=list(p=~1))

  # Add built-in covars to the data frame
  data$.time <- as.factor(1:nocc)
  data$.Time <- 1:nocc
  pDf <- as.data.frame(data)
  # Standardise!!
  
  # Do the model matrix for p:
  pModMat <- model.matrix(model$p, pDf)
  pK <- ncol(pModMat)
  K <- pK + 1
  
  beta.mat <- matrix(NA_real_, K, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("psi",
    paste("p:", colnames(pModMat)))
  lp.mat <- matrix(NA_real_, nocc + 1, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c("psi", paste0("p", 1:nocc))
  logLik <- NA_real_
  varcov <- NULL    # ????
  
  if(ncol(DH) > 1 && sum(DH, na.rm=TRUE) > 0)  {
    # Negative log-likelihood function:
    nll <- function(params) {
      psi <- plogis(params[1])
      pBeta <- params[-1]
      p <- plogis(pModMat %*% pBeta)
      p.DH <- sweep(DH, 2, p, "*") + sweep((1-DH), 2, (1-p), "*")
      llh <- sum(log(psi * apply(p.DH, 1, prod, na.rm=TRUE) + 
          (1 - psi) * (rowSums(DH, na.rm=TRUE) == 0)))
      return(min(-llh, .Machine$double.xmax)) # min(..) stops Inf being returned
    }
    params <- rep(0, K)
    res <- nlm(nll, params, hessian=TRUE)
    if(res$code < 3)  {  # exit code 1 or 2 is ok.
      beta.mat[,1] <- res$estimate
      lp.mat[, 1] <- c(beta.mat[1], pModMat %*% beta.mat[-1,1])
      varcov0 <- try(solve(res$hessian), silent=TRUE)
      if (!inherits(varcov0, "try-error") && all(diag(varcov0) > 0)) {
        varcov <- varcov0
        SE <- sqrt(diag(varcov))
        beta.mat[, 2] <- SE
        beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
        SElp <- c(sqrt(varcov[1,1]),
                sqrt(diag(pModMat %*% varcov[-1,-1] %*% t(pModMat))))
        lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
        logLik <- -res$minimum
      }
    }
    # Do the plot
    if(plot) {
      real.p <- plogis(lp.mat[-1, ])
      ylim <- range(0, real.p, na.rm=TRUE)
      plot(1:nocc, real.p[, 1], type='l', ylim=ylim,
        xlab="Time", ylab="Probability of detection")
      lines(1:nocc, real.p[, 2], lty=3)
      lines(1:nocc, real.p[, 3], lty=3)
    }
  }
  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(lp.mat),
              logLik=c(logLik=logLik, df=K, nobs=nrow(DH)))
  class(out) <- c("wiqid", "list")
  return(out)
}
