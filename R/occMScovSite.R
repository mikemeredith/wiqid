
# Multiseason occupancy
# This version allows for site covars and season differences,
#   but not (yet) seasonal covariates
#   eg psi(hab) gamma(hab+season) ...

# See MacKenzie et al (2006) "Occupancy..." p194ff

# function Prh1A is defined in the file occMSseason.R

occMScovSite <- function(DH, occsPerSeason,
             model=list(psi1~1, gamma~1, epsilon~1, p~1),
             data=NULL, ci=0.95) {    
  # ** DH is detection data in a 1/0/NA matrix or data frame, sites in rows, 
  #    detection occasions in columns..
  # ** occsPerSeason is a scalar or vector with the number of occasions per season
  # ci is the required confidence interval.
             
  if(ci > 1 | ci < 0.5)
    stop("ci must be between 0.5 and 1")
  alf <- (1 - ci[1]) / 2
  crit <- qnorm(c(alf, 1 - alf))

  # Check for all-NA rows (eg, Grand Skinks data set!)
  #  remove them (alternative would be to flag up with, say, last=0 below.
  # allNA <- rowSums(!is.na(DH)) == 0
  # if(any(allNA))  {
    # DH <- DH[!allNA, ]
    # if (!is.null(data))
      # data <- data[!allNA, ]
  # }

  # Deal with occsPerSeason
  nOcc <- ncol(DH)
  if(length(occsPerSeason) == 1)
    occsPerSeason <- rep(occsPerSeason, nOcc/occsPerSeason)
  if(sum(occsPerSeason) != nOcc)
    stop("Detection data do not match occasions per season.")
  nseas <- length(occsPerSeason)
  seasonID <- rep(1:nseas, occsPerSeason)
  # find last season with data
  getLast <- function(dh, grp) {
    if(all(dh==0)) {  # Check for all-NA rows (eg, Grand Skinks data set!)
      return(NA)
    } else {
      return(max(which(rowsum(dh, grp) > 0)))
    }
  }
  last <- as.vector(apply((!is.na(DH))*1, 1, getLast, grp=factor(seasonID)))
  DHplus <- as.matrix(cbind(last, DH))

  # Standardise the model:
  model <- stdModel(model, defaultModel=list(psi1=~1, gamma=~1, epsilon=~1, p=~1))

  # Check data file
  nSites <- nrow(DH)
  siteNames <- rownames(DH)
  if (is.null(siteNames))
    siteNames <- rownames(data)
  if (is.null(siteNames))
    siteNames <- 1:nSites 

  if(!is.null(data)) {
    if(nrow(data) != nSites)
      stop("data must have a row for each site")
    rownames(data) <- NULL
  } else {
    data <- data.frame(.dummy = rep(NA, nSites))
  }
  # cat("Preparing design matrices...") ; flush.console()
  GEseason <- rep(1:(nseas-1), each=nSites)
  ddfGE <- as.data.frame(cbind(data, season=as.factor(GEseason)))
  Pseason <- rep(1:nseas, each=nSites)
  ddfP <- as.data.frame(cbind(data, season=as.factor(Pseason)))

  # Build model matrices  
  psi1Mat <- model.matrix(model$psi1, as.data.frame(data))
  psi1K <- ncol(psi1Mat)
  gamMat <- model.matrix(model$gamma, ddfGE)
  gamK <- ncol(gamMat)
  epsMat <- model.matrix(model$epsilon, ddfGE)
  epsK <- ncol(epsMat)
  pMat <- model.matrix(model$p, ddfP)
  pK <- ncol(pMat)
  K <- psi1K + gamK + epsK + pK
  parID <- rep(1:4, c(psi1K, gamK, epsK, pK))
  
  
  beta.mat <- matrix(NA_real_, K, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c(
    paste("psi1:", colnames(psi1Mat)),
    paste("gam:", colnames(gamMat)),
    paste("eps:", colnames(epsMat)),
    paste("p:", colnames(pMat)))
  lp.mat <- matrix(NA_real_, nSites*(3*nseas-1), 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- c(
    paste0("psi:", siteNames),
    paste0("gamma:", siteNames, ",", GEseason),
    paste0("epsilon:", siteNames, ",", GEseason),
    paste0("p:", siteNames, ",", Pseason))
  logLik <- NA_real_
  varcov <- NULL
  
  nll <- function(param){
    psi1Beta <- param[parID==1]
    gamBeta <- param[parID==2]
    epsBeta <- param[parID==3]
    pBeta <- param[parID==4]
    psi1Prob <- plogis(psi1Mat %*% psi1Beta)
    gamProb <- matrix(plogis(gamMat %*% gamBeta), nrow=nSites)
    epsProb <- matrix(plogis(epsMat %*% epsBeta), nrow=nSites)
    pProb <- matrix(plogis(pMat %*% pBeta), nrow=nSites)
    Prh <- numeric(nSites)
    for(i in 1:nSites) {
      PHI0 <- c(psi1Prob[i], 1-psi1Prob[i])
      PHIt <- array(0, c(2, 2, nseas-1))
      PHIt[1, 1, ] <- 1 - epsProb[i, ]
      PHIt[1, 2, ] <- epsProb[i, ]
      PHIt[2, 1, ] <- gamProb[i, ]
      PHIt[2, 2, ] <- 1 - gamProb[i, ]
      p <- pProb[i, seasonID]
      Prh[i] <- Prh1A(DHplus[i, ], p=p, PHI0=PHI0, PHIt=PHIt, seasonID)
    }
    return(min(-sum(log(Prh)), .Machine$double.xmax))
  }

  # cat("done\n")
  # cat("Maximizing likelihood...") ; flush.console()
  start <- rep(0, K)
  res <- nlm(nll, start, hessian=TRUE)
  # cat("done\n")
  # cat("Organizing output...") ; flush.console()
  if(res$code > 2)   # exit code 1 or 2 is ok.
    warning(paste("Convergence may not have been reached (code", res$code, ")"))
  beta.mat[,1] <- res$estimate
  lp.mat[, 1] <- c(psi1Mat %*% beta.mat[parID==1, 1],
                   gamMat %*% beta.mat[parID==2, 1],
                   epsMat %*% beta.mat[parID==3, 1],
                   pMat %*% beta.mat[parID==4, 1])
  varcov0 <- try(solve(res$hessian), silent=TRUE)
  if (!inherits(varcov0, "try-error") && all(diag(varcov0) > 0)) {
    varcov <- varcov0
    SE <- sqrt(diag(varcov))
    beta.mat[, 2] <- SE  # tidy later
    beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, res$estimate, "+")
    temp <- c(
       diag(psi1Mat %*% varcov[parID==1, parID==1] %*% t(psi1Mat)),
       diag(gamMat %*% varcov[parID==2, parID==2] %*% t(gamMat)),
       diag(epsMat %*% varcov[parID==3, parID==3] %*% t(epsMat)),
       diag(pMat %*% varcov[parID==4, parID==4] %*% t(pMat)))
    if(all(temp >= 0))  {
      SElp <- sqrt(temp)
      lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")
      logLik <- -res$minimum
    }
  }
  # cat("done\n") ; flush.console()

  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=nrow(DH)))
  class(out) <- c("wiqid", "list")
  return(out)
}
