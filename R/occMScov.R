
# Multiseason occupancy
# This version allows for site covars and season differences,
#   but not (yet) seasonal covariates
#   eg psi(hab) gamma(hab+season) ...

# See MacKenzie et al (2006) "Occupancy..." p194ff

occMScov <- function(DH, occsPerSeason,
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

  # Deal with occsPerSeason
  nOcc <- ncol(DH)
  if(length(occsPerSeason) == 1)
    occsPerSeason <- rep(occsPerSeason, nOcc/occsPerSeason)
  if(sum(occsPerSeason) != nOcc)
    stop("Detection data do not match occasions per season.")
  nseas <- length(occsPerSeason)
  seasonID <- rep(1:nseas, occsPerSeason)

  # Standardise the model:
  if(inherits(model, "formula"))
    model <- list(model)
  model <- stdform (model)
  model0 <- list(psi1=~1, gamma=~1, epsilon=~1, p=~1)
  model <- replace (model0, names(model), model)

  # Check data file
  nSites <- nrow(DH)
  siteNames <- 1:nSites  # Why not use names in DH?
  if(!is.null(data)) {
    if(nrow(data) != nSites)
      stop("data must have a row for each sites")
    rownames(data) <- NULL
  } else {
    data <- data.frame(.dummy = rep(NA, nSites))
  }
  cat("Preparing design matrices...") ; flush.console()
  # ddfPSI <- data
  # ddfPSI$.dummy <- rep(NA, nSites)
  # ddfGE <- data
  GEseason <- rep(1:(nseas-1), each=nSites)
  # ddfGE$season <- as.factor(GEseason)
  ddfGE <- as.data.frame(cbind(data, season=as.factor(GEseason)))
  # ddfP <- data
  Pseason <- rep(1:nseas, each=nSites)
  # ddfP$season <- as.factor(Pseason)
  ddfP <- as.data.frame(cbind(data, season=as.factor(Pseason)))

  # Build model matrices  
  # psi1Mat <- model.matrix(model$psi1, as.data.frame(ddfPSI))
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
    paste("ps1:", colnames(psi1Mat)),
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

  nll <- function(param){
    cat("+") ; flush.console()
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
      Prh[i] <- Prh1A(DH[i, ], p=p, PHI0=PHI0, PHIt=PHIt, seasonID)
    }
    return(min(-sum(log(Prh)), .Machine$double.xmax))
  }

  cat("done\n")
  cat("Maximizing likelihood...\n|") ; flush.console()
  start <- rep(0, K)
  res <- nlm(nll, start, hessian=TRUE)
  cat("| done\n")
  cat("Organizing output...") ; flush.console()
  if(res$code < 3)  {  # exit code 1 or 2 is ok.
    beta.mat[,1] <- res$estimate
    lp.mat[, 1] <- c(psi1Mat %*% beta.mat[parID==1, 1],
                     gamMat %*% beta.mat[parID==2, 1],
                     epsMat %*% beta.mat[parID==3, 1],
                     pMat %*% beta.mat[parID==4, 1])
    varcov <- try(solve(res$hessian), silent=TRUE)
    if (!inherits(varcov, "try-error") && all(diag(varcov) > 0)) {
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
  }
  cat("done\n") ; flush.console()

  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=nrow(DH)))
  class(out) <- c("wiqid", "list")
  return(out)
}

# ....................................................

# A function to get Pr(dh) for a single detection history,
#   ie, one row of DH. This version has a 3-D array for PHIt

# Not exported

# dh is a 0/1/NA of length equal total no. of surveys
# p is a scalar, or vector of detection probs of length equal to 
#   dh
# PHI0 is the vector c(psi1, 1-psi1)
# PHIt is a 2 x 2 x (nseas-1) array, where
#   PHIt[,,t] = matrix(c(1-eps[t], gam[t], eps[t], 1-gam[t]), 2)
# seasonID is a vector of length equal to dh, identifying the season.

Prh1A <- function(dh, p, PHI0, PHIt, seasonID) {
  if(all(is.na(dh)))
    return(1)
  # last season with data:
  last <- max(which(tapply(!is.na(dh), seasonID, sum) > 0)) #####
  stopifnot(last > 1) # TODO deal with this properly!
  pvec <- p * dh + (1-p)*(1-dh)
  res <- PHI0
  for(t in 1:(last-1)) {
    # if(t == 0) break    # happens if last = 1
    if(!all(is.na(pvec[seasonID==t])))  {
      D <- diag(c(prod(pvec[seasonID==t], na.rm=TRUE),
                  1-max(dh[seasonID==t], na.rm=TRUE)))
      res <- res %*% D
    }
    res <- res %*% PHIt[, , t]
  }
  PT <- c(prod(pvec[seasonID==last], na.rm=TRUE), 1-max(dh[seasonID==last], na.rm=TRUE))
  res <- res %*% PT
  return(res)
}
  