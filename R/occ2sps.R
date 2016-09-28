
# Single-season 2-species occupancy function, based on Richmond et al 2010.

# Richmond et al (2010) interpret "species A detected" to mean
#  "species A detected at the site on the occasion in question."
# See the example in Eqn2 p2038.

# DHA is detection history of the dominant species
# DHB is detection history of the subordinate species
# model is now a list of two-sided formulae (new 2015-09-19)
#  The default is list(psiA~1, psiBa~1, pA~1, pB~1). If not included in the model
#  list, the remaining parameters are assigned the following values:
#  psiBA <- psiBa, rA <- pA, rBa <- pB, rBA <- rBa.

occ2sps <- function(DHA, DHB, model=NULL, data=NULL, ci=0.95, verify=TRUE)  {

  DHA <- as.matrix(DHA)
  DHB <- as.matrix(DHB)
  if(verify) {
    stopifnot(all.equal(dim(DHA), dim(DHB)))
    DHA <- verifyDH(DHA, allowNA = TRUE)
    DHB <- verifyDH(DHB, allowNA = TRUE)
    # Check that the NAs match up
    stopifnot(all.equal(is.na(DHA), is.na(DHB), check.attributes=FALSE))
  }
  
  # Standardise the model:
  model <- stdModel(model, list(psiA=~1, psiBa=~1, pA=~1, pB=~1))
  # Check for invalid submodels in 'model':
  parNames <- c("psiA", "psiBa", "psiBA", "pA", "pB", "rA", "rBa", "rBA")
  ok <- names(model) %in% parNames
  if(any(!ok))
    stop("Invalid submodels for: ", paste(names(model)[!ok], collapse=", "))
  # modPars is a vector of length 8 which maps the submodels needed to the
  #   elements of 'model':
  modPars <- pmatch(parNames, names(model))
  names(modPars) <- parNames
  if(is.null(model$psiBA))
    modPars[3] <- modPars[2]  # psiBA <- psiBa
  if(is.null(model$rA))
    modPars[6] <- modPars[4]  # rA <- pA
  if(is.null(model$rBa))
    modPars[7] <- modPars[5]  # rBa <- pB
  if(is.null(model$rBA))
    modPars[8] <- modPars[7]  # rBA <- rBa

  if(is.null(data))
    return(occ2sps0(DHA, DHB, modPars, ci=ci))

  M <- length(model)  # Number of elements in the model
  nSites <- nrow(DHA)
  site.names <- rownames(data)
  if(is.null(site.names))
    site.names <- 1:nSites
  crit <- fixCI(ci)

  data <- as.data.frame(stddata(data, nocc=NULL))

  # Build model matrices
  modMatList <- vector('list', M)
  for(i in 1:M)
    modMatList[[i]] <- model.matrix(model[[i]], data)
  parK <- sapply(modMatList, ncol)    # Number of parameters for each model matrix
  K <- sum(parK)  # total number of parameters
  idK <- rep(1:M, parK)  # specifies which of the K parameters belongs to each model matrix
  # Get coefficient names
  coefNames <- paste(rep(names(model), parK),
      unlist(lapply(modMatList, colnames)), sep=":")

  # function to get the occupancy matrix
  # psiX is a matrix with columns psiA, psiBa, psiBA
  getlogPHI <- function(psiX) {
    log(cbind(psiX[, 1] * psiX[, 3],             # both
      psiX[, 1] * (1 - psiX[, 3]),       # A only
      (1 - psiX[, 1]) * psiX[, 2],       # B only
      (1 - psiX[, 1]) * (1 - psiX[, 2]))) # neither
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
  beta.mat <- matrix(NA_real_, K, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- coefNames
  logLik <- NA_real_
  varcov <- NULL
  lp.mat <- matrix(NA_real_, nSites * 8, 3)
  colnames(lp.mat) <- c("est", "lowCI", "uppCI")
  rownames(lp.mat) <- as.vector(t(outer(parNames, site.names, paste, sep=":")))

  # Do the neg log lik function:
  real0 <- matrix(NA, nSites, M)
  nll <- function(params) {
    for(i in 1:M) {
      betas <- params[idK == i]
      real0[, i] <- plogis(modMatList[[i]] %*% betas)
    }
    real <- real0[, modPars]
    logPHI <- getlogPHI(real[, 1:3])
    loglik <- numeric(nSites)
    for(i in 1:nSites)
      loglik[i] <- logSumExp(logPHI[i, ] + getlogP(dhA=DHA[i, ], dhB=DHB[i, ], pX=real[i, 4:8]) )
    return(min(-sum(loglik), .Machine$double.xmax))
  }

  # Run mle estimation with optim:
  params <- rep(0, K)
  res <- optim(params, nll, method="L-BFGS-B", lower=-10, upper=10, hessian=TRUE)
  if(res$convergence > 0) {
    warning(paste("Convergence may not have been reached.", res$message))
  } else {
    logLik <- -res$value
  }

  beta.mat[,1] <- res$par
  lp.mat0 <- matrix(NA, nSites, M)
  for(i in 1:M) {
    betas <- res$par[idK == i]
    lp.mat0[, i] <- modMatList[[i]] %*% betas
  }
  lp.mat[, 1] <- as.vector(lp.mat0[, modPars])

  varcov0 <- try(chol2inv(chol(res$hessian)), silent=TRUE)
  # if (!inherits(varcov0, "try-error") && all(diag(varcov0) > 0)) {
  if (!inherits(varcov0, "try-error")) {
    varcov <- varcov0
    SE <- suppressWarnings(sqrt(diag(varcov)))
    beta.mat[, 2] <- SE
    beta.mat[, 3:4] <- sweep(outer(SE, crit), 1, beta.mat[, 1], "+")
  }
  SElp0 <- matrix(NA, nSites, M)
  for(i in 1:M) {
    varcov1 <- varcov[idK == i, idK == i]
    SElp0[, i] <- sqrt(diag(modMatList[[i]] %*% varcov1 %*% t(modMatList[[i]])))
  }
  SElp <- as.vector(SElp0[, modPars])
  lp.mat[, 2:3] <- sweep(outer(SElp, crit), 1, lp.mat[, 1], "+")

  out <- list(call = match.call(),
              beta = beta.mat,
              beta.vcv = varcov,
              real = plogis(lp.mat),
              logLik = c(logLik=logLik, df=K, nobs=nSites))
  class(out) <- c("wiqid", "list")
  return(out)
}
