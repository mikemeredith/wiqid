
# Single-season 2-species occupancy function, based on Richmond et al 2010.

# Richmond et al (2010) interpret "species A detected" to mean
#  "species A detected at the site on the occasion in question."
# See the example in Eqn2 p2038.

# DHA is detection history of the dominant species
# DHB is detection history of the subordinate species
# modelSpec is a 3-digit number, where
  # the first digit specifies occupancy model (1: psiBA = psiBa, 2: psiBA != psiBa),
  # the second, detection probability of species A (1: pA = rA, 2: pA != rA),
  # the third, detection probability of species B (1: pB = rBA = rBa,
  #   2: pB != rBA = rBa, 3: pB = rBa != rBA, 4: pB != rBA != rBa)

occ2sps <- function(DHA, DHB, modelSpec=111, ci=0.95)  {
  DHA <- as.matrix(DHA)
  DHB <- as.matrix(DHB)
  # Check that the NAs match up
  stopifnot(all.equal(is.na(DHA), is.na(DHB), check.attributes=FALSE))
  # Check for rows with all NAs ### this probably doesn't matter
  # bad <- rowSums(!is.na(DHA)) == 0
  # if(sum(bad) > 0) {
    # warning()
    # DHA <- DHA[!bad, ]
    # DHB <- DHB[!bad, ]
  # }
  nSites <- nrow(DHA)

  crit <- fixCI(ci)

  # modPars is a vector of length 8 which maps the coefficients estimated to the
  #   model parameters; order is
  # psiA, psiBA, psiBa, pA, pB, rA, rBA, rBa.
  modPars <- 1:8
  modDesc <- c("psiBa != psiBA", "rA != pA", NA)
  if(modelSpec < 200) {
    modPars[3] <- 2
    modDesc[1] <- "psiBa = psiBA"
  }
  if(modelSpec %% 100 < 20)  {
    modPars[6] <- 4
    modDesc[2] <- "rA = pA"
  }
  modPars[c(5, 7, 8)] <- switch(as.character(modelSpec %% 10),
    "1" = c(5, 5, 5),   # 1: pB = rBA = rBa
    "2" = c(5, 7, 7),   # 2: pB != rBA = rBa
    "3" = c(5, 7, 5),   # 3: pB = rBa != rBA
    "4" = c(5, 7, 8),   # 4: pB != rBA != rBa)
    stop("modelSpec ", modelSpec, " not recognized."))
  modPars <- as.integer(as.factor(modPars))
  modDesc[3] <- switch(as.character(modelSpec %% 10),
    "1" = "pB = rBA = rBa",
    "2" = "pB != rBA = rBa",
    "3" = "pB = rBa != rBA",
    "4" = "pB != rBA != rBa")

  # get the occupancy vector
  # psiX is a vector with elements psiA, psiBA, psiBa
  getPHI <- function(psiX) {
    c(psiX[1] * psiX[2],             # both
      psiX[1] * (1 - psiX[2]),       # A only
      (1 - psiX[1]) * psiX[3],       # B only
      (1 - psiX[1]) * (1 - psiX[3])) # neither
  }

  # Do the detection vector for one site
  # pX is a vector with elements pA, pB, rA, rBA, rBa
  getP <- function(dhA, dhB, pX)  {
    # prob of detecting B if both present conditional on detection of A
    probCapB <- dhA * pX[4] + (1 - dhA) * pX[5]
    c(
      # Both sps present, use the r's
      prod(dhA * pX[3] + (1 - dhA) * (1 - pX[3]),      # A
        dhB * probCapB + (1 - dhB) * (1 - probCapB), na.rm=TRUE),  # B
      # Sps A present, B absent, use pA
      if(sum(dhB, na.rm=TRUE) > 0) { 0 } else {
        prod(dhA * pX[1] + (1 - dhA) * (1 - pX[1]), na.rm=TRUE) # A
      },
      # Sps A absent, B present, use pB
      if(sum(dhA, na.rm=TRUE) > 0) { 0 } else {  
        prod(dhB * pX[2] + (1 - dhB) * (1 - pX[2]), na.rm=TRUE) # B
      },
      # Neither present
      if(sum(dhA, dhB, na.rm=TRUE) > 0) { 0 } else { 1 } )
  }

  # objects to hold output
  beta.mat <- matrix(NA_real_, 8, 4)
  colnames(beta.mat) <- c("est", "SE", "lowCI", "uppCI")
  rownames(beta.mat) <- c("psiA", "psiBA", "psiBa",
    "pA", "pB", "rA", "rBA", "rBa")
  logLik <- NA_real_
  varcov <- NULL

  # Do the neg log lik function:
  nll <- function(params) {
    real <- plogis(params[modPars])
    PHI <- getPHI(real[1:3])
    lik <- numeric(nSites)
    for(i in 1:nSites)
      lik[i] <- PHI %*% getP(dhA=DHA[i, ], dhB=DHB[i, ], pX=real[4:8])  
    return(min(-sum(log(lik)), .Machine$double.xmax))
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
              logLik = c(logLik=logLik, df=length(unique(modPars)), nobs=nSites))
  class(out) <- c("wiqid", "list")
  return(out)
}
