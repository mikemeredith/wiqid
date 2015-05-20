
# Bayesian modelling of normal distribution with JAGS
# ===================================================
# Essentially the same as one-sample version of BEST::BESTmcmc
#  but with normal instead of t-distribution.

Bnormal <- function(y, priors=NULL, doPriorsOnly=FALSE,
    numSavedSteps=1e5, thinSteps=1, burnInSteps = 1000,
    verbose=TRUE, rnd.seed=NULL, parallel=NULL) {

  startTime <- Sys.time()

  if(doPriorsOnly && verbose)
    cat("Warning: The output shows the prior distributions,
      NOT the posterior distributions for your data.\n")
  # Parallel processing check
  nCores <- detectCores()
  if(!is.null(parallel) && parallel && nCores < 4)  {
    if(verbose)
      warning("Not enough cores for parallel processing, running chains sequentially.")
    parallel <- FALSE
  }
  if(is.null(parallel))
    parallel <- nCores > 3

  # Data checks
  if(!all(is.finite(y)))
    stop("The input data include NA or Inf.")
  if(length(unique(y)) < 2 &&      # sd(y) will be 0 or NA; ok if priors specified.
        (is.null(priors) ||
          is.null(priors$muSD) ||
          is.null(priors$sigmaMode) ||
          is.null(priors$sigmaSD)))
  stop("If priors are not specified, data must include at least 2 (non-equal) values.")

  # Prior checks:
  if(!is.null(priors))  {
    if(!is.list(priors)) {
        stop("'priors' must be a list (or NULL).")
    }
    nameOK <- names(priors) %in%
          c("muM", "muSD", "sigmaMode", "sigmaSD")
    if(!all(nameOK))
      stop("Invalid items in prior specification: ",
          paste(sQuote(names(priors)[!nameOK]), collapse=", "))
    if(!all(sapply(priors, is.numeric)))
      stop("All items in 'priors' must be numeric.")
    if(!is.null(priors$muSD) && priors$muSD <= 0)
      stop("muSD must be > 0")
  }
  if(is.null(rnd.seed))
    rnd.seed <- floor(runif(1,1,10000))

  # THE PRIORS
  if(is.null(priors)) {   # use the safe prior specification
    dataForJAGS <- list(
      muM = mean(y) ,
      muP = 0.000001 * 1/sd(y)^2 ,
      sigmaLow = sd(y) / 1000 ,
      sigmaHigh = sd(y) * 1000
    )
  } else {    # use gamma priors
    priors0 <- list(  # default priors
      muM = mean(y) ,
      muSD = sd(y)*5 ,
      sigmaMode = sd(y),
      sigmaSD = sd(y)*5)
    priors0 <- modifyList(priors0, priors)  # user's priors take prior-ity (duh!!)
    if(mean(y) > priors0$muM + priors0$muSD || mean(y) < priors0$muM - priors0$muSD)
      warning("Sample mean is outside the prior range muM \u00B1 muSD.")
    # Convert to Shape/Rate
    sigmaShRa <- gammaShRaFromModeSD(mode=priors0$sigmaMode, sd=priors0$sigmaSD)
    dataForJAGS <- list(
      muM = priors0$muM,
      muP = 1/priors0$muSD^2,  # convert SD to precision
      Sh = sigmaShRa$shape,
      Ra = sigmaShRa$rate)
  }

  # THE MODEL.
  modelFile <- file.path(tempdir(), "BESTmodel.txt")
  if(is.null(priors)) {  # use old broad priors
    modelString = "
    model {
      for ( i in 1:Ntotal ) {
        y[i] ~ dnorm(mu, tau)
      }
      mu ~ dnorm(muM, muP)
      tau <- pow(sigma, -2)
      sigma ~ dunif(sigmaLow, sigmaHigh)
    }
    " # close quote for modelString
  } else {    # use gamma priors
    modelString = "
    model {
      for ( i in 1:Ntotal ) {
        y[i] ~ dnorm(mu, tau)
      }
      mu ~ dnorm(muM, muP)
      tau <- pow(sigma, -2)
      sigma ~ dgamma(Sh, Ra)
    }
    " # close quote for modelString
  }
  # Write out modelString to a text file
  writeLines( modelString , con=modelFile )

  # THE DATA.
  # dataForJAGS already has the priors, add the data:
  if(!doPriorsOnly)
    dataForJAGS$y <- y
  dataForJAGS$Ntotal <- length(y)

  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  initsList0 <- list(mu=mean(y), sigma=sd(y), .RNG.seed=rnd.seed)
  initsList <- list(
                c(initsList0, .RNG.name="base::Wichmann-Hill"),
                c(initsList0, .RNG.name="base::Marsaglia-Multicarry"),
                c(initsList0, .RNG.name="base::Super-Duper") )

  # RUN THE CHAINS
  codaSamples <- jags.basic(
    data = dataForJAGS,
    inits = initsList,
    parameters.to.save = c( "mu" , "sigma" ),     # The parameters to be monitored
    model.file = modelFile,
    n.chains = 3,    # Do not change this without also changing initsList.
    n.adapt = 500,
    n.iter = ceiling( ( numSavedSteps * thinSteps) / 3  + burnInSteps ),
    n.burnin = burnInSteps,
    n.thin = thinSteps,
    modules = NULL,
    parallel = parallel,
    DIC = FALSE,
    seed = rnd.seed,
    verbose = verbose)

  Rhat <- try(gelman.diag(codaSamples, autoburnin=FALSE)$psrf[, 1], silent=TRUE)
  if(inherits(Rhat, "try-error") || !all(is.finite(Rhat)))
    Rhat <- NULL

  out <- as.Bwiqid(codaSamples,
      header = "Model fitted in JAGS with jagsUI::jags.basic",
      defaultPlot = names(codaSamples)[1])
  attr(out, "call") <- match.call()
  attr(out, "n.chains") <- 3
  attr(out, "n.eff") <- effectiveSize(codaSamples)
  attr(out, "Rhat") <- Rhat
  attr(out, "doPriorsOnly") <- doPriorsOnly
  if(!is.null(priors))
    attr(out, "priors") <- priors0

  attr(out, "timetaken") <- Sys.time() - startTime
  return(out)
}
