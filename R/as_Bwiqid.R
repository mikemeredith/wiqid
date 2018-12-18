

# This file has the S3 generic 'as.Bwiqid' function and a series of methods.

as.Bwiqid <- function(object, ...) UseMethod("as.Bwiqid")

as.Bwiqid.default <- function(object, ...) {
    stop(paste("No applicable method for class", class(object)))
}
# ...................................................................

# Class Bwiqid (catches errors and allows header, defaultPlot to be changed)
as.Bwiqid.Bwiqid <- function(object, header, defaultPlot, ...) {
  out <- object
  if(!missing(header))
    attr(out, "header") <- header
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}
# ................................................................

# Class data.frame
as.Bwiqid.data.frame <- function(object, header, defaultPlot, n.chains=1,
    Rhat=TRUE, n.eff=TRUE, ...) {
  npar <- length(object)
  out <- object
  class(out) <- c("Bwiqid", class(out))
  if(!missing(header))
    attr(out, "header") <- header
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  attr(out, "n.chains") <- n.chains

  if(n.chains > 1) {
    if(is.logical(Rhat)) {
      if(Rhat)
        attr(out, "Rhat") <- simpleRhat(out, n.chains=n.chains)
    } else if(is.numeric(Rhat) && length(Rhat) == npar) {
      attr(out, "Rhat") <- Rhat
    } else {
      warning("'Rhat' must be logical, or a vector with a value for each column.",
        call.=FALSE)
    }
  }
  if(is.logical(n.eff)) {
    if(n.eff)
      attr(out, "n.eff") <- safeNeff(out)
  } else if(is.numeric(n.eff) && length(n.eff) == npar) {
    attr(out, "n.eff") <- n.eff
  } else {
    warning("'n.eff' must be logical, or a vector with a value for each column.",
      call.=FALSE)
  }
  return(out)
}
# .......................................................................

# Class mcmc.list from (inter alia) rjags package
as.Bwiqid.mcmc.list <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(as.matrix(object))
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(!missing(header))
    attr(out, "header") <- header
  attr(out, "n.chains") <- length(object)
  attr(out, "n.eff") <- safeNeff(out)
  if(length(object) > 1) {
    attr(out, "Rhat") <- simpleRhat(out, n.chains=length(object))
  }
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Class mcmc
as.Bwiqid.mcmc <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(as.matrix(object))
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(!missing(header))
    attr(out, "header") <- header
  attr(out, "n.chains") <- 1
  attr(out, "n.eff") <- safeNeff(out)
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Class bugs from R2WinBUGS package and R2OpenBUGS
as.Bwiqid.bugs <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(object$sims.matrix)
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(missing(header))
    header <- paste("Model fitted in", object$program)
  attr(out, "header") <- header
  attr(out, "n.chains") <- object$n.chains
  if(object$n.chains > 1)
    attr(out, "Rhat") <- simpleRhat(out, n.chains=object$n.chains)
  attr(out, "n.eff") <- safeNeff(out)
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Class rjags from R2jags package
as.Bwiqid.rjags <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(object$BUGSoutput$sims.matrix)
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(missing(header))
    header <- "Model fitted in JAGS with R2jags"
  attr(out, "header") <- header
  attr(out, "n.chains") <- object$BUGSoutput$n.chains
  if(object$BUGSoutput$n.chains > 1)
    attr(out, "Rhat") <- simpleRhat(out, n.chains=object$BUGSoutput$n.chains)
  attr(out, "n.eff") <- safeNeff(out)
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Class jagsUI from jagsUI package
as.Bwiqid.jagsUI <- function(object, header, defaultPlot, ...) {
  stopifnot(class(object$samples) == "mcmc.list")
  out <- as.data.frame(as.matrix(object$samples))
  names(out) <- fixNames(names(out))
  n.chains <- length(object$samples)
  class(out) <- c("Bwiqid", class(out))
  if(missing(header))
    header <- "Model fitted in JAGS with jagsUI"
  attr(out, "header") <- header
  attr(out, "n.chains") <- n.chains
  attr(out, "n.eff") <- safeNeff(out)
  if(n.chains > 1)
    attr(out, "Rhat") <- simpleRhat(out, n.chains=n.chains)
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  attr(out, "timetaken") <- as.difftime(object$mcmc.info$elapsed.mins, units="mins")
  return(out)
}
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Class runjags from runjags package
as.Bwiqid.runjags <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(as.matrix(object$mcmc))
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(missing(header))
    header <- "Model fitted in JAGS with runjags"
  attr(out, "header") <- header
  attr(out, "n.chains") <- length(object$mcmc)
  attr(out, "n.eff") <- safeNeff(out)
  # attr(out, "MCerror") <- as.data.frame(unclass(object$mcse))$mcse
  if(length(object$mcmc) > 1)
    attr(out, "Rhat") <- simpleRhat(out, n.chains=length(object$mcmc))
  if(!missing(defaultPlot))
    attr(out, "defaultPlot") <- defaultPlot
  attr(out, "timetaken") <- object$timetaken
  return(out)
}
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# An error-catching wrapper for coda::effectiveSize
safeNeff <- function(x) {
  # x is a data frame or matrix with a column for each parameter
  safe1 <- function(v) {
    tmp <- try(coda::effectiveSize(v), silent=TRUE)
    if(inherits(tmp, "try-error"))
      return(NA)
    return(tmp)
  }
  apply(x, 2, safe1)
}


