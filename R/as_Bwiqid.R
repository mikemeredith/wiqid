

# This file has the S3 generic 'as.Bwiqid' function and a series of methods.

as.Bwiqid <- function(object, ...) UseMethod("as.Bwiqid")

as.Bwiqid.default <- function(object, ...) {
    stop(paste("No applicable method for class", class(object)))
}

# Class Bwiqid (catches errors and allows header, defaultPlot to be changed)
as.Bwiqid.Bwiqid <- function(object, header, defaultPlot, ...) {
  out <- object
  if(!missing(header))
    attr(out, "header") <- header
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}

# Class data.frame 
as.Bwiqid.data.frame <- function(object, header, defaultPlot, ...) {
  out <- object
  class(out) <- c("Bwiqid", class(out))
  if(!missing(header))
    attr(out, "header") <- header
  attr(out, "n.chains") <- 1
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}

# Class mcmc.list from (inter alia) rjags package
as.Bwiqid.mcmc.list <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(as.matrix(object))
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(!missing(header))
    attr(out, "header") <- header
  attr(out, "n.chains") <- length(object)
  attr(out, "n.eff") <- effectiveSize(object)
  if(length(object) > 1)
    attr(out, "Rhat") <- gelman.diag(object, autoburnin=FALSE, multivariate=FALSE)$psrf[, 1]
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}
# Class mcmc
as.Bwiqid.mcmc <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(as.matrix(object))
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(!missing(header))
    attr(out, "header") <- header
  attr(out, "n.chains") <- 1
  attr(out, "n.eff") <- effectiveSize(object)
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}

# Class bugs from R2WinBUGS package and R2OpenBUGS
as.Bwiqid.bugs <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(object$sims.matrix)
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(missing(header))
    header <- paste("Model fitted in", object$program)
  attr(out, "header") <- header
  attr(out, "n.chains") <- object$n.chains
  if(object$n.chains > 1) {
    if(ncol(out) > 1) {
      attr(out, "n.eff") <- object$summary[, 'n.eff']
      attr(out, "Rhat") <- object$summary[, 'Rhat']
    } else {
      attr(out, "n.eff") <- object$summary['n.eff'] # summary is a vector
      attr(out, "Rhat") <- object$summary['Rhat']
    }
  }
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}

# Class rjags from R2jags package
as.Bwiqid.rjags <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(object$BUGSoutput$sims.matrix)
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(missing(header))
    header <- "Model fitted in JAGS with R2jags"
  attr(out, "header") <- header
  attr(out, "n.chains") <- object$BUGSoutput$n.chains
  if(object$BUGSoutput$n.chains > 1) {
    attr(out, "n.eff") <- object$BUGSoutput$summary[, 'n.eff'] ## not calculated if 1 chain
    attr(out, "Rhat") <- object$BUGSoutput$summary[, 'Rhat']
  }
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}

# Class jagsUI from jagsUI package
as.Bwiqid.jagsUI <- function(object, header, defaultPlot, ...) {
  stopifnot(class(object$samples) == "mcmc.list")
  out <- as.data.frame(as.matrix(object$samples))
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(missing(header))
    header <- "Model fitted in JAGS with jagsUI"
  attr(out, "header") <- header
  attr(out, "n.chains") <- length(object$samples)
  if(length(object$samples) > 1) {
    attr(out, "n.eff") <- unlist(object$n.eff)
    if(length(object) > 1)
      attr(out, "Rhat") <- unlist(object$Rhat)
  }
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  attr(out, "timetaken") <- as.difftime(object$mcmc.info$elapsed.mins, units="mins")
  return(out)
}

# Class runjags from runjags package
as.Bwiqid.runjags <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(as.matrix(object$mcmc))
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(missing(header))
    header <- "Model fitted in JAGS with runjags"
  attr(out, "header") <- header
  attr(out, "n.chains") <- length(object$mcmc)
  attr(out, "n.eff") <- as.data.frame(unclass(object$mcse))$sseff
  attr(out, "MCerror") <- as.data.frame(unclass(object$mcse))$mcse
  if(length(object$mcmc) > 1)
    attr(out, "Rhat") <- object$psrf$psrf[, 1]
  if(!missing(defaultPlot))
    attr(out, "defaultPlot") <- defaultPlot
  attr(out, "timetaken") <- object$timetaken
  return(out)
}




