

# This file has the S3 generic 'as.Bwiqid' function and a series of methods.

as.Bwiqid <- function(object, ...) UseMethod("as.Bwiqid")

as.Bwiqid.default <- function(object, ...) {
    stop(paste("No applicable method for class", class(object)))
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
  attr(out, "Rhat") <- gelman.diag(object, multivariate=FALSE)$psrf[, 1]
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}

# Class bugs from R2WinBUGS package (and R2OpenBUGS too?)
as.Bwiqid.bugs <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(object$sims.matrix)
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(missing(header))
    header <- paste("Model fitted in", object$program)
  attr(out, "header") <- header
  attr(out, "n.chains") <- object$n.chains
  attr(out, "n.eff") <- object$summary[, 'n.eff']
  attr(out, "Rhat") <- object$summary[, 'Rhat']
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
  attr(out, "n.eff") <- object$BUGSoutput$summary[, 'n.eff']
  attr(out, "Rhat") <- object$BUGSoutput$summary[, 'Rhat']
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}

# Class runjags from runjags package
as.Bwiqid.runjags <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(as.matrix(as.mcmc.list(object)))
  names(out) <- fixNames(names(out))
  class(out) <- c("Bwiqid", class(out))
  if(missing(header))
    header <- "Model fitted in JAGS with runjags"
  attr(out, "header") <- header
  attr(out, "n.chains") <- length(object$mcmc)
  attr(out, "n.eff") <- as.data.frame(unclass(object$mcse))$sseff
  attr(out, "MCerror") <- as.data.frame(unclass(object$mcse))$mcse
  attr(out, "Rhat") <- object$psrf$psrf[, 1]
  if(!missing(defaultPlot))
    attr(out, "defaultPlot") <- defaultPlot
  attr(out, "timetaken") <- object$timetaken
  return(out)
}




