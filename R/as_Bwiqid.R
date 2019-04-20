

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
  names(out) <- make.names(names(out), unique=TRUE)
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
  names(out) <- make.names(names(out), unique=TRUE)
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

# Core function for classes 'bugs', 'rjags' and 'jagsUI'
#   handles MCMC in matrix plus summary of results already calculated
asBwiqidCore <- function(mcmat, summary, nChains) {
  out <- as.data.frame(mcmat)
  BUGSnames <- colnames(out)
  if(!all(BUGSnames == rownames(summary)))
    stop("'summary' names do not match MCMC names.")
  names(out) <- make.names(names(out), unique=TRUE)
  # Try to recover pre-calculated Rhat and n.eff, else roll our own
  Rhat <- try(summary[, 'Rhat'], silent=TRUE)
  if(inherits(Rhat, "try-error")) {
    if(nChains > 1) {
      Rhat <- simpleRhat(out, n.chains=nChains)
    } else {
      Rhat <- NULL
    }
  }
  n.eff <- try(summary[, 'n.eff'], silent=TRUE)
  if(inherits(n.eff, "try-error")) {
    n.eff <- safeNeff(out)
  }
  # Main attributes
  class(out) <- c("Bwiqid", class(out))
  attr(out, "BUGSnames") <- BUGSnames
  attr(out, "n.chains") <- nChains
  if(!is.null(Rhat))
    attr(out, "Rhat") <- Rhat
  attr(out, "n.eff") <- n.eff
  return(out)
} # ''''''''''''''''''''''''''''''''''''''''''''''


# Class bugs from R2WinBUGS package and R2OpenBUGS
as.Bwiqid.bugs <- function(object, header, defaultPlot, ...) {
  out <- asBwiqidCore(object$sims.matrix, object$summary, object$n.chains)
  if(missing(header))
    header <- paste("Model fitted in", object$program)
  attr(out, "header") <- header
  attr(out, "modelFile") <- object$model.file
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Class rjags from R2jags package
as.Bwiqid.rjags <- function(object, header, defaultPlot, ...) {
  object <- object$BUGSoutput
  out <- asBwiqidCore(object$sims.matrix, object$summary, object$n.chains)
  if(missing(header))
    header <- "Model fitted in JAGS with R2jags"
  attr(out, "header") <- header
  attr(out, "modelFile") <- object$model.file
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  return(out)
}
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Class jagsUI from jagsUI package ## updated 2019-04-17
as.Bwiqid.jagsUI <- function(object, header, defaultPlot, ...) {
  stopifnot(class(object$samples) == "mcmc.list")
  out <- asBwiqidCore(as.matrix(object$samples), object$summary, object$mcmc.info$n.chains)
  if(missing(header))
    header <- "Model fitted in JAGS with jagsUI"
  attr(out, "header") <- header
  attr(out, "modelFile") <- object$modfile
  if(!missing("defaultPlot"))
    attr(out, "defaultPlot") <- defaultPlot
  attr(out, "timetaken") <- as.difftime(object$mcmc.info$elapsed.mins, units="mins")
  return(out)
}
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Class runjags from runjags package
as.Bwiqid.runjags <- function(object, header, defaultPlot, ...) {
  out <- as.data.frame(as.matrix(object$mcmc))
  names(out) <- make.names(names(out), unique=TRUE)
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

# safeNeff is now in the simpleRhat.R file
