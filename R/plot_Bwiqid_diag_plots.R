# tracePlot, densityPlot and acfPlot functions for class Bwiqid, ie. MCMC output

# diagPlot and density0 moved to plot_diagPlot.R
# crosscorrPlot moved plot_crosscorrPlot.R

tracePlot <- function(object, ask=TRUE, ...)  {
  if(!inherits(object, "Bwiqid"))
    stop("object is not a valid Bwiqid object")
  n.chains <- attr(object, "n.chains")
  if(is.null(n.chains) || nrow(object) %% n.chains != 0) {
    warning("Invalid number of chains, treating data as a single chain")
    n.chains <- 1
  }
  old.ask <- devAskNewPage(ask) ; on.exit(devAskNewPage(old.ask))
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  defaultArgs <- list(xlab="Iterations", type='l', lty=1)
  useArgs <- modifyList(defaultArgs, dots)
  mainStem <- useArgs$main

  for( i in 1:ncol(object)) {
    mat <- matrix(object[, i], ncol=n.chains)
    useArgs$ylab <- names(object[i])
    useArgs$y <- mat
    useArgs$main <- names(object)[i]
    do.call(matplot, useArgs)
    abline(h=mean(mat))
  }

  return(invisible(NULL))
}
# ..........................................................

densityPlot <- function(object, ask=TRUE, ...)  {
  if(!inherits(object, "Bwiqid"))
    stop("object is not a valid Bwiqid object")
  n.chains <- attr(object, "n.chains")
  if(is.null(n.chains) || nrow(object) %% n.chains != 0) {
    warning("Invalid number of chains, treating data as a single chain")
    n.chains <- 1
  }
  old.ask <- devAskNewPage(ask) ; on.exit(devAskNewPage(old.ask))
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  defaultArgs <- list(ylab="Density", type='l', lty=1, xlab="")
  useArgs <- modifyList(defaultArgs, dots)
  # mainStem <- useArgs$main

  for( i in 1:ncol(object)) {
    mat <- matrix(object[, i], ncol=n.chains)
    useArgs$main <- names(object)[i]
    density0(mat, useArgs)
  }

  return(invisible(NULL))
}
# ..........................................................

acfPlot <- function(object, lag.max=NULL, ask=TRUE, ...)  {
  if(!inherits(object, "Bwiqid"))
    stop("object is not a valid Bwiqid object")
  n.chains <- attr(object, "n.chains")
  if(is.null(n.chains) || nrow(object) %% n.chains != 0) {
    warning("Invalid number of chains, treating data as a single chain")
    n.chains <- 1
  }
  old.ask <- devAskNewPage(ask) ; on.exit(devAskNewPage(old.ask))
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  defaultArgs <- list(ylab="ACF", xlab="Lag", type='h', lty=1)
  useArgs <- modifyList(defaultArgs, dots)
  # mainStem <- useArgs$main

  for( i in 1:ncol(object)) {
    mat <- matrix(object[, i], ncol=n.chains)
    acor <- apply(mat, 2, function(x) acf(x, lag.max=lag.max, plot=FALSE)$acf)
    if(any(is.na(acor))) {
      plot(1, 1, type = "n", ann = FALSE, axes = FALSE)
      text(1,1, "No ACF calculated.")
      title(main=names(object[i]))
    } else {
      lags <- 0:(nrow(acor) - 1)
      if (n.chains > 1) {
        jitt <- seq(-0.2, 0.2, length=n.chains)
      } else {
        jitt <- 0
      }
      lags.mat <- outer(lags, jitt, "+")
      useArgs$main <- names(object[i])
      useArgs$x <- lags.mat
      useArgs$y <- acor
      do.call(matplot, useArgs)
      abline(h=0)
    }
  }

  return(invisible(NULL))
}
