# trace plot, density plot and acf plot functions for class Bwiqid, ie. MCMC output

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
  # Recommended colours for colour-blind people:
  cbCol <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7")
  defaultArgs <- list(xlab="Iterations", type='l', lty=1, col=cbCol)
  useArgs <- modifyList(defaultArgs, dots)
  mainStem <- useArgs$main
  
  for( i in 1:ncol(object)) {
    mat <- matrix(object[, i], ncol=n.chains)
    useArgs$ylab <- names(object[i])
    useArgs$y <- mat
    do.call(matplot, useArgs)
    abline(h=mean(mat))
  }
    
  return(invisible(NULL))
}

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
  # Recommended colours for colour-blind people:
  cbCol <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7")
  defaultArgs <- list(ylab="Density", type='l', lty=1, col=cbCol)
  useArgs <- modifyList(defaultArgs, dots)
  # mainStem <- useArgs$main
  
  for( i in 1:ncol(object)) {
    mat <- matrix(object[, i], ncol=n.chains)
    # bw <- bw.SJ(mat) # bw.SJ is slow and often fails: "sample is too sparse..."
    bw <- bw.nrd0(mat)
    from <- min(mat) - 3*bw
    to <- max(mat) + 3*bw
    dens <- apply(mat, 2, function(x) density(x, bw=bw, n=128, from=from, to=to)$y)
    useArgs$xlab <- names(object[i])
    useArgs$x <- seq(from, to, length=128)
    useArgs$y <- dens
    do.call(matplot, useArgs)
    abline(v=mean(mat))
  }
    
  return(invisible(NULL))
}

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
  # Recommended colours for colour-blind people:
  cbCol <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7")
  defaultArgs <- list(ylab="ACF", xlab="Lag", type='h', lty=1, col=cbCol)
  useArgs <- modifyList(defaultArgs, dots)
  # mainStem <- useArgs$main
  
  for( i in 1:ncol(object)) {
    mat <- matrix(object[, i], ncol=n.chains)
    acor <- apply(mat, 2, function(x) acf(x, lag.max=lag.max, plot=FALSE)$acf)
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
    
  return(invisible(NULL))
}

