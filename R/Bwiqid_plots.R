# diagPlot, tracePlot, densityPlot and acfPlot functions for class Bwiqid, ie. MCMC output

# Function to do multiple trace and density plots
diagPlot <- function(object, which, ask=TRUE, maxRow=4, RhatBad=1.05, ...) {
  if(!inherits(object, "Bwiqid"))
    stop("object is not a valid Bwiqid object")
  npars <- ncol(object)
  n.chains <- attr(object, "n.chains")
  if(is.null(n.chains) || nrow(object) %% n.chains != 0) {
    warning("Invalid number of chains, treating data as a single chain")
    n.chains <- 1
  }
  Rhat <- attr(object, "Rhat")
  if(is.null(Rhat))
    Rhat <- rep(NA, npars)
  Rhat <- round(Rhat, 2)

  if(!missing(which)) {
    if(is.character(which))
      which <- pmatch(which, names(object))
    which <- which[!is.na(which)]
    if(length(which) == 0)
      stop("No parameters selected.")
    if(any(which > npars))
      stop("Invalid parameters selected.")
    object <- object[which]
    npars <- ncol(object)
    Rhat <- Rhat[which]
  }

  old.par <- par(mar = c(2,2,2, 0)+0.1, oma=c(1,1,1,1), "mfrow")
    on.exit(par(old.par))

  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  # Recommended colours for colour-blind people:
  cbCol <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7")
  defaultArgs <- list(xlab="Iterations", ylab="Density",
    type='l', lty=1, col=cbCol)
  useArgsT <- useArgsD <- modifyList(defaultArgs, dots)

  if(npars > maxRow) {
    old.ask <- devAskNewPage(ask)
    on.exit(devAskNewPage(old.ask), add=TRUE)
  }
  nrows <- min(npars, maxRow)
  layout(matrix(1:(nrows*2), ncol=2, byrow=TRUE), widths=2:1)
  for(i in 1:npars) {
    redFlag <- !is.na(Rhat[i]) && Rhat[i] > RhatBad
    mat <- matrix(object[, i], ncol=n.chains)
    # do trace plot
    useArgsT$ylab <- names(object[i])
    useArgsT$y <- mat
    do.call(matplot, useArgsT)
    abline(h=mean(mat))
    title( main = paste0(names(object[i]), ", Rhat = ", Rhat[i]),
      line=1, adj=1, col.main=1 + redFlag)
    if(redFlag)
      box(col='red', lwd=2)
    # do density plot
    density0(mat, useArgsD)
    if(redFlag)
      box(col='red', lwd=2)
  }
}
# ..........................................................

# Helper function to do density plot (or histogram) for 1 parameter
# mat : matrix of MCMC output with 1 column per chain
# plotArgs : list with plotting parameters
density0 <- function(mat, plotArgs, ...)  {

  bw <- bw.nrd0(mat)
  # histogram or density plot?
  if(max(abs(mat - floor(mat))) == 0 || bw == 0 || length(unique(mat)) == 1) {
    breaks <- seq(min(mat)-1, max(mat), by=1)
    breaks2 <- hist(mat, plot=FALSE)$breaks
    if(length(breaks2) < length(breaks)) {
      breaks <- breaks2
      breaks[1] <- min(mat) - 1
      breaks[length(breaks)] <- max(mat)
    }
    nb <- length(breaks)
    hgt <- apply(mat, 2, function(x) hist(x, breaks=breaks, plot=FALSE)$density)
    if(!is.matrix(hgt))
      hgt <- matrix(hgt, nrow=1)
    plot(breaks, rep(0, nb), type='n', ylim=c(0, max(hgt)))
    for(i in 1:ncol(mat))
      rect(breaks[-nb], rep(0, nb-1), breaks[-1], hgt[,i], border=plotArgs$col[i])
    abline(h=0)
  } else {
    # deal with folding for probability and non-negative values
    meanMat <- mean(mat) # do this before folding
    if (min(mat) >= 0 && min(mat) < 2 * bw) {
      if (max(mat) <= 1 && 1 - max(mat) < 2 * bw) { # it's a probability
        constrain <- 2
        mat <- rbind(mat, -mat, 2-mat)
      } else {                                      # it's non-negative
        constrain <- 1
        mat <- rbind(mat, -mat)
      }
    } else {
      constrain <- 0
    }

    # All columns must have same x coordinates; see ?density for details
    from <- min(mat) - 3*bw
    to <- max(mat) + 3*bw
    n <- 512
    x <- seq(from, to, length.out=n)
    dens <- apply(mat, 2, function(x) density(x, bw=bw, from=from, to=to, n=n)$y)

    if(constrain == 2) {         # probability
      dens <- dens[x >=0 & x <=1, ] * 3
      x <- x[x >=0 & x <=1]
    } else if(constrain == 1) {  # non-negative
      dens <- dens[x >=0, ] * 2
      x <- x[x >=0]
    }                            # if(constrain == 0) do nothing

    plotArgs$x <- x
    plotArgs$y <- dens
    do.call(matplot, plotArgs)
    abline(v=meanMat)
  }
}
# ..........................................................

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
# ..........................................................

# FIXME: use density0
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

