# diagPlot, tracePlot, densityPlot and acfPlot functions for class Bwiqid, ie. MCMC output

# Function to do multiple trace and density plots
diagPlot <- function(object, which, howMany, ask=TRUE, maxRow=4, RhatBad=1.05, ...) {
  if(!inherits(object, "Bwiqid"))
    stop("object is not a valid Bwiqid object")
  npars <- ncol(object)
  n.chains <- attr(object, "n.chains")
  if(is.null(n.chains) || nrow(object) %% n.chains != 0) {
    warning("Invalid number of chains, treating data as a single chain")
    n.chains <- 1
  }
  n.iter <- nrow(object) / n.chains
  Rhat <- attr(object, "Rhat")
  if(is.null(Rhat))
    Rhat <- rep(NA, npars)
  Rhat <- round(Rhat, 2)

  n.eff <- attr(object, "n.eff")
  if(is.null(n.eff))
    n.eff <- rep(NA, npars)
  n.eff <- round(n.eff)

  if(!missing(howMany) && abs(howMany) < n.iter) {
    if(howMany > 0) {
      object <- window(object, end = howMany)
    } else {
      object <- window(object, start = n.iter + howMany)
    }
  }

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
    n.eff <- n.eff[which]
  }

  old.par <- par(mar = c(2,2,2, 0)+0.1, oma=c(1,1,1,1), "mfrow")
    on.exit(par(old.par))

  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  defaultArgs <- list(xlab="Iterations", ylab="Density",
    type='l', lty=1)
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
    title(main = paste0(names(object[i]), ", Rhat = ", Rhat[i]),
      line=1, adj=1, col.main=1 + redFlag)
    if(redFlag)
      box(col=2, lwd=2)
    # do density plot
    density0(mat, useArgsD)
    title(main = paste0("n.eff = ", n.eff[i]),
        line=1, adj=0, col.main=1 + redFlag)
    if(redFlag)
      box(col=2, lwd=2)
  }
}
# ..........................................................

# Helper function to do density plot (or histogram) for 1 parameter
# mat : matrix of MCMC output with 1 column per chain
# plotArgs : list with plotting parameters
density0 <- function(mat, plotArgs, ...)  {

  bw <- bw.nrd0(mat)
  # histogram or density plot?
  unik <- unique.default(mat)
  if(length(unik) == 1) {
    plot(1, 1, type = "n", ann = FALSE, axes = FALSE)
    text(1,1, paste("All values are the same:\n", signif(unik, 4)))
  } else if(all(mat %% 1 == 0) && all(mat >= 0) && diff(range(mat)) < 50) {
    # "histogram"
    t1 <- apply(mat+1, 2, tabulate, nbins=max(mat)+1)/nrow(mat) # +1 cos tabulate ignores 0
    if(min(mat) > 0)
      t1 <- t1[-(1:min(mat)), ]
    ymax <- apply(t1, 1, max)
    xx <- min(mat):max(mat)
    xlim <- c(min(mat)-0.5, max(mat)+0.5)
    plot(xx, ymax, type='h', col='grey', xlim=xlim)
    abline(h=0)
    abline(v=colMeans(mat), col=1:ncol(mat), lwd=2, lty=3)
    segments(x0 = rep(xx - 0.4, ncol(mat)),
             y0 = t1,
             x1 = rep(xx + 0.4, ncol(mat)),
             y1 = t1,
             col = col(t1))
  } else if (bw == 0){
    plot(1, 1, type = "n", ann = FALSE, axes = FALSE)
    text(1,1, "Bandwidth is zero.")
  } else {
    # density plot
    # deal with folding for probability and non-negative values
    # meanMat <- mean(mat) # do this before folding
    meanMat <- colMeans(mat) # do this before folding
    # use these values if folding is not needed:
    from <- min(mat) - 3*bw
    to <- max(mat) + 3*bw
    mult <- 1
    xx <- mat

    if (min(mat) >= 0 && min(mat) < 2 * bw) {  # it's non-negative
      from <- 0
      xx <- rbind(mat, -mat)
      mult <- 2
    }
    if (min(mat) >= 0 && max(mat) <= 1 &&
          (min(mat) < 2 * bw || 1 - max(mat) < 2 * bw)) { # it's a probability
      to <- min(to, 1)
      xx <- rbind(mat, -mat, 2-mat)
      mult <- 3
    }

    # fit density to each column
    n <- 512
    dens <- apply(xx, 2, function(x) density(x, bw=bw, from=from, to=to, n=n)$y)

    plotArgs$x <- seq(from, to, length.out=n)
    plotArgs$y <- dens * mult
    do.call(matplot, plotArgs)
    abline(h=0, col='grey')
    abline(v=meanMat, col=1:ncol(mat), lwd=2, lty=3)
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

crosscorrPlot <- function(object, which, ...)  {
  if(!inherits(object, "Bwiqid"))
    stop("object is not a valid Bwiqid object")

  if(!missing(which)) {
    if(is.character(which))
      which <- pmatch(which, names(object))
    which <- which[!is.na(which)]
    if(length(which) == 0)
      stop("No parameters selected.")
    if(any(which > ncol(object)))
      stop("Invalid parameters selected.")
    object <- object[which]
  }

  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  defaultArgs <- list(method='color', type='lower', cl.pos='r', tl.col='black',
      col=topo.colors(10), addgrid.col='black')
  useArgs <- modifyList(defaultArgs, dots)

  crosscorr <- suppressWarnings(cor(object))
  useArgs$corr <- crosscorr
  do.call(corrplot::corrplot, useArgs)

  return(invisible(crosscorr))
}


