# Just diagPlot, liberated to work with other classes than Bwiqid

# Do helper functions to produce 3D array with MCMC output
#  itersPerChain x chains x parameters
#  plus matrix with Rhat and MCEpc

# This works with most R-to-BUGS/JAGS output
getArrA <- function(x, summary, nChains) {
  if(is.array(x) && length(dim(x)) == 3) {
    mcmc3d <- x
    if(!is.null(nChains) && dim(x)[2] != nChains)
      stop("nChains does not match array dimensions.")
    mcmat <- matrix(x, dim(x)[1] * dim(x)[2], dim(x)[3])
  } else {
    mcmat <- as.matrix(x)
    if(is.null(nChains) || nrow(mcmat) %% nChains != 0) {
      warning("Invalid number of chains, treating data as a single chain")
      nChains <- 1
    }
    mcmc3d <- array(mcmat, c(nrow(mcmat)/nChains, nChains, ncol(mcmat)))
    dimnames(mcmc3d) <- list(iter=NULL, chain=1:nChains, param=colnames(mcmat))
  }
  if(!all(dimnames(mcmc3d)[[3]] == rownames(summary)))
    stop("'summary' names do not match MCMC names.")
  # Try to recover pre-calculated Rhat and n.eff, else roll our own
  Rhat <- try(summary[, 'Rhat'], silent=TRUE)
  if(inherits(Rhat, "try-error")) {
    if(nChains > 1) {
      Rhat <- simpleRhat(mcmat, n.chains=nChains)
    } else {
      Rhat <- NA
    }
  }
  # n.eff <- try(summary[, 'n.eff'], silent=TRUE)
  # if(inherits(n.eff, "try-error")) {
    # n.eff <- safeNeff(mcmat)
  # }
  MCerror <- attr(x, "MCerror")
  if(is.null(MCerror))
    MCerror <- getMCerror(mcmat, n.chains=nChains)
  SD <- apply(mcmat, 2, sd)
  MCEpc <- MCerror / SD * 100

  return(list(mcmc3d = mcmc3d,
              stats = cbind(Rhat=Rhat, MCEpc=MCEpc)))
} # ''''''''''''''''''''''''''''''''''''''''''''''

# This works with mcmc.list, Bwiqid and data frames
getArrB <- function(x, nChains=1) {
  mcmat <- as.matrix(x)
  if(!is.numeric(mcmat))
    stop("Sorry, can't find the numbers in your input.", call.=FALSE)
  names <- colnames(mcmat)
  # Look for attributes, else roll our own
  nc <- attr(x, "n.chains")
  if(!is.null(nc))
    nChains <- nc
  if(is.null(nChains) || nrow(mcmat) %% nChains != 0) {
    warning("Invalid number of chains, treating data as a single chain")
    nChains <- 1
  }
  Rhat <- attr(x, "Rhat")
  if(is.null(Rhat) && nChains > 1)
    Rhat <- simpleRhat(mcmat, n.chains=nChains)
  if(is.null(Rhat))
    Rhat <- NA
  # n.eff <- attr(x, "n.eff")
  # if(is.null(n.eff))
    # n.eff <- safeNeff(mcmat)
  MCerror <- attr(x, "MCerror")
  if(is.null(MCerror))
    MCerror <- getMCerror(mcmat, n.chains=nChains)
  SD <- apply(mcmat, 2, sd)
  MCEpc <- MCerror / SD * 100

  dim(mcmat) <- c(nrow(mcmat)/nChains, nChains, ncol(mcmat))
  dimnames(mcmat) <- list(iter=NULL, chain=1:nChains, param=names)

  return(list(mcmc3d = mcmat,
              stats = cbind(Rhat=Rhat, MCEpc=MCEpc)))
} # ''''''''''''''''''''''''''''''''''''''''''''''


# Function to do multiple trace and density plots
diagPlot <- function(x, params=NULL, howMany, chains, ask=NULL,
  maxRow=4, RhatBad=1.05, ...) {

  if(is.null(ask))
    ask <- dev.interactive(orNone=TRUE)

  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  # catch the old 'which' argument
  if(!is.null(dots$which) && is.null(params)) {
    warning("Argument 'which' is deprecated, please use 'params'.", call.=FALSE)
    params <- dots$which
    dots$which <- NULL
  }
  mainTitle <- dots$main # this goes in outer margin
  dots$main <- NULL
  if(is.null(mainTitle))
    mainTitle <- paste("Diagnostics for", deparse(substitute(x)))
  defaultArgs <- list(xlab="", ylab="", type='l', lty=1)
  useArgsT <- useArgsD <- modifyList(defaultArgs, dots)
  selPlot <- names(useArgsT) %in%
    c(names(as.list(args(title))), names(par(no.readonly=TRUE)))
  titleArgs <- useArgsT[selPlot]

  # Deal with numerical input
  arrList <- switch(class(x)[1],
      jagsUI  = getArrA(x$samples, x$summary, x$mcmc.info$n.chains),
      bugs    = getArrA(x$sims.array, x$summary, x$n.chains),
      rjags   = getArrA(x$BUGSoutput$sims.array, x$BUGSoutput$summary, x$BUGSoutput$n.chains),
      mcmc.list = getArrB(x, length(x)),
      runjags = getArrB(x$mcmc, length(x$mcmc)),
      getArrB(x))

  mcmc3d <- arrList$mcmc3d
  if(!is.numeric(mcmc3d) || length(dim(mcmc3d)) != 3)
    stop("Could not extract numeric MCMC chains from your input.", call.=FALSE)
  Rhat <- round(arrList$stats[, 'Rhat'], 3)
  # n.eff <- round(arrList$stats[, 'n.eff'])
  MCEpc <- round(arrList$stats[, 'MCEpc'], 2)

  niter <- dim(mcmc3d)[1]
  nchains <- dim(mcmc3d)[2]
  npars <- dim(mcmc3d)[3]
  parnames <- dimnames(mcmc3d)[[3]]
  if(is.null(parnames))
    parnames <- paste0("V", 1:npars)

  # Deal with subsetting
  if(!missing(params)) {
    params <- matchStart(params, parnames)
    if(length(params) == 0)
      stop("No columns match the specification in 'params'.", call.=FALSE)
    mcmc3d <- mcmc3d[, , params, drop=FALSE]
    Rhat <- Rhat[params]
    # n.eff <- n.eff[params]
    MCEpc <- MCEpc[params]
    parnames <- dimnames(mcmc3d)[[3]]
    npars <- dim(mcmc3d)[3]
  }
  if(!missing(howMany) && abs(howMany) < niter) {
    if(howMany > 0) {
      mcmc3d <- mcmc3d[1:howMany, , , drop=FALSE]
    } else {
      mcmc3d <- mcmc3d[(niter+howMany+1):niter, , , drop=FALSE]
    }
  }
  if(!missing(chains) && all(chains <= nchains)) {
    mcmc3d <- mcmc3d[, chains, , drop=FALSE]
    nchains <- length(chains)
  }

  # Do the plots
  # ------------
  old.par <- par(mar = c(2,2,2,0)+0.1, oma=c(1,1,1,1), "mfrow")
    on.exit(par(old.par))
  if(!is.null(mainTitle))
    par(oma=c(1,1,3,1))

  if(npars > maxRow) {
    old.ask <- devAskNewPage(ask)
    on.exit(devAskNewPage(old.ask), add=TRUE)
  }
  nrows <- min(npars, maxRow)
  layout(matrix(1:(nrows*2), ncol=2, byrow=TRUE), widths=2:1)
  for(i in 1:npars) {
    redFlag <- !is.na(Rhat[i]) && Rhat[i] > RhatBad
    mat <- matrix(mcmc3d[, , i], ncol=nchains) # need 'matrix' if 1 chain
    if(any(is.na(mat))) {
      warning("The chain '", parnames[i], "' contains NAs and cannot be plotted.", call.=FALSE)
      next
    }
    # do trace plot
    useArgsT$ylab <- parnames[i]
    useArgsT$y <- mat
    do.call(matplot, useArgsT)
    abline(h=mean(mat))
    titleArgs$main <- paste0(parnames[i], ": Rhat = ", Rhat[i])
    titleArgs$line <- 0.3
    titleArgs$adj <- 1
    titleArgs$col.main <- 1 + redFlag
    titleArgs$outer <- FALSE
    do.call(title, titleArgs)
    if(redFlag)
      box(col=2, lwd=2)
    # do density plot
    density0(mat, useArgsD)
    # titleArgs$main <- paste0("n.eff = ", n.eff[i])
    titleArgs$main <- paste0("MCE% = ", MCEpc[i])
    titleArgs$adj <- 0
    do.call(title, titleArgs)
    if(redFlag)
      box(col=2, lwd=2)
    if(!is.null(mainTitle)) {
      titleArgs$main <- mainTitle
      titleArgs$line <- dots$line
      titleArgs$adj <- dots$adj
      titleArgs$col.main <- dots$col.main
      titleArgs$outer <- TRUE
      do.call(title, titleArgs)
    }
  }
}
# ..........................................................

# Helper function to do density plot (or lollipops) for 1 parameter
# mat : matrix of MCMC output with 1 column per chain
# plotArgs : list with plotting parameters
density0 <- function(mat, plotArgs, ...)  {

  bw <- bw.nrd0(mat)
  # lollipops or density plot?
  # unik <- unique.default(mat)
  # if(length(unik) == 1) {
  if(bw == 0 || diff(range(mat)) < sqrt(.Machine$double.eps)) {
    plot(1, 1, type = "n", ann = FALSE, axes = FALSE)
    text(1,1, paste("All values are the same:\n", signif(mat[1], 4)))
  } else if(all(mat %% 1 == 0) && all(mat >= 0) && diff(range(mat)) < 50) {
    # "lollipops"
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

