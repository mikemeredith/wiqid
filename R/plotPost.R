# Function taken from package BEST, original code by John Kruschke.
# Modified by Mike to make best use of ... argument and various other
#  enhancements, including shading of the HDI in showCurve=TRUE plots.

plotPost <-
function( paramDraws, credMass=0.95, compVal=NULL, ROPE=NULL,
           HDItextPlace=0.7, showMode=FALSE, showCurve=FALSE,
           shadeHDI=NULL, ... ) {

  # Does a plot for a single parameter. Called by plot.Bwiqid but also exported.
  # Returns a histogram object invisibly.
  # This stuff should be in the ... argument:
  #   yaxt="n", ylab="", xlab="Parameter", main="", cex.lab=1.5, cex=1.4,
  #   xlim=range(compVal, paramDraws), col="skyblue", border="white",
  #   breaks=NULL

  # Deal with ... argument:
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]

  if(!is.null(dots$paramSampleVec)) {
    message("*The 'paramSampleVec' argument is deprecated, please use 'paramDraws'.*")
    paramDraws <- dots$paramSampleVec
  }
  if(!is.numeric(paramDraws))
    stop("The first argument must be a vector of numbers.")

  defaultArgs <- list(xlab=deparse(substitute(paramDraws)),
    yaxt="n", ylab="", main="", cex.lab=1.5,
    cex=1.4, col="skyblue", border="white", bty="n", lwd=5, freq=FALSE,
    xlim=range(compVal, hdi(paramDraws, 0.99)))
  useArgs <- modifyList(defaultArgs, dots)
  # Get breaks argument
  breaks <- dots$breaks
  if (is.null(breaks)) {
    if (all(paramDraws == round(paramDraws))) { # all integers
      breaks <- seq(min(paramDraws), max(paramDraws) + 1) - 0.5
    } else {
      nbreaks <- ceiling(diff(range(paramDraws)) /
                          diff(hdi(paramDraws)) * 18)
      breaks <- seq( from=min(paramDraws), to=max(paramDraws),
                     length.out=nbreaks)
    }
  }
  histinfo <- hist(paramDraws, breaks=breaks, plot=FALSE)
  histinfo$xname <- useArgs$xlab

  if (showCurve) {
    densCurve <- densityFolded( paramDraws, adjust=2 )
    cenTendHt <- 0.9 * max(densCurve$y)  # For plotting
    selPlot <- names(useArgs) %in%
      c(names(as.list(args(plot.default))), names(par(no.readonly=TRUE)))
    plotArgs <- useArgs[selPlot]
    plotArgs$x <- densCurve$x
    plotArgs$y <- densCurve$y
    plotArgs$type <- "l"
    do.call(plot, plotArgs)
    abline(h=0, col='grey')
    # Display the HDI.
    if(!is.null(credMass)) {
      HDI <- hdi(densCurve, credMass, allowSplit=TRUE)
      ht <- attr(HDI, "height")
      if(nrow(HDI) == 1)  # hdi is not split
        HDI <- matrix(hdi(paramDraws, credMass), nrow=1)
      if(!is.null(shadeHDI))  {
        for (i in 1:nrow(HDI)) {
          inHDI <- which(densCurve$x >= HDI[i, 1] & densCurve$x <= HDI[i, 2])
          polyx <- c(HDI[i, 1], HDI[i, 1], densCurve$x[inHDI], HDI[i, 2], HDI[i, 2])
          polyy <- c(0, ht, densCurve$y[inHDI], ht, 0)
          polygon(polyx, polyy, border=NA, col=shadeHDI)
        }
      } else {
        segments(HDI, 0, HDI, ht, lty=2)
      }
      do.call(lines, plotArgs)
      segments(HDI[, 1], ht, HDI[, 2], ht, lwd=4, lend='butt')
      text( mean(HDI), ht, bquote(.(100*credMass) * "% HDI" ),
            adj=c(.5,-1.7), cex=useArgs$cex, xpd=TRUE )
      # text( HDI, ht, bquote(.(signif(HDI, 3))),
      text( HDI, ht, signifish(HDI, 3),
            pos=3, cex=useArgs$cex, xpd=TRUE )
    }
  } else {
    cenTendHt <- 0.9 * max(histinfo$density)  # For plotting
    plot.histogram.args.names <- c("freq", "density", "angle", "border",
      "main", "sub", "xlab", "ylab", "xlim", "ylim", "axes", "labels",
      "add") # plot.histogram not exported, so need to cheat!
    selPlot <- names(useArgs) %in%
      c(plot.histogram.args.names, names(par(no.readonly=TRUE)))
    plotArgs <- useArgs[selPlot]
    plotArgs$lwd <- 1
    plotArgs$x <- histinfo
    do.call(plot, plotArgs)
    # Display the HDI.
    if(!is.null(credMass)) {
      HDI <- hdi( paramDraws, credMass )
      lines(HDI, c(0,0), lwd=4, lend='butt')
      text( mean(HDI), 0, bquote(.(100*credMass) * "% HDI" ),
            adj=c(.5,-1.7), cex=useArgs$cex, xpd=TRUE )
      text( HDI[1], 0, signifish(HDI[1],3),
            adj=c(HDItextPlace,-0.5), cex=useArgs$cex, xpd=TRUE )
      text( HDI[2], 0, signifish(HDI[2],3),
            adj=c(1.0-HDItextPlace,-0.5), cex=useArgs$cex, xpd=TRUE )
    }
  }

  # Display mean or mode:
  if ( showMode==FALSE ) {
      meanParam <- mean( paramDraws )
      text( meanParam, cenTendHt,
            bquote(mean==.(signifish(meanParam,3))), adj=c(.5,0), cex=useArgs$cex, xpd=TRUE )
  } else {
      dres <- density( paramDraws )
      modeParam <- dres$x[which.max(dres$y)]
      text( modeParam, cenTendHt,
            bquote(mode==.(signifish(modeParam,3))), adj=c(.5,0), cex=useArgs$cex, xpd=TRUE )
  }
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    cvHt <- 0.7 * max(histinfo$density)
    cvCol <- "darkgreen"
    pcgtCompVal <- round( 100 * sum( paramDraws > compVal )
                          / length( paramDraws ) , 1 )
     pcltCompVal <- 100 - pcgtCompVal
     lines( c(compVal,compVal), c(0.96*cvHt,0),
            lty="dotted", lwd=1, col=cvCol )
     text( compVal, cvHt,
           bquote( .(pcltCompVal)*"% < " *
                   .(signifish(compVal,3)) * " < "*.(pcgtCompVal)*"%" ),
           adj=c(pcltCompVal/100,0), cex=0.8*useArgs$cex, col=cvCol, xpd=TRUE )
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    ROPEtextHt <- 0.55 * max(histinfo$density)
    ropeCol <- "darkred"
     pcInROPE <- ( sum( paramDraws > ROPE[1] & paramDraws < ROPE[2] )
                          / length( paramDraws ) )
     lines( c(ROPE[1],ROPE[1]), c(0.96*ROPEtextHt,0), lty="dotted", lwd=2,
            col=ropeCol )
     lines( c(ROPE[2],ROPE[2]), c(0.96*ROPEtextHt,0), lty="dotted", lwd=2,
            col=ropeCol)
     text( mean(ROPE), ROPEtextHt,
           bquote( .(round(100*pcInROPE))*"% in ROPE" ),
           adj=c(.5,0), cex=1, col=ropeCol, xpd=TRUE )
  }

  return(invisible(histinfo))
}
