# Print and plot methods for class Bwiqid, ie. MCMC output

print.Bwiqid <- function(x, digits=4, ...)  {
  if(!inherits(x, "data.frame"))
    stop("x is not a valid Bwiqid object")
  header <- attr(x, "header")
  MCerror <- attr(x, "MCerror")
  Rhat <- attr(x, "Rhat")
  n.eff <- attr(x, "n.eff")
  timetaken <- attr(x, "timetaken")

  toPrint <- cbind(
    mean = colMeans(x),
    sd = apply(x, 2, sd),
    median = apply(x, 2, median), 
    t(hdi(x)))
  colnames(toPrint)[4:5] <- c("HDIlo", "HDIup")
  if(!is.null(MCerror))
    toPrint <- cbind(toPrint, 'MCE%' = round(100 * MCerror/toPrint[, 'sd'], 1))
  if(!is.null(Rhat))
    toPrint <- cbind(toPrint, Rhat = Rhat)
  if(!is.null(n.eff))
    toPrint <- cbind(toPrint, n.eff = round(n.eff))

  toPrint0 <- unique(toPrint)
    
  if(is.null(header))
    header <- "MCMC fit results:"
  cat(header, "\n")
  cat(nrow(x), "simulations saved.\n")
  if(nrow(toPrint0) < nrow(toPrint))
    cat("(Duplicate rows removed.)\n")
  print(toPrint0, digits = digits)
  cat("\n'HDIlo' and 'HDIup' are the limits of a 95% HDI credible interval.\n")
  if(!is.null(MCerror))
    cat("'MCE%' is the Monte Carlo error as a %age of the SD (should be less than 5%).\n")
  if(!is.null(Rhat))
    cat("'Rhat' is the potential scale reduction factor (at convergence, Rhat=1).\n")
  if(!is.null(n.eff))
    cat("'n.eff' is a crude measure of effective sample size.\n")
  if(!is.null(timetaken)) {
    took <- format(round(timetaken, 1))
    cat("MCMC sample generation:", took, "\n")
  }
}
# .........................................................

plot.Bwiqid <-
function(x, which=NULL, credMass=0.95,
          ROPE=NULL, compVal=NULL, showCurve=FALSE,  showMode=FALSE,
          shadeHDI=NULL, ...) {
  # This function plots the posterior distribution for one selected item. 
  # Description of arguments:
  # x is mcmc.list object of the type returned by B* functions in 'wiqid'.
  # which indicates which item should be displayed; if NULL, looks for a 'toPlot' attribute in x; if missing does first column.
  # ROPE is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE.
  # compVal is a scalar specifying the value for comparison.
  # showCurve if TRUE the posterior should be displayed as a fitted density curve
  #   instead of a histogram (default).

  # TODO additional sanity checks.
  # Sanity checks:
  if(!inherits(x, "data.frame"))
    stop("x is not a valid Bwiqid object")
    
  # Deal with ... argument
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]

  if(is.null(which)) # && !is.null(attr(x, "defaultPlot")))
      which <- attr(x, "defaultPlot")
  if(is.null(which))
    which <- colnames(x)[1]
  if(is.na(match(which, colnames(x))))
    stop(paste("Could not find", which, "in the output"))  
  if(is.null(dots$xlab))
    dots$xlab <- which
  # Plot posterior distribution of selected item:
  out <- plotPost(x[[which]], credMass=credMass, ROPE=ROPE, compVal=compVal,
                  showCurve=showCurve, showMode=showMode, shadeHDI=shadeHDI,
                  graphicPars=dots)

  return(invisible(out))
}

