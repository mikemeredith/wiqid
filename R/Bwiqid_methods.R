# Print method for class Bwiqid, ie. MCMC output

print.Bwiqid <- function(x, digits=4, ...)  {
  if(!inherits(x, "data.frame"))
    stop("x is not a valid Bwiqid object")
  header <- attr(x, "header")
  MCerror <- attr(x, "MCerror")
  Rhat <- attr(x, "Rhat")
  n.eff <- attr(x, "n.eff")
  timing <- attr(x, "timing")

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
  if(!is.null(timing)) {
    took <- format(round(diff(timing), 1))
    cat("Run completed at ", format(timing[2]), ", took ", took, ".\n", sep="")
  }
}
# .........................................................

plot.Bwiqid <-
function(x, which=NULL, credMass=0.95,
                    ROPE=NULL, compVal=NULL, showCurve=FALSE, ...) {
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

  if(is.null(which) && !is.null(attr(x, "defaultPlot")))
      which <- attr(x, "defaultPlot")
  if(is.na(match(which, colnames(x)))) {
    warning(paste("Could not find", which, "in the output"))
    which <- colnames(x)[1]
  }   

  # Plot posterior distribution of selected item:
  plotPost(x[[which]], col="skyblue", credMass=credMass, ROPE=ROPE, showCurve=showCurve,
                  xlab=which , cex.lab = 1.75 , showMode=FALSE,
                  compVal=compVal, main=which, ...) 

  return(invisible(NULL))
}

