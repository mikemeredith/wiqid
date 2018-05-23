# Print, summary, plot, window, head and tail methods for class Bwiqid, ie. MCMC output


print.Bwiqid <- function(x, digits=4, ...)  {
  if(!inherits(x, "data.frame"))
    stop("x is not a valid Bwiqid object")
  call <- attr(x, "call")
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

  # if(!is.null(call))
    # cat("Call:", call, "\n")
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


summary.Bwiqid <- function(object, digits=3, ...)  {
  if(!inherits(object, "data.frame"))
    stop("object is not a valid Bwiqid object")
  call <- attr(object, "call")
  header <- attr(object, "header")
  MCerror <- attr(object, "MCerror")
  Rhat <- attr(object, "Rhat")
  n.eff <- attr(object, "n.eff")
  timetaken <- attr(object, "timetaken")

  toPrint <- cbind(
    mean = colMeans(object),
    sd = apply(object, 2, sd),
    median = apply(object, 2, median),
    t(hdi(object)))
  colnames(toPrint)[4:5] <- c("HDIlo", "HDIup")
  if(!is.null(MCerror))
    toPrint <- cbind(toPrint, 'MCE%' = round(100 * MCerror/toPrint[, 'sd'], 1))
  if(!is.null(Rhat))
    toPrint <- cbind(toPrint, Rhat = Rhat)
  if(!is.null(n.eff))
    toPrint <- cbind(toPrint, n.eff = round(n.eff))

  if(is.null(header))
    header <- "MCMC fit results:"
  cat(header, "\n")
  cat(nrow(object), "simulations saved.\n")
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
  return(invisible(round(toPrint, digits=digits)))
}
# .........................................................

plot.Bwiqid <-
function(x, which=NULL, credMass=0.95,
          ROPE=NULL, compVal=NULL, showCurve=FALSE,  showMode=FALSE,
          shadeHDI=NULL, ...) {
  # This function plots the posterior distribution for one selected item.
  # Description of arguments:
  # x is mcmc.list object of the type returned by B* functions in 'wiqid'.
  # which indicates which item should be displayed; if NULL, looks for a 'toPlot'
  #   attribute in x; if missing does first column.
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
  if(!is.character(which))
    stop("'which' must be an object of class 'character'.") # Added 2017-09-27
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

# .........................................................

window.Bwiqid <- function(x, start=NULL, end=NULL, thin=1, ...)  {
  if(!inherits(x, "Bwiqid"))
    stop("x is not a valid Bwiqid object")
  n.chains <- attr(x, "n.chains")
  if(is.null(n.chains) || nrow(x) %% n.chains != 0)
    stop("n.chains attribute of x is missing or invalid.")
  chain.length <- nrow(x) / n.chains
  if(is.null(start) || start > chain.length)
    start <- 1
  if(is.null(end) || end > chain.length || end <= start)
    end <- chain.length
  if(start < 1 || end < 1 || thin < 1)
    stop("Arguments start, end, and thin must be integers > 1")

  x_new <- vector('list', ncol(x))
  for( i in 1:ncol(x)) {
    mat <- matrix(x[, i], ncol=n.chains)
    mat_new <- mat[seq(start, end, by=thin), ]
    x_new[[i]] <- as.vector(mat_new)
  }
  names(x_new) <- colnames(x)
  x_df <- as.data.frame(x_new)
  rownames(x_df) <- NULL
  class(x_df) <- class(x)
  attr(x_df, "header") <- attr(x, "header")
  attr(x_df, "n.chains") <- attr(x, "n.chains")
  attr(x_df, "defaultPlot") <- attr(x, "defaultPlot")
  attr(x_df, "timetaken") <- attr(x, "timetaken")

  return(x_df)
}

# .........................................................

head.Bwiqid <- function(x, n=6L, ...) {
  head(as.data.frame(x), n=n, ...)
}

tail.Bwiqid <- function(x, n=6L, ...) {
  tail(as.data.frame(x), n=n, ...)
}


