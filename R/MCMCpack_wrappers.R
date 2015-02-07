

# This file has wrappers for certain functions in the MCMCpack package

# The objective is to have the same calling syntax as other wiqid functions
#  and to return objects of class Bwiqid.


Bregress <- function(formula, data = NULL, priors=list(),
              n.iter=11000, n.burnin=1000, ...)  {
  startTime <- Sys.time()

  priorsDefault <- list(b0=0, B0=0, c0=0.001, d0=0.001)
  priors <- replace(priorsDefault, names(priors), priors)

  mcmc <- MCMCregress(formula=formula, data=data,
            burnin = n.burnin, mcmc = n.iter-n.burnin,
            b0=priors$b0, B0=priors$B0, c0=priors$c0, d0=priors$d0, ...)

  if(any(abs(geweke.diag(mcmc)$z) > 1.96))
    warning("Geweke diagnostic > 1.96; chain may not have converged.")

  out <- as.Bwiqid(mcmc,
    header="Regression line fitted by MCMCpack::MCMCregress")
  attr(out, "call") <- match.call()
  attr(out, "timetaken") <- Sys.time() - startTime
  return(out)
}
# ...........................................................................

Blogit <- function(formula, data = NULL, priors=list(),
              n.iter=11000, n.burnin=1000, ...)  {
  startTime <- Sys.time()

  priorsDefault <- list(b0=0, B0=0, c0=0.001, d0=0.001)
  priors <- replace (priorsDefault, names(priors), priors)

  
  mcmc <- MCMClogit(formula=formula, data=data,
            burnin = n.burnin, mcmc = n.iter-n.burnin,
            b0=priors$b0, B0=priors$B0, c0=priors$c0, d0=priors$d0, ...)

  if(any(abs(geweke.diag(mcmc)$z) > 1.96))
    warning("Geweke diagnostic > 1.96; chain may not have converged.")

  out <- as.Bwiqid(mcmc,
    header="Logistic regression fitted by MCMCpack::MCMClogit")
  attr(out, "call") <- match.call()
  attr(out, "timetaken") <- Sys.time() - startTime
  return(out)
}
# ...........................................................................


