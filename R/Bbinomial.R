
# Simulated draws from posterior of a binomial likelihood with beta prior
# =======================================================================

# The beta prior is specified by mode and concentration.

Bbinom <- function(y, n, priors=NULL, sample=50000) {

  if(!is.null(priors$conc) && priors$conc < 2)
    stop("priors$conc must not be less than 2.")
  if(!is.null(priors$mode) && (priors$mode < 0 || priors$mode > 1 ))
    stop("priors$mode must be between 0 and 1.")
  if(y > n)
    stop("Number of successes (y) cannot be greater than the number of trials (n).")

  if(!is.null(priors$conc) && !is.null(priors$mode)) {
    pr1 <- priors$mode * (priors$conc - 2) + 1
    pr2 <- (1 - priors$mode) * (priors$conc - 2) + 1
  } else {
    pr1 <- pr2 <- 1
  }

  po1 <- pr1 + y
  po2 <- pr2 + n - y

  post <- rbeta(sample, po1, po2)

  out <- as.Bwiqid(data.frame(pi = post),
      header = "Sample drawn from beta posterior distribution",
      defaultPlot = "pi")
  attr(out, "call") <- match.call()
  attr(out, "n.chains") <- 1
  return(out)
}
