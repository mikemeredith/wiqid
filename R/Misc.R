# Function to calculate the MARK-style confidence intervals for N
# See help for Closed Captures
# Not exported

getMARKci <- function(beta, SE.beta, ci) {
  f0.hat <- exp(beta)
  crit <- qnorm((1 - ci[1]) / 2, lower.tail=FALSE)
  C <- exp(crit * sqrt(log(1 + SE.beta^2))) # See the Burnham et al reference!
  return(c(f0.hat, f0.hat/C, f0.hat*C))
}

# Creates a table from a vector of AIC-type criterion values
# Exported

makeAICtable <- function(x) {
  xlab = deparse(substitute(x))
  delta <- x - min(x, na.rm=TRUE)
  ModelLik <- exp( - delta / 2)
  ModelWt <- ModelLik / sum(ModelLik, na.rm=TRUE)
  out <- cbind(x, delta, ModelLik, ModelWt)
  colnames(out)[1] <- xlab
  if(!is.null(rownames(out))) { # only sort if rows are named
    ord <- order(x)
    out <- out[ord, ]
  }
  return(out)
}