# Function to calculate the MARK-style confidence intervals for N
# See help for Closed Captures
# Not exported

getMARKci <- function(beta, SE.beta, ci) {
  f0.hat <- exp(beta)
  crit <- qnorm((1 - ci[1]) / 2, lower.tail=FALSE)
  C <- exp(crit * sqrt(log(1 + SE.beta^2))) # See the Burnham et al reference!
  return(c(f0.hat, f0.hat/C, f0.hat*C))
}

# Creates a table from output of AIC or AICc
# Exported

AICtable <- function(x) {
  if(is.vector(x)) {
    name <- deparse(substitute(x))
    IC <- x
    x <- data.frame(x)
    colnames(x) <- name
  } else if (is.data.frame(x) && ncol(x) > 1 )  {
    IC <- x[, 2] 
  } else {
    stop("x must be a vector or a data frame with > 1 column")
  }
  Delta <- IC - min(IC, na.rm=TRUE)
  ModelLik <- exp( - Delta / 2)
  ModelWt <- ModelLik / sum(ModelLik, na.rm=TRUE)
  out <- cbind(x, Delta, ModelLik, ModelWt)
  if(!is.null(rownames(out))) { # only sort if rows are named
    ord <- order(IC)
    out <- out[ord, ]
  }
  return(out)
}