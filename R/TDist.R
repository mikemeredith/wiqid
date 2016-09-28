
# The functions 'dt', 'pt', etc in R calculate the (cumulative) probability
#   density for a t-distribution given the t-statistic and the number of
#   degrees of freedom (df). That means you have to calculate the t-stat first,
#   unlike 'dnorm'/'pnorm', where you put in the mean and sd.

# The function 'dt2' and 'pt2' calculate these values with given mean, se,
#   and df, just like 'dnorm'/'dnorm' with the addition of 'df'.

dt2 <- function(x, mean, sd, df) {
  tstat <- (x - mean)/sd
  return(dt(tstat, df))
}

pt2 <- function(x, mean, sd, df, lower.tail=TRUE, log.p=FALSE) {
  tstat <- (x - mean)/sd
  return(pt(tstat, df, lower.tail=lower.tail, log.p=log.p))
}

qt2 <- function(p, mean, sd, df, lower.tail=TRUE, log.p=FALSE) {
  tstat <- qt(p, df, lower.tail=lower.tail, log.p=log.p)
  return(tstat * sd + mean)
}

rt2 <- function(n, mean, sd, df) {
  tstat <- rt(n, df)
  return(tstat * sd + mean)
}

