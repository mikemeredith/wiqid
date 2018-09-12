
# Wrappers for dbeta, pbeta, etc which use mean and sd as parameters.

# See https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance


getBetaPar <- function(mean, sd) {
  if(mean < 0 || mean > 1)
    stop("'mean' must be between 0 and 1.")
  if(sd < 0 || sd > 0.5)
    stop("'sd' must be between 0 and 0.5.")
  alpha <- ((1 - mean) / sd^2 - 1 / mean) * mean ^ 2
  beta <- alpha * (1 / mean - 1)
  c(alpha, beta)
}

dbeta2 <- function(x, mean, sd) {
  shapes <- getBetaPar(mean,sd)
  return(dbeta(x, shapes[1], shapes[2]))
}

pbeta2 <- function(q, mean, sd, lower.tail=TRUE, log.p=FALSE) {
  shapes <- getBetaPar(mean,sd)
  return(pbeta(q, shapes[1], shapes[2], lower.tail=lower.tail, log.p=log.p))
}
qbeta2 <- function(p, mean, sd, lower.tail=TRUE, log.p=FALSE) {
  shapes <- getBetaPar(mean,sd)
  return(qbeta(p, shapes[1], shapes[2], lower.tail=lower.tail, log.p=log.p))
}

rbeta2 <- function(n, mean, sd) {
  shapes <- getBetaPar(mean,sd)
  return(rbeta(n, shapes[1], shapes[2]) * sd + mean)
}

