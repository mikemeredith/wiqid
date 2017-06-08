
# Function to estimate parameters of Robust Design model
#  with "ad hoc" strategy:
# N and p are estimated beforehand with a closedCap* function


survRDah <- function(DH, freq=1, occsPerSeason, N, pStar)  {

  # Do sanity checks here
  stopifnot(length(occsPerSeason) == 1)  # For the moment!

  K <- ncol(DH) / occsPerSeason # Number of seasons
  stopifnot(length(N) == K)
  stopifnot(length(pStar) == K)
  seasonID <- rep(1:K, each=occsPerSeason)
  stopifnot(length(seasonID) == ncol(DH))

  # turn the DH into a season-wise DH and do m-array
  getDHseason <- function(dh)
    tapply(dh, as.factor(seasonID), max)
  DHseason <- t(apply(DH, 1, getDHseason))
  mMat <- ch2mArray(DHseason, freq=freq)

  param <- rep(0, K-1)
  # Log likelihood function
  nll <- function(param)  {
    log_phi <- plogis(param, log.p=TRUE)
    nll <- -sum(mMat * log_qArray(log_phi, log(pStar[-1]), log(1 - pStar[-1])))
    return(min(nll, .Machine$double.xmax))
  }

  res <- nlm(nll, param)  ## , hessian=TRUE)
  phiHat <- plogis(res$estimate)

  bHat <- N[-1] / N[-K] - phiHat   # return this

  return(list(
    phiHat = phiHat,
    bHat = bHat,
    pStarHat = pStar,
    Nhat = N))
}
