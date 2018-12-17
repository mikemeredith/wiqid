
# test_that code for BoccSS0 and BoccSS

# library(testthat)
require(wiqid)

context("Bayesian occupancy models")

test_that("BoccSS0 gives same answers",  {
  # Use the salamanders data
  data(salamanders)
  y <- rowSums(salamanders)
  n <- rowSums(!is.na(salamanders))
  set.seed(123)
  Bout <- BoccSS0(y, n)
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(30000, 2)))
  expect_that(names(Bout), equals(c("psi", "p")))
  expect_that(attr(Bout, "header"),
    equals("Model fitted in R with a Gibbs sampler"))
  expect_that(attr(Bout, "n.chains"), equals(3))
  expect_equivalent(round(attr(Bout, "n.eff")), c(4440, 5999))
  expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1, 1))
  expect_equal(as.character(attr(Bout, "call")), c("BoccSS0", "y", "n"))

  expect_equivalent(round(colMeans(Bout), 4), c(0.6152, 0.2580))
  expect_equivalent(round(c(hdi(Bout)), 4),
    c(0.3770, 0.8620, 0.1504, 0.3660))

  set.seed(234)
  Bout <- BoccSS0(y, n, psiPrior=c(5,5), pPrior=c(3,3),
                    chains=4, sample=4000, burnin=10)
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(4000, 2)))
  expect_that(attr(Bout, "n.chains"), equals(4))
  expect_equivalent(round(attr(Bout, "n.eff")), c(1419, 1654))
  expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1.000, 1.001))
  expect_equal(as.character(attr(Bout, "call")),
    c("BoccSS0", "y", "n", "c(5, 5)", "c(3, 3)", "4", "4000", "10" ))

  expect_equivalent(round(colMeans(Bout), 4), c(0.5622, 0.2814))
  expect_equivalent(round(c(hdi(Bout)), 4),
    c(0.3998, 0.7544, 0.1771, 0.3789))
})
# ............................................................

if(parallel::detectCores() > 3) {
  test_that("BoccSS parallel gives same answers",  {
    # Use the weta data
    data(weta)
    DH <- weta[, 1:5]
    weta.covs <- weta[, 6:11]

    expect_message({Bout <- BoccSS(DH, sample=3000, burnin=100, seed=123)},
      "Starting MCMC run for 3 chains with 1100 iterations")
    expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
    expect_that(dim(Bout), equals(c(3000, 2)))
    expect_that(names(Bout), equals(c("psi_(Intercept)", "p_(Intercept)")))
    expect_that(attr(Bout, "header"),
      equals("Model fitted in R with a Gibbs sampler"))
    expect_that(attr(Bout, "n.chains"), equals(3))
    expect_equal(as.character(attr(Bout, "call")),
      c("BoccSS", "DH", "3000", "100", "123"))
    expect_equivalent(round(attr(Bout, "n.eff")), c(256, 351))
    expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1.020, 1.010))
    expect_equivalent(round(colMeans(Bout), 4), c(0.3479, -0.3970))
    expect_equivalent(round(c(hdi(Bout)), 4),
      c(-0.1277, 0.8399, -0.6700, -0.1181))
    expect_message({Bout <- BoccSS(DH, model=list(psi~Browsed-1, p~.Time),
      data=weta,
      priors=list(sigmaPsi=c(1,1)), chains=2, sample=2000, burnin=100,
      seed=234)}, " ")
    expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
    expect_that(dim(Bout), equals(c(2000, 4)))
    expect_that(attr(Bout, "n.chains"), equals(2))
    expect_equivalent(round(attr(Bout, "n.eff")), c(332, 159, 355, 943))
    expect_equivalent(round(attr(Bout, "Rhat"), 3),
      c(1.004, 1.051, 1.014, 1.001))
    expect_equal(as.character(attr(Bout, "call")),
      c( "BoccSS", "DH", "list(psi ~ Browsed - 1, p ~ .Time)",
      "weta", "list(sigmaPsi = c(1, 1))", "2", "2000", "100", "234"))

    expect_equivalent(round(colMeans(Bout), 4),
      c(0.0150, 0.7826, -0.4222, 0.3359))
    expect_equivalent(round(c(hdi(Bout)), 3),
      c(-0.515, 0.583, 0.053, 1.734, -0.679, -0.144, -0.034, 0.756))
  })
}
# ........................................................

test_that("BoccSS sequential gives same answers",  {
  # Use the weta data
  data(weta)
  DH <- weta[, 1:5]
  weta.covs <- weta[, 6:11]

  expect_message({Bout <- BoccSS(DH, sample=3000, burnin=100, seed=123, parallel=FALSE)},
  "Starting MCMC run for 3 chains with 1100 iterations")
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(3000, 2)))
  expect_that(names(Bout), equals(c("psi_(Intercept)", "p_(Intercept)")))
  expect_that(attr(Bout, "header"), equals("Model fitted in R with a Gibbs sampler"))
  expect_that(attr(Bout, "n.chains"), equals(3))
  expect_equal(as.character(attr(Bout, "call")),
    c("BoccSS", "DH", "3000", "100", "FALSE", "123"))
  expect_equivalent(round(attr(Bout, "n.eff")), c(193, 313))
  expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1.047, 1.002))
  expect_equivalent(round(colMeans(Bout), 4), c(0.3593, -0.4058))
  expect_equivalent(round(c(hdi(Bout)), 4),
    c(-0.1667,  0.8787, -0.7005, -0.1146))
  expect_message({Bout <- BoccSS(DH, model=list(psi~Browsed-1, p~.Time), data=weta,
    priors=list(sigmaPsi=c(1,1)), chains=1, sample=1000, burnin=100,
    seed=234)}, " ")
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(1000, 4)))
  expect_that(attr(Bout, "n.chains"), equals(1))
  expect_equivalent(round(attr(Bout, "n.eff")), c(172, 77, 112, 460))
  expect_true(is.null(attr(Bout, "Rhat")))
  expect_equal(as.character(attr(Bout, "call")),
    c( "BoccSS", "DH", "list(psi ~ Browsed - 1, p ~ .Time)",
    "weta", "list(sigmaPsi = c(1, 1))", "1", "1000", "100", "234"))

  expect_equivalent(round(colMeans(Bout), 4),
    c(0.0316, 0.7624, -0.4214, 0.3181))
  expect_equivalent(round(c(hdi(Bout)), 4),
    c(-0.5437, 0.5730, -0.0126, 1.5376, -0.6816, -0.1360, -0.0682, 0.7056))
})
