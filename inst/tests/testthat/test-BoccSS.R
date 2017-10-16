
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
  expect_that(attr(Bout, "header"), equals("Model fitted in R with a Gibbs sampler"))
  expect_that(attr(Bout, "n.chains"), equals(3))
  expect_equivalent(round(attr(Bout, "n.eff")), c(4440, 5999))
  expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1, 1))
  expect_equal(as.character(attr(Bout, "call")), c("BoccSS0", "y", "n"))

  expect_equivalent(round(colMeans(Bout), 4), c(0.6152, 0.2580))
  expect_equivalent(round(c(hdi(Bout)), 4), c(0.3770, 0.8620, 0.1504, 0.3660))

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
  expect_equivalent(round(c(hdi(Bout)), 4), c(0.3998, 0.7544, 0.1771, 0.3789))
})
# ............................................................

if(parallel::detectCores() > 3) {
  test_that("BoccSS parallel gives same answers",  {
    # Use the weta data
    data(weta)
    DH <- weta[, 1:5]
    weta.covs <- weta[, 6:11]

    expect_output({Bout <- BoccSS(DH, sample=3000, burnin=100, seed=123)},
    "Starting MCMC run for 3 chains with 1100 iterations.\n\nMCMC run complete.")
    expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
    expect_that(dim(Bout), equals(c(3000, 2)))
    expect_that(names(Bout), equals(c("psi_(Intercept)", "p_(Intercept)")))
    expect_that(attr(Bout, "header"), equals("Model fitted in R with a Gibbs sampler"))
    expect_that(attr(Bout, "n.chains"), equals(3))
    expect_equal(as.character(attr(Bout, "call")), c("BoccSS", "DH", "3000", "100", "123"))
    if(getRversion() < "3.3.0") {
      expect_equivalent(round(attr(Bout, "n.eff")), c(299, 399))
      expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1.009, 1.010))
      expect_equivalent(round(colMeans(Bout), 4), c(0.3426, -0.4035))
      expect_equivalent(round(c(hdi(Bout)), 4),
        c(-0.1232, 0.8500, -0.6895, -0.1363))
    } else {
      expect_equivalent(round(attr(Bout, "n.eff")), c(248, 361))
      expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1.003, 1.002))
      expect_equivalent(round(colMeans(Bout), 4), c(0.3927, -0.4286))
      expect_equivalent(round(c(hdi(Bout)), 4),
        c(-0.1096, 0.9313, -0.7059, -0.1553))
    }
    expect_output({Bout <- BoccSS(DH, model=list(psi~Browsed-1, p~.Time), data=weta,
      priors=list(sigmaPsi=c(1,1)), chains=2, sample=2000, burnin=100,
      seed=234)}, " ")
    expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
    expect_that(dim(Bout), equals(c(2000, 4)))
    expect_that(attr(Bout, "n.chains"), equals(2))
    expect_equivalent(round(attr(Bout, "n.eff")), c(398, 175, 228, 844))
    expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1.007, 1.006, 1.004, 1.001))
    expect_equal(as.character(attr(Bout, "call")),
      c( "BoccSS", "DH", "list(psi ~ Browsed - 1, p ~ .Time)",
      "weta", "list(sigmaPsi = c(1, 1))", "2", "2000", "100", "234"))

    expect_equivalent(round(colMeans(Bout), 4), c(0.0207, 0.7862, -0.4192, 0.3369))
    expect_equivalent(round(c(hdi(Bout)), 3),
      c(-0.497, 0.563, 0.059, 1.589, -0.689, -0.160, -0.053, 0.760))
  })
}
# ........................................................

test_that("BoccSS sequential gives same answers",  {
  # Use the weta data
  data(weta)
  DH <- weta[, 1:5]
  weta.covs <- weta[, 6:11]

  expect_output({Bout <- BoccSS(DH, sample=3000, burnin=100, seed=123, parallel=FALSE)},
  "Starting MCMC run for 3 chains with 1100 iterations.\n\nMCMC run complete.")
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(3000, 2)))
  expect_that(names(Bout), equals(c("psi_(Intercept)", "p_(Intercept)")))
  expect_that(attr(Bout, "header"), equals("Model fitted in R with a Gibbs sampler"))
  expect_that(attr(Bout, "n.chains"), equals(3))
  expect_equal(as.character(attr(Bout, "call")), c("BoccSS", "DH", "3000", "100", "FALSE", "123"))
  if(getRversion() < "3.3.0") {
    expect_equivalent(round(attr(Bout, "n.eff")), c(155, 224))
    expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1.642, 1.228))
    expect_equivalent(round(colMeans(Bout), 4), c(1.0043, -0.4708))
    expect_equivalent(round(c(hdi(Bout)), 4),
      c(-0.1668, 5.2056, -0.8808, -0.1557))
  } else {
    expect_equivalent(round(attr(Bout, "n.eff")), c(174, 288))
    expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1.014, 1.003))
    expect_equivalent(round(colMeans(Bout), 4), c(0.3788, -0.4145))
    expect_equivalent(round(c(hdi(Bout)), 4),
      c(-0.1052,  1.0259, -0.7200, -0.1142))
  }
  expect_output({Bout <- BoccSS(DH, model=list(psi~Browsed-1, p~.Time), data=weta,
    priors=list(sigmaPsi=c(1,1)), chains=1, sample=1000, burnin=100,
    seed=234)}, " ")
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(1000, 4)))
  expect_that(attr(Bout, "n.chains"), equals(1))
  expect_equivalent(round(attr(Bout, "n.eff")), c(173, 49, 96, 341))
  expect_true(is.null(attr(Bout, "Rhat")))
  expect_equal(as.character(attr(Bout, "call")),
    c( "BoccSS", "DH", "list(psi ~ Browsed - 1, p ~ .Time)",
    "weta", "list(sigmaPsi = c(1, 1))", "1", "1000", "100", "234"))

  expect_equivalent(round(colMeans(Bout), 4), c(0.0569, 0.9520, -0.4501, 0.3449))
  expect_equivalent(round(c(hdi(Bout)), 4),
    c(-0.4745, 0.7187, 0.0238, 2.1089, -0.7293, -0.1732, -0.0374, 0.7599))
})
