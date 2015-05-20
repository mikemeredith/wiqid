
# test_that code for Bnormal

# library(testthat)
require(wiqid)

context("Bayesian normal models")

test_that("Bnormal gives same answers",  {
  # Generate data
  set.seed(123)
  x <- rnorm(10, 1, 0.15)
  Bout <- Bnormal(x, rnd.seed=345, verbose=FALSE)
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(100002, 2)))
  expect_that(names(Bout), equals(c("mu", "sigma")))
  expect_that(attr(Bout, "header"),
    equals("Model fitted in JAGS with jagsUI::jags.basic"))
  expect_that(attr(Bout, "n.chains"), equals(3))
  expect_equivalent(round(attr(Bout, "n.eff")), c(103943, 22727 ))
  expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1, 1.003))
  expect_equal(as.character(attr(Bout, "call")), c("Bnormal", "x", "FALSE", "345"))

  expect_equivalent(round(colMeans(Bout), 4), c(1.0114, 0.1685))
  expect_equivalent(round(c(hdi(Bout)), 4), c(0.9010, 1.1222, 0.0932, 0.2665))
  xx <- x * 1000
  expect_warning(Bout <- Bnormal(xx, priors=list(muM=0, muSD=10)), # silly prior
      "Sample mean is outside the prior range")
  Bout <- Bnormal(xx, priors=list(), rnd.seed=345, verbose=FALSE)  # minimally informative gamma priors
  expect_equivalent(round(colMeans(Bout), 1), c(1011.4, 167.1))
  expect_equivalent(round(c(hdi(Bout)), 1), c(900.4, 1119.4,   92.1,  261.3))
  Bout <- Bnormal(xx, priors=list(muM=1000), rnd.seed=345, verbose=FALSE)  # sensible priors
  expect_equivalent(round(colMeans(Bout), 1), c(1011.2, 167.2))
  expect_equivalent(round(c(hdi(Bout)), 1), c(900.4, 1119.0, 92.4, 261.4))

})