
# test_that code for Bnormal

# library(testthat)
require(wiqid)

context("Bayesian normal models")

test_that("Bnormal gives same answers",  {
  # Generate data
  set.seed(123)
  x <- rnorm(10, 1, 0.15)
  Bout <- Bnormal(x)
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(30000, 2)))
  expect_that(names(Bout), equals(c("mu", "sigma")))
  expect_that(attr(Bout, "header"), equals("Model fitted in R with a Gibbs sampler"))
  expect_that(attr(Bout, "n.chains"), equals(3))
  expect_equivalent(round(attr(Bout, "n.eff")), c(29897, 23590))
  expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1, 1))
  expect_equal(as.character(attr(Bout, "call")), c("Bnormal", "x"))

  expect_equivalent(round(colMeans(Bout), 4), c(1.0112, 0.1572))
  expect_equivalent(round(c(hdi(Bout)), 4), c(0.9083, 1.1122, 0.0906, 0.2419))
  xx <- x * 1000
  expect_warning(Bout <- Bnormal(xx),
      "Sample mean is outside the prior range")
  set.seed(345)
  Bout <- Bnormal(xx, priors=NULL)  # flat priors
  expect_equivalent(round(colMeans(Bout), 1), c(1010.5, 156.3))
  expect_equivalent(round(c(hdi(Bout)), 1), c(912.5, 1114.7, 90.0, 241.0))
  set.seed(456)
  Bout <- Bnormal(xx, priors=list(mMean=1000))  # sensible priors
  expect_equivalent(round(colMeans(Bout), 1), c(1009.0, 154.2))
  expect_equivalent(round(c(hdi(Bout)), 1), c(921.8, 1097.2, 90.0, 234.4))

  set.seed(123)
  x <- rnorm(5, 1, 0.15)
  Bout <- Bnormal(x, n.iter=500, n.burnin=50, n.chains=4)
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(1800, 2)))
  expect_equivalent(round(colMeans(Bout), 4), c(1.0263, 0.1601))
  expect_equivalent(round(c(hdi(Bout)), 4), c(0.8606, 1.1920, 0.0648, 0.3159))
})