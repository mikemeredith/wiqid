
# test_that code for Bnormal

# library(testthat)
require(wiqid)

context("Bayesian normal models")

test_that("Bnormal gives same answers",  {
  # Generate data
  set.seed(123)
  x <- rnorm(10, 1, 0.15)
  # With default (flat) priors:
  Bout <- Bnormal(x)
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(50001, 2)))
  expect_that(names(Bout), equals(c("mu", "sigma")))
  expect_that(attr(Bout, "header"),
    equals("Model fitted in R with a Gibbs sampler"))
  expect_that(attr(Bout, "n.chains"), equals(3))
  expect_equal(as.character(attr(Bout, "call")), c("Bnormal", "x"))
  expect_equivalent(round(attr(Bout, "n.eff")), c(48528, 39372))
  expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1, 1))

  expect_equivalent(round(colMeans(Bout), 4), c(1.0113, 0.1563))
  expect_equivalent(round(c(hdi(Bout)), 4), c(0.9073, 1.1098, 0.0904, 0.2416))
  xx <- x * 1000
  expect_warning(Bout <- Bnormal(xx, priors=list(muMean=0, muSD=10)), # silly prior
      "Sample mean is outside the prior range")
  set.seed(123)
  Bout <- Bnormal(x, priors=list(muMean=0, muSD=10)) # informative prior for mu
  expect_equivalent(round(colMeans(Bout), 4), c(1.0112, 0.1563))
  expect_equivalent(round(c(hdi(Bout)), 4), c(0.9072, 1.1098, 0.0904, 0.2416))
  set.seed(234)
  Bout <- Bnormal(xx, priors=list(muMean=1000))  # sensible priors
  expect_equivalent(round(colMeans(Bout), 1), c(1011.4, 156.6))
  expect_equivalent(round(c(hdi(Bout)), 1), c(909.1, 1113.8, 91.2, 243.0))

})

test_that("Bnormal2 gives same answers",  {
  # Generate data
  set.seed(123)
  x <- rnorm(10, 1, 0.15)
  Bout <- Bnormal2(x, rnd.seed=345, verbose=FALSE)  # default prior
  expect_that(class(Bout), equals(c("Bwiqid", "data.frame")))
  expect_that(dim(Bout), equals(c(100002, 2)))
  expect_that(names(Bout), equals(c("mu", "sigma")))
  expect_that(attr(Bout, "header"),
    equals("Model fitted in JAGS with jagsUI::jags.basic"))
  expect_that(attr(Bout, "n.chains"), equals(3))
  expect_equal(as.character(attr(Bout, "call")), c("Bnormal2", "x", "FALSE", "345"))
  if(packageVersion("rjags") < "4.0.0") {
    expect_equivalent(round(attr(Bout, "n.eff")), c(103943, 22727 ))
    expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1, 1.003))
    expect_equivalent(round(colMeans(Bout), 4), c(1.0114, 0.1685))
    expect_equivalent(round(c(hdi(Bout)), 4), c(0.9010, 1.1222, 0.0932, 0.2665))
  } else {
    expect_equivalent(round(attr(Bout, "n.eff")), c(99939, 26041 ))
    expect_equivalent(round(attr(Bout, "Rhat"), 3), c(1, 1))
    expect_equivalent(round(colMeans(Bout), 4), c(1.0114, 0.1683))
    expect_equivalent(round(c(hdi(Bout)), 4), c(0.9002, 1.1216, 0.0930, 0.2664))
  }
  xx <- x * 1000
  expect_warning(Bout <- Bnormal2(xx, priors=list(muMean=0, muSD=10)), # silly prior
      "Sample mean is outside the prior range")
  Bout <- Bnormal2(x, priors=list(muMean=0, muSD=10), rnd.seed=345, verbose=FALSE)  # informative prior for mu, default for sigma
  if(packageVersion("rjags") < "4.0.0") {
    expect_equivalent(round(colMeans(Bout), 4), c(1.0111, 0.1687))
    expect_equivalent(round(c(hdi(Bout)), 4), c(0.9003, 1.1221, 0.0918, 0.2659))
  } else {
    expect_equivalent(round(colMeans(Bout), 4), c(1.0112, 0.1684))
    expect_equivalent(round(c(hdi(Bout)), 4), c(0.8987, 1.1202, 0.0934, 0.2677))
  }
  Bout <- Bnormal2(x, priors=list(muMean=0, muSD=10, sigmaMode=0.1, sigmaSD=0.2), rnd.seed=345, verbose=FALSE)  # informative prior for mu and sigma
  expect_equivalent(round(colMeans(Bout), 4), c(1.0112, 0.1622))
  if(packageVersion("rjags") < "4.0.0") {
    expect_equivalent(round(c(hdi(Bout)), 4), c(0.9032, 1.1156, 0.0917, 0.2473))
  } else {
    expect_equivalent(round(c(hdi(Bout)), 4), c(0.9028, 1.1144, 0.0927, 0.2474))
  }
})
