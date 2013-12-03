
# More test_that code for occSS series of functions

# library(testthat)
# test_file("test-occSScov.R")
# test_file("wiqid/inst/tests/test-occSScov.R")

context("Single-season occupancy + covars")

test_that("occSScovSite gives right answers",  {
  # Data set (weta)
  require(wiqid)
  data(weta)
  DH <- weta[, 1:5]
  y <- rowSums(DH, na.rm=TRUE)
  n <- rowSums(!is.na(DH))
  weta.covs <- weta[, 6, drop=FALSE]
  weta1 <- occSScovSite(y, n)

  expect_that(class(weta1), equals(c("occupancy", "list"))) 
  expect_that(names(weta1), equals(c("call", "beta", "real", "logLik"))) 
  expect_that(is.call(weta1$call), is_true())

  expect_that(dim(weta1$beta), equals(c(2, 4)))
  expect_that(colnames(weta1$beta),
    equals(c("est", "SE",  "lowCI", "uppCI")))
  expect_that(rownames(weta1$beta),
    equals(c("psi: (Intercept)", "p: (Intercept)")))
  expect_that(round(as.vector(weta1$beta), 4), 
      equals(c(0.4751, -0.6218,  0.3746,  0.2364, -0.2590, -1.0851,  1.2093, -0.1585)))

  expect_that(dim(weta1$real), equals(c(144, 3)))
  expect_that(colnames(weta1$real),
    equals(c("est", "lowCI", "uppCI")))
  # Check against PRESENCE results:
  expect_that(round(as.vector(weta1$real[1, ]), 4), 
      equals(c(0.6166, 0.4356, 0.7702)))
  expect_that(round(as.vector(weta1$real[73, ]), 4), 
      equals(c(0.3494, 0.2525, 0.4604)))
  expect_that(round(AIC(weta1), 4), 
      equals(265.7872))

  weta1a <- occSScovSite(y, n, ci=0.85)
  expect_that(round(as.vector(weta1a$real[1, ]), 4), 
      equals(c(0.6166, 0.4840, 0.7339)))
  expect_that(round(as.vector(weta1a$real[73, ]), 4), 
      equals(c(0.3494, 0.2765, 0.4301)))

  weta2 <- occSScovSite(y, n, psi ~ Browsed, data=weta.covs)

  expect_that(names(weta2), equals(c("call", "beta", "real", "logLik"))) 

  expect_that(dim(weta2$beta), equals(c(3, 4)))
  expect_that(colnames(weta2$beta),
    equals(c("est", "SE",  "lowCI", "uppCI")))
  expect_that(rownames(weta2$beta),
    equals(c("psi: (Intercept)", "psi: Browsed", "p: (Intercept)")))
  expect_that(round(as.vector(weta2$beta), 4), 
      equals(c(0.5196, 0.6167, -0.6223,  0.4188,  0.3627,  0.2367, -0.3012,
              -0.0942, -1.0861, 1.3403,  1.3276, -0.1584)))

  expect_that(dim(weta2$real), equals(c(144, 3)))
  expect_that(colnames(weta2$real),
    equals(c("est", "lowCI", "uppCI")))
  # Check against PRESENCE results:
  expect_that(round(as.vector(weta2$real[1, ]), 4), 
      equals(c(0.7594, 0.4660, 0.9194)))
  expect_that(round(as.vector(weta2$real[4, ]), 4), 
      equals(c(0.4810, 0.2843, 0.6837)))
  expect_that(round(as.vector(weta2$real[73, ]), 4), 
      equals(c(0.3493,  0.2524, 0.4605)))
  expect_that(round(AIC(weta2), 4), 
      equals(264.2643))

  weta3 <- occSScovSite(y, n, list(psi ~ Browsed, p ~ Browsed), data=weta.covs)

  expect_that(names(weta3), equals(c("call", "beta", "real", "logLik"))) 

  expect_that(dim(weta3$beta), equals(c(4, 4)))
  expect_that(colnames(weta3$beta),
    equals(c("est", "SE",  "lowCI", "uppCI")))
  expect_that(rownames(weta3$beta),
    equals(c("psi: (Intercept)", "psi: Browsed", "p: (Intercept)", "p: Browsed")))
  expect_that(round(as.vector(weta3$beta), 4), 
      equals(c( 0.5180,  0.6030, -0.6260,  0.0149,  0.4168,  0.4236,  0.2446,
                0.2446, -0.2989, -0.2273, -1.1054, -0.4645,  1.3349,  1.4333,
               -0.1465,  0.4943)))
  expect_that(dim(weta3$real), equals(c(144, 3)))
  expect_that(colnames(weta3$real),
    equals(c("est", "lowCI", "uppCI")))
  # Check against PRESENCE results:
  expect_that(round(as.vector(weta3$real[1, ]), 4), 
      equals(c(0.7565, 0.4438, 0.9237))) # PRESENCE last value is 0.9236
  expect_that(round(as.vector(weta3$real[4, ]), 4), 
      equals(c(0.4839, 0.2692, 0.7048))) 
  expect_that(round(as.vector(weta3$real[73, ]), 4), 
      equals(c(0.3519, 0.2309, 0.4954)))
  expect_that(round(as.vector(weta3$real[76, ]), 4), 
      equals(c(0.3452, 0.2000, 0.5264))) # PRESENCE last value is 0.5263
  expect_that(round(AIC(weta3), 4), 
      equals(266.2606))
})


test_that("occSScov gives right answers",  {
  # Data set (weta)
  require(wiqid)
  data(weta)
  DH <- weta[, 1:5]
  weta.covs <- weta[, 6, drop=FALSE]
  weta.obs <- weta[, 7:11]
  weta.dat <- as.list(weta.covs)
  weta.dat$obs <- as.matrix(weta.obs)

  weta4 <- occSScov(DH)
  expect_that(class(weta4), equals(c("occupancy", "list"))) 
  expect_that(names(weta4), equals(c("call", "beta", "real", "logLik"))) 
  expect_that(is.call(weta4$call), is_true())
  expect_that(dim(weta4$beta), equals(c(2, 4)))
  expect_that(colnames(weta4$beta),
    equals(c("est", "SE",  "lowCI", "uppCI")))
  expect_that(rownames(weta4$beta),
    equals(c("psi: (Intercept)", "p: (Intercept)")))
  expect_that(round(as.vector(weta4$beta), 4), 
      equals(c(0.4751, -0.6218,  0.3746,  0.2364, -0.2590, -1.0851,  1.2093, -0.1585)))

  expect_that(dim(weta4$real), equals(c(334, 3)))
  expect_that(colnames(weta4$real),
    equals(c("est", "lowCI", "uppCI")))
  # Check against PRESENCE results:
  expect_that(round(as.vector(weta4$real[1, ]), 4), 
      equals(c(0.6166, 0.4356, 0.7702)))
  expect_that(round(as.vector(weta4$real[73, ]), 4), 
      equals(c(0.3494, 0.2525, 0.4604)))
  expect_that(round(AIC(weta4), 4), 
      equals(265.7872))
  weta4a <- occSScov(DH, ci=0.85)
  expect_that(round(as.vector(weta4a$real[1, ]), 4), 
      equals(c(0.6166, 0.4840, 0.7339)))
  expect_that(round(as.vector(weta4a$real[73, ]), 4), 
      equals(c(0.3494, 0.2765, 0.4301)))


  weta5 <- occSScov(DH, psi= ~ Browsed, data=weta.dat)

  expect_that(names(weta5), equals(c("call", "beta", "real", "logLik"))) 

  expect_that(dim(weta5$beta), equals(c(3, 4)))
  expect_that(colnames(weta5$beta),
    equals(c("est", "SE",  "lowCI", "uppCI")))
  expect_that(rownames(weta5$beta),
    equals(c("psi: (Intercept)", "psi: Browsed", "p: (Intercept)")))
  expect_that(round(as.vector(weta5$beta), 4), 
      equals(c( 0.5196,  0.6167, -0.6223,  0.4188,  0.3627,  0.2367, -0.3012,
               -0.0942, -1.0861,  1.3403,  1.3276, -0.1584)))

  expect_that(dim(weta5$real), equals(c(334, 3)))
  expect_that(colnames(weta5$real),
    equals(c("est", "lowCI", "uppCI")))
  # Check against PRESENCE results:
  expect_that(round(as.vector(weta5$real[1, ]), 4), 
      equals(c(0.7594, 0.4660, 0.9194)))
  expect_that(round(as.vector(weta5$real[4, ]), 4), 
      equals(c(0.4810, 0.2843, 0.6837)))
  expect_that(round(as.vector(weta5$real[73, ]), 4), 
      equals(c(0.3493,  0.2524, 0.4605)))
  expect_that(round(as.vector(weta5$real[334, ]), 4), 
      equals(c(0.3493,  0.2524, 0.4605)))
  expect_that(round(AIC(weta5), 4), 
      equals(264.2643))

  weta6 <- occSScov(DH, psi= ~ Browsed, p= ~ obs, data=weta.dat)

  expect_that(names(weta6), equals(c("call", "beta", "real", "logLik"))) 

  expect_that(dim(weta6$beta), equals(c(5, 4)))
  expect_that(colnames(weta6$beta),
    equals(c("est", "SE",  "lowCI", "uppCI")))
  expect_that(rownames(weta6$beta),
    equals(c("psi: (Intercept)", "psi: Browsed", "p: (Intercept)", 
             "p: obsB",  "p: obsC")))
  expect_that(round(as.vector(weta6$beta), 4), 
      equals(c( 0.5119,  0.5935, -1.2453,  0.7498,  1.0294,  0.4061,  0.3532,
                0.3615,  0.4427,  0.4324, -0.2841, -0.0988, -1.9537, -0.1179,
                0.1819,  1.3079,  1.2858, -0.5368,  1.6175,  1.8770)))

  expect_that(dim(weta6$real), equals(c(334, 3)))
  expect_that(colnames(weta6$real),
    equals(c("est", "lowCI", "uppCI")))
  # Check against PRESENCE results:
  expect_that(round(as.vector(weta6$real[1, ]), 4), 
      equals(c(0.7536, 0.4723, 0.9127)))
  expect_that(round(as.vector(weta6$real[4, ]), 4), 
      equals(c(0.4847, 0.2865, 0.6877)))
  expect_that(round(as.vector(weta6$real[73, ]), 4), 
      equals(c(0.2235, 0.1241, 0.3689))) # PRESENCE 2nd value is 0.1242
  expect_that(round(as.vector(weta6$real[94, ]), 4), 
      equals(c(0.3786, 0.2383, 0.5426)))
  expect_that(round(as.vector(weta6$real[334, ]), 4), 
      equals(c(0.4462, 0.2980, 0.6047)))

  expect_that(round(AIC(weta6), 4), 
      equals(262.0421))
})

