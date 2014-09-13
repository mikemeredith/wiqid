
# More test_that code for occSS series of functions

context("Single-season 2-species occupancy")

test_that("occ2sps gives right answers",  {
  require(wiqid)
  data(railSims)
  DHA <- railSims[, 1:3]
  DHB <- railSims[, 4:6]
  # Default model (no interaction)
  rail1 <- occ2sps(DHA, DHB)
  expect_that(class(rail1), equals(c("wiqid", "list"))) 
  expect_that(names(rail1), equals(c("call", "beta", "beta.vcv", "real", "logLik"))) 
  expect_that(is.call(rail1$call), is_true())

  expect_that(dim(rail1$beta), equals(c(8, 4)))
  expect_that(colnames(rail1$beta),
    equals(c("est", "SE",  "lowCI", "uppCI")))
  expect_that(rownames(rail1$beta),
    equals(c("psiA", "psiBA", "psiBa", "pA", "pB", "rA", "rBA", "rBa")))
  expect_that(round(as.vector(rail1$beta), 4), 
      equals(c(-0.3794, 0.1628, 0.1628, 0.4268, 1.5220, 0.4268, 1.5220,
          1.5220, 0.1732, 0.1597, 0.1597, 0.1753, 0.1687, 0.1753, 0.1687,
          0.1687, -0.7189, -0.1502, -0.1502,  0.0832, 1.1913, 0.0832,
          1.1913,  1.1913, -0.0399,  0.4758,  0.4758, 0.7703,  1.8527,
          0.7703,  1.8527,  1.8527)))
  expect_that(dim(rail1$real), equals(c(8, 3)))
  expect_that(colnames(rail1$real),
    equals(c("est", "lowCI", "uppCI")))
  # Check against PRESENCE results:
  expect_that(round(as.vector(rail1$real[, 1]), 4), 
      equals(c(0.4063, 0.5406, 0.5406, 0.6051, 0.8208, 0.6051, 0.8208, 0.8208)))
  expect_that(round(as.vector(rail1$real[, 2]), 4), 
      equals(c(0.3276, 0.4625, 0.4625, 0.5208, 0.767, 0.5208, 0.767, 0.767)))
  expect_that(round(as.vector(rail1$real[, 3]), 4), 
      equals(c(0.4900, 0.6168, 0.6168, 0.6836, 0.8644, 0.6836, 0.8644, 0.8644)))
  expect_that(round(AIC(rail1), 4), 
      equals(911.0552))

  # Model with full interaction
  rail2 <- occ2sps(DHA, DHB, modelSpec=224)
  # Check against PRESENCE results:
  expect_that(round(as.vector(rail2$real[, 1]), 4), 
      equals(c(0.4093, 0.7844, 0.3723, 0.7258, 0.7866, 0.5663, 0.8825, 0.7925)))
  expect_that(round(as.vector(rail2$real[, 2]), 4), 
      equals(c(0.3294, 0.6642, 0.2751, 0.5559, 0.6853, 0.4679, 0.7943, 0.6625)))
  expect_that(round(as.vector(rail2$real[, 3]), 4), 
      equals(c(0.4942, 0.8699, 0.4811, 0.8484, 0.8618, 0.6597, 0.9359, 0.8814)))
  expect_that(round(AIC(rail2), 4), 
      equals(890.4437))
})


