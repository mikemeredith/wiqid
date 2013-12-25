
# test_that code for survCJS series of functions

# library(testthat)
# test_file("test-cjs.R")
# test_file("wiqid/inst/tests/test-cjs.R")

context("CJS survival models")


test_that("survCJS gives right answers",  {
  # Data set (Lereton et al dippers)
  dippers <- matrix(c(
    11,2,0,0,0,0,9,
    0,24,1,0,0,0,35,
    0,0,34,2,0,0,42,
    0,0,0,45,1,2,32,
    0,0,0,0,51,0,37,
    0,0,0,0,0,52,46),nrow=6,byrow=TRUE)

  res <- survCJS(dippers)  # default is a phi(.) p(.) model
  expect_that(class(res), equals(c("wiqid", "list")))
  expect_that(names(res), equals(c("call", "beta", "beta.vcv", "real", "logLik")))
  expect_that(colnames(res$real), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(res$beta), equals(c("phi: (Intercept)", "p: (Intercept)")))
  expect_that(round(as.vector(t(res$real[6:7,])), 4), 
      equals(c(0.5602, 0.5105, 0.6088, 0.9026, 0.8305, 0.9460))) # MARK output
  expect_that(round(AIC(res), 4), equals(670.8377)) # MARK output

  res <- survCJS(dippers, phi ~ time) 
  expect_that(rownames(res$beta), equals(c("phi: (Intercept)", "phi: time2",
          "phi: time3", "phi: time4", "phi: time5", "phi: time6",
          "p: (Intercept)")))
  expect_that(round(as.vector(t(res$real[1:7,])), 4), 
      equals(c(0.6258, 0.3965, 0.8098,
               0.4542, 0.3295, 0.5849,
               0.4784, 0.3669, 0.5921,
               0.6244, 0.5079, 0.7281,
               0.6079, 0.4970, 0.7088,
               0.5833, 0.4688, 0.6895,
               0.9021, 0.8286, 0.9461)))# MARK output
  expect_that(round(AIC(res), 4), equals(673.7301)) # MARK output

  dd <- data.frame(flood = c(0, 1, 1, 0, 0, 0))
  res <- survCJS(dippers, phi ~ flood, data=dd)
  expect_that(rownames(res$beta), equals(c("phi: (Intercept)", "phi: flood",
          "p: (Intercept)")))
  expect_that(round(as.vector(t(res$real[c(1,2,7),])), 4), 
      equals(c(0.6071, 0.5451, 0.6658,
               0.4688, 0.3858, 0.5537,
               0.8998, 0.8262, 0.9443)))   # MARK output
  expect_that(round(AIC(res), 4), equals(666.1028)) # MARK output

  res <- survCJS(dippers, phi ~ flood, data=dd, ci=0.85)
  expect_that(round(as.vector(t(res$real[c(1,2,7),])), 4), 
      equals(c(0.6071, 0.5618, 0.6506,
               0.4688, 0.4074, 0.5312,
               0.8998, 0.8491, 0.9348)))
  expect_that(round(AIC(res), 4), equals(666.1028)) # MARK output
}  )
# .........................................................................



