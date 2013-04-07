
# test_that code for closedCap series of functions

# library(testthat)
# test_file("test-occSS.R")

context("Closed captures")

test_that("closedCapM0 gives right answers",  {
  # Rabbit data from Edwards and Eberhart (1967)
  freq2 <- c(43, 16, 8, 6, 0, 2, 1) ; t2 <- 18 
  res <- closedCapM0(freq2, t2)
  expect_that(colnames(res), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(res), equals(c("Nhat", "pHat")))
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(96.2589, 0.0820, 86.4162, 0.0663, 115.4025, 0.1010)))
  expect_that(round(attr(res, "AIC"), 4), equals(379.5941))
      # These are almost the same as the values returned by MARK.
  res <- closedCapM0(freq2, t2, ci=0.85)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(96.2589, 0.0820, 88.4286, 0.0701, 109.0225, 0.0956)))

   # All zeros:
  res <- closedCapM0(rep(0, 18))
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 6)))
  expect_that(round(attr(res, "AIC"), 4), equals(NA_real_))
  # Lots of captures, no recaptures
  res <- closedCapM0(30, 18)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 6)))
  expect_that(round(attr(res, "AIC"), 4), equals(NA_real_))
  # Just 1 animal recaptured
  res <- closedCapM0(c(0,1), 18)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 6)))
  expect_that(round(attr(res, "AIC"), 4), equals(NA_real_))
  # Just 2 animals recaptured
  res <- closedCapM0(c(0,2), 18)
  expect_that(signif(as.vector(res)[-5], 5), # 5th value is nonsense
      is_equivalent_to(c(2.0000e+00, 1.1111e-01,  2.0000e+00, 4.2336e-02,
          2.6114e-01)))
  expect_that(round(attr(res, "AIC"), 4), equals(27.7296))
}  )
# .........................................................................

test_that("closedCapMh2 gives right answers",  {
  # Rabbit data from Edwards and Eberhart (1967)
  freq2 <- c(43, 16, 8, 6, 0, 2, 1) ; t2 <- 18 
  res <- closedCapMh2(freq2, t2)
  expect_that(colnames(res), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(res), equals(c("Nhat", "piHat","p1hat", "p2hat")))
  expect_that(signif(as.vector(res), 4),
      is_equivalent_to(c(135.5, 0.1548,   0.1795,   0.03602,  93.85,   0.0536,   0.1027,
                         0.01222, 274.2,   0.3718,   0.2949,   0.1015)))
  expect_that(round(attr(res, "AIC"), 4), equals(369.6155))
      # These are almost the same as the values returned by MARK.
  res <- closedCapMh2(freq2, t2, 0.85)
  expect_that(signif(as.vector(res), 4),
      is_equivalent_to(c(135.5,   0.1548,  0.1795,   0.03602, 100.6,   0.07179,   0.1197,
   0.01632, 220.0,  0.3024,   0.2604,   0.07765)))
   # All zeros:
  res <- closedCapMh2(rep(0, 18))
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 12)))
  expect_that(round(attr(res, "AIC"), 4), equals(NA_real_))
  # Lots of captures, no recaptures
  res <- closedCapMh2(30, 18)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 12)))
  expect_that(round(attr(res, "AIC"), 4), equals(NA_real_))
  # Just 1 animal recaptured
  res <- closedCapMh2(c(0,1), 18)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 12)))
  expect_that(round(attr(res, "AIC"), 4), equals(NA_real_))
  # Kanha tiger data
  res <- closedCapMh2(c(10,6,6,2,2), 10)
  expect_that(signif(as.vector(res), 4),
      is_equivalent_to(c(3.152e+01, 4.920e-01, 2.644e-01, 1.061e-01, 2.629e+01,
        1.477e-02, 1.130e-01, 5.821e-03, 1.294e+02, 9.843e-01, 5.036e-01, 7.066e-01)))
  expect_that(round(attr(res, "AIC"), 4), equals(158.6416))
}  )
# .........................................................................


test_that("closedCapMhJK gives right answers",  {
  # Rabbit data from Edwards and Eberhart (1967)
  freq2 <- c(43, 16, 8, 6, 0, 2, 1) ; t2 <- 18 
  res <- closedCapMhJK(freq2, t2)
  expect_that(colnames(res), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(res), equals(c("Nhat", "pHat")))
  expect_that(attr(res, "AIC"), equals(NA))
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(143.8743, 0.0548, 112.9113, 0.0393, 200.8106, 0.0699)))
      # These are almost the same as the rounded values returned by CAPTURE.
  res <- closedCapMhJK(freq2, t2, 0.85)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(143.8743, 0.0548, 119.3915, 0.0433, 182.1710, 0.0661)))

   # All zeros:
  res <- closedCapMhJK(rep(0, 18))
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 6)))
  expect_that(attr(res, "AIC"), equals(NA))
  # Lots of captures, no recaptures
  res <- closedCapMhJK(30, 18)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 6)))
  # Just 1 animal recaptured
  res <- closedCapMhJK(c(0,1), 18)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 6)))
  # Kanha tiger data
  res <- closedCapMhJK(c(10,6,6,2,2), 10)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(33.3188, 0.1741, 27.8792, 0.1064, 54.5044, 0.2080)))
}  )
# .........................................................................
