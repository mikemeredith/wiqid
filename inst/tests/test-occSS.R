
# test_that code for occSS series of functions

# library(testthat)
# test_file("test-occSS.R")

context("Single-season occupancy")

test_that("occSS0 gives right answers",  {
  # Data set (Blue Ridge Salamanders)
  BRS <- matrix(0, 39, 5)
  BRS[c(4,10,13,14,41,42,43,82,83,84,85,86,87,89,90,92,118,121,128,
          129,130,131,153,154,157,167,168,169,194,195)] <- 1
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  expect_that(colnames(occSS0(y, n)), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(occSS0(y, n)), equals(c("psiHat", "pHat")))
  expect_that(round(as.vector(occSS0(y, n)), 4), 
      is_equivalent_to(c(0.5946, 0.2587, 0.3512, 0.1622, 0.7990, 0.3863)))
  expect_that(round(attr(occSS0(y, n), "AIC"), 4), equals(165.7586))
      # These are the values returned by PRESENCE
  # Put in some NAs
  BRS[c(6,167,130,123,89,154,32,120,127,147)] <- NA
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  expect_that(round(as.vector(occSS0(y, n)), 4), 
      is_equivalent_to(c(0.5861, 0.2445, 0.3238, 0.1422, 0.8073, 0.3872)))
  expect_that(round(attr(occSS0(y, n), "AIC"), 4), equals(149.7430))
  # Put in a row of NAs
  BRS[3,] <- NA
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  expect_that(round(as.vector(occSS0(y, n)), 4), 
      is_equivalent_to(c(0.5558, 0.2531, 0.3095, 0.1475, 0.7775, 0.3989)))
  expect_that(round(attr(occSS0(y, n), "AIC"), 4), equals(144.1232))
  # Put in a column of NAs
  BRS[, 3] <- NA
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  expect_that(round(as.vector(occSS0(y, n)), 4), 
      is_equivalent_to(c(0.3798, 0.3160, 0.1918, 0.1644, 0.6124, 0.5204)))
  expect_that(round(attr(occSS0(y, n), "AIC"), 4), equals(101.1297))
  # All ones:
  expect_that(round(as.vector(occSS0(n, n)), 4), 
      is_equivalent_to(c(1, 1, NA_real_, NA_real_, NA_real_, NA_real_)))
  expect_that(round(attr(occSS0(n, n), "AIC"), 4), equals(4))
  # All zeros:
  expect_that(round(as.vector(occSS0(rep(0, length(n)), n)), 4), 
      is_equivalent_to(c(0.0000, 0.0078, 0.0000, 0.0000, 1.0000, 1.0000)))
  expect_that(round(attr(occSS0(rep(0, length(n)), n), "AIC"), 4), equals(4))
  # All NAs:
  expect_that(as.vector(occSS0(rep(0, length(n)), rep(0, length(n)))), 
      is_equivalent_to(rep(NA_real_, 6)))
  expect_that(attr(occSS0(rep(0, length(n)), rep(0, length(n))), "AIC"),
      equals(NA_real_))
}  )
# .........................................................................

test_that("occSS.t gives right answers with time=FALSE",  {
  # Data set (Blue Ridge Salamanders)
  BRS <- matrix(0, 39, 5)
  BRS[c(4,10,13,14,41,42,43,82,83,84,85,86,87,89,90,92,118,121,128,
          129,130,131,153,154,157,167,168,169,194,195)] <- 1
  res <- occSS.t(BRS, time=FALSE)
  expect_that(colnames(res), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(res), equals(c("psiHat", "pHat")))
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.5946, 0.2587, 0.3512, 0.1622, 0.7990, 0.3863)))
  expect_that(round(attr(res, "AIC"), 4), equals(165.7586))
      # These are the values returned by PRESENCE
  # Put in some NAs
  BRS[c(6,167,130,123,89,154,32,120,127,147)] <- NA
  res <- occSS.t(BRS, time=FALSE)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.5861, 0.2445, 0.3238, 0.1422, 0.8073, 0.3872)))
  expect_that(round(attr(res, "AIC"), 4), equals(149.7430))
  # Put in a row of NAs
  BRS[3,] <- NA
  res <- occSS.t(BRS, time=FALSE)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.5558, 0.2531, 0.3095, 0.1475, 0.7775, 0.3989)))
  expect_that(round(attr(res, "AIC"), 4), equals(144.1232))
  # Put in a column of NAs
  BRS[, 3] <- NA
  res <- occSS.t(BRS, time=FALSE)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.3798, 0.3160, 0.1918, 0.1644, 0.6124, 0.5204)))
  expect_that(round(attr(res, "AIC"), 4), equals(101.1297))
  # All ones:
  tst <- matrix(1, 39, 5)
  res <- occSS.t(tst, time=FALSE)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(1, 1, NA_real_, NA_real_, NA_real_, NA_real_)))
  expect_that(round(attr(res, "AIC"), 4), equals(4))
  # All zeros:
  tst <- matrix(0, 39, 5)
  res <- occSS.t(tst, time=FALSE)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 6)))
  expect_that(attr(res, "AIC"), equals(NA_real_))
  # All NAs:
  tst <- matrix(NA, 39, 5)
  res <- occSS.t(tst, time=FALSE)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 6)))
  expect_that(attr(res, "AIC"),
      equals(NA_real_))
}  )
# ....................................................................

test_that("occSS.t gives right answers with time=default",  {
  # Data set (Blue Ridge Salamanders)
  BRS <- matrix(0, 39, 5)
  BRS[c(4,10,13,14,41,42,43,82,83,84,85,86,87,89,90,92,118,121,128,
          129,130,131,153,154,157,167,168,169,194,195)] <- 1
  res <- occSS.t(BRS)
  expect_that(colnames(res), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(res),
       equals(c("psiHat", "p1Hat", "p2Hat", "p3Hat", "p4Hat", "p5Hat")))
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.5799, 0.1769, 0.1327, 0.3980, 0.3537, 0.2653, 0.3490,
          0.0644, 0.0415, 0.1998, 0.1712, 0.1156, 0.7804, 0.4013, 0.3506, 0.6364,
          0.5920, 0.4993)))
  expect_that(round(attr(res, "AIC"), 4), equals(167.7144))
      # These are the values returned by PRESENCE

  # Put in some NAs
  BRS[c(6,167,130,123,89,154,32,120,127,147)] <- NA
  res <- occSS.t(BRS)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.5637, 0.1930, 0.1365, 0.3812, 0.3450, 0.2383, 0.3223,
        0.0690, 0.0421, 0.1773, 0.1424, 0.0938, 0.7783, 0.4354, 0.3621, 0.6378,
        0.6257, 0.4861)))
  expect_that(round(attr(res, "AIC"), 4), equals(153.1581))
  # Put in a row of NAs
  BRS[3,] <- NA
  res <- occSS.t(BRS)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.5316, 0.2107, 0.0990, 0.4166, 0.3596, 0.2604, 0.3067,
        0.0758, 0.0238, 0.1952, 0.1514, 0.1031, 0.7444, 0.4650, 0.3308, 0.6778,
        0.6387, 0.5188)))
  expect_that(round(attr(res, "AIC"), 4), equals(145.6360))
  # Put in a column of NAs
  BRS[, 3] <- NA
  res <- occSS.t(BRS)
  expect_that(round(res[, 1], 4), 
      is_equivalent_to(c(0.3579, 0.3017, 0.1471, NA_real_, 0.5434, 0.3969)))
  expect_that(round(as.vector(res[-4, 2:3]), 4), 
      is_equivalent_to(c(0.1882, 0.1078, 0.0349, 0.2149, 0.1511, 0.5726, 0.6070,
        0.4511, 0.8380, 0.7086)))
  expect_that(res[4, ], 
      is_equivalent_to(rep(NA_real_, 3)))
  expect_that(round(attr(res, "AIC"), 4), equals(102.4882))
  # All ones:
  tst <- matrix(1, 39, 5)
  res <- occSS.t(tst)
  expect_that(round(as.vector(res[,1]), 4), 
      is_equivalent_to(rep(1, 6)))
  expect_that(as.vector(res[, 2:3]), 
      is_equivalent_to(rep(NA_real_, 12)))
  expect_that(round(attr(res, "AIC"), 4), equals(12))
  # All zeros:
  tst <- matrix(0, 39, 5)
  res <- occSS.t(tst)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 18)))
  expect_that(attr(res, "AIC"), equals(NA_real_))
  # All NAs:
  tst <- matrix(NA, 39, 5)
  res <- occSS.t(tst)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 18)))
  expect_that(attr(res, "AIC"),
      equals(NA_real_))
}  )
# ..............................................................

test_that("occSS.T gives right answers",  {
  # Data set (Blue Ridge Salamanders)
  BRS <- matrix(0, 39, 5)
  BRS[c(4,10,13,14,41,42,43,82,83,84,85,86,87,89,90,92,118,121,128,
          129,130,131,153,154,157,167,168,169,194,195)] <- 1
  res <- occSS.T(BRS)
  expect_that(colnames(res), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(res),
       equals(c("psiHat", "beta0", "beta1")))
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.5899, -1.4732, 0.2056, 0.3505, -2.3376,
          -0.0954, 0.7931, -0.6088, 0.5067)))
  expect_that(round(attr(res, "AIC"), 4), equals(165.9228))

  # Put in some NAs
  BRS[c(6,167,130,123,89,154,32,120,127,147)] <- NA
  res <- occSS.T(BRS)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.5800, -1.4246, 0.1524, 0.3234, -2.3391,
          -0.1632, 0.7997, -0.5101, 0.4680)))
  expect_that(round(attr(res, "AIC"), 4), equals(150.8355))
  # Put in a row of NAs
  BRS[3,] <- NA
  res <- occSS.T(BRS)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.5494, -1.4417, 0.1827, 0.3087, -2.3797,
          -0.1435, 0.7689, -0.5036, 0.5089)))
  expect_that(round(attr(res, "AIC"), 4), equals(144.8956))
  # Put in a column of NAs
  BRS[, 3] <- NA
  res <- occSS.T(BRS)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.3666, -1.2386, 0.2549, 0.1889, -2.3419,
          -0.1262, 0.5899, -0.1353, 0.6361)))
  expect_that(round(attr(res, "AIC"), 4), equals(101.3259))
  # All ones:
  tst <- matrix(1, 39, 5)
  res <- occSS.T(tst)
  expect_that(round(as.vector(res[1,1]), 4), 
      is_equivalent_to(1))
  expect_that(as.vector(res)[-1], 
      is_equivalent_to(rep(NA_real_, 8)))
  expect_that(round(attr(res, "AIC"), 4), equals(NA_real_))
  # All zeros:
  tst <- matrix(0, 39, 5)
  res <- occSS.T(tst)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 9)))
  expect_that(attr(res, "AIC"), equals(NA_real_))
  # All NAs:
  tst <- matrix(NA, 39, 5)
  res <- occSS.T(tst)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 9)))
  expect_that(attr(res, "AIC"),
      equals(NA_real_))
}  )
# .............................................

test_that("occSSnr gives right answers",  {
  # Data set (Blue Ridge Salamanders)
  BRS <- matrix(0, 39, 5)
  BRS[c(4,10,13,14,41,42,43,82,83,84,85,86,87,89,90,92,118,121,128,
          129,130,131,153,154,157,167,168,169,194,195)] <- 1
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  res <- occSSrn(y, n)
  expect_that(colnames(res), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(res), equals(c("psiHat", "lambdaHat", "rHat")))
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.6810, 1.1425, 0.1475, 0.3772, 0.4735, 0.0581, 0.9365,
          2.7568, 0.3266)))
  expect_that(round(attr(res, "AIC"), 4), equals(164.0016))
      # These are the values returned by PRESENCE
  # Put in some NAs
  BRS[c(6,167,130,123,89,154,32,120,127,147)] <- NA
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  res <- occSSrn(y, n)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.6739, 1.1204, 0.1383, 0.3402, 0.4159, 0.0473, 0.9511,
          3.0186, 0.3414)))
  expect_that(round(attr(res, "AIC"), 4), equals(148.5542))
  # Put in a row of NAs
  BRS[3,] <- NA
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  res <- occSSrn(y, n)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.6349, 1.0076, 0.1510, 0.3244, 0.3922, 0.0545, 0.9249,
          2.5884, 0.3540)))
  expect_that(round(attr(res, "AIC"), 4), equals(142.9799))
  # Put in a column of NAs
  BRS[, 3] <- NA
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  res <- occSSrn(y, n)
  expect_that(round(as.vector(res), 4), 
      is_equivalent_to(c(0.4030, 0.5158, 0.2464, 0.1915, 0.2126, 0.0957, 0.7139,
          1.2513, 0.5026)))
  expect_that(round(attr(res, "AIC"), 4), equals(101.3649))
  # All ones:
  y <- n <- rep(5, 39)
  res <- occSSrn(y, n)
  expect_that(round(as.vector(res[, 1]), 4), 
      is_equivalent_to(c(1.0000, 26.9883,  0.8662)))
  expect_that(as.vector(res[, 2:3]), 
      is_equivalent_to(rep(NA_real_, 6)))
  expect_that(round(attr(res, "AIC"), 4), equals(4))
  # All zeros:
  n <- rep(5, 39)
  y <- rep(0, 39)
  res <- occSSrn(y, n)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 9)))
  expect_that(attr(res, "AIC"), equals(NA_real_))
  # All NAs:
  y <- n <- rep(0, 39)
  res <- occSSrn(y, n)
  expect_that(as.vector(res), 
      is_equivalent_to(rep(NA_real_, 9)))
  expect_that(attr(res, "AIC"),
      equals(NA_real_))
}  )
# ....................................................................

