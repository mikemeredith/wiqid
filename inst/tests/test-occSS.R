
# test_that code for occSS series of functions

# library(testthat)
# test_file("test-occSS.R")

context("Single-season occupancy")

test_that("occSS0 gives right answers",  {
  # Data set (Blue Ridge Salamanders)
  require(wiqid)
  data(salamanders)
  BRS <- salamanders
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  brs1 <- occSS0(y, n)
  expect_that(class(brs1), equals(c("occupancy", "list")))
  expect_that(names(brs1), equals(c("call", "beta", "real", "AIC")))
  expect_that(is.call(brs1$call), is_true())
  expect_that(colnames(brs1$real), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(brs1$real), equals(c("psiHat", "pHat")))
  expect_that(round(as.vector(brs1$real), 4), 
      equals(c(0.5946, 0.2587, 0.3512, 0.1622, 0.7990, 0.3863)))
  expect_that(round(brs1$AIC, 4), equals(165.7586))
      # These are the values returned by PRESENCE
  # Put in some NAs
  BRS[c(6,167,130,123,89,154,32,120,127,147)] <- NA
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  brs2 <-  occSS0(y, n)
  expect_that(round(as.vector(brs2$real), 4), 
      equals(c(0.5861, 0.2445, 0.3238, 0.1422, 0.8073, 0.3872)))
  expect_that(round(brs2$AIC, 4), equals(149.7430))
  # Put in a row of NAs
  BRS[3,] <- NA
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  brs3 <- occSS0(y, n)
  expect_that(round(as.vector(brs3$real), 4), 
      equals(c(0.5558, 0.2531, 0.3095, 0.1475, 0.7775, 0.3989)))
  expect_that(round(brs3$AIC, 4), equals(144.1232))
  # Put in a column of NAs
  BRS[, 3] <- NA
  n <- rowSums(!is.na(BRS))
  y <- rowSums(BRS > 0, na.rm=TRUE)
  brs4 <-  occSS0(y, n)
  expect_that(round(as.vector(brs4$real), 4), 
      equals(c(0.3798, 0.3160, 0.1918, 0.1644, 0.6124, 0.5204)))
  expect_that(round(brs4$AIC, 4), equals(101.1297))
  # All ones:
  brs5 <- occSS0(n, n)
  expect_that(round(as.vector(brs5$real), 4), 
      equals(c(1, 1, NA_real_, NA_real_, NA_real_, NA_real_)))
  expect_that(round(brs5$AIC, 4), equals(4))
  # All zeros:
  brs6 <- occSS0(rep(0, length(n)), n)
  expect_that(round(as.vector(brs6$real), 4), 
      equals(c(0.0000, 0.0078, 0.0000, 0.0000, 1.0000, 1.0000)))
  expect_that(round(brs6$AIC, 4), equals(4))
  # All NAs:
  brs7 <- occSS0(rep(0, length(n)), rep(0, length(n)))
  expect_that(as.vector(brs7$real), 
      equals(rep(NA_real_, 6)))
  expect_that(brs7$AIC, equals(NA_real_))
}  )
# .........................................................................


test_that("occSS.t gives right answers",  {
  # Data set (Blue Ridge Salamanders)
  require(wiqid)
  data(salamanders)
  BRS <- salamanders
  res <- occSS.t(BRS)
  expect_that(class(res), equals(c("occupancy", "list")))
  expect_that(names(res), equals(c("call", "beta", "real", "AIC")))
  expect_that(is.call(res$call), is_true())
  expect_that(colnames(res$real), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(res$real),
       equals(c("psiHat", "p1Hat", "p2Hat", "p3Hat", "p4Hat", "p5Hat")))
  expect_that(round(as.vector(res$real), 4), 
      equals(c(0.5799, 0.1769, 0.1327, 0.3980, 0.3537, 0.2653, 0.3490,
          0.0644, 0.0415, 0.1998, 0.1712, 0.1156, 0.7804, 0.4013, 0.3506, 0.6364,
          0.5920, 0.4993)))
  expect_that(round(res$AIC, 4), equals(167.7144))
      # These are the values returned by PRESENCE

  # Put in some NAs
  BRS[c(6,167,130,123,89,154,32,120,127,147)] <- NA
  res <- occSS.t(BRS)
  expect_that(round(as.vector(res$real), 4), 
      equals(c(0.5637, 0.1930, 0.1365, 0.3812, 0.3450, 0.2383, 0.3223,
        0.0690, 0.0421, 0.1773, 0.1424, 0.0938, 0.7783, 0.4354, 0.3621, 0.6378,
        0.6257, 0.4861)))
  expect_that(round(res$AIC, 4), equals(153.1581))
  # Put in a row of NAs
  BRS[3,] <- NA
  res <- occSS.t(BRS)
  expect_that(round(as.vector(res$real), 4), 
      equals(c(0.5316, 0.2107, 0.0990, 0.4166, 0.3596, 0.2604, 0.3067,
        0.0758, 0.0238, 0.1952, 0.1514, 0.1031, 0.7444, 0.4650, 0.3308, 0.6778,
        0.6387, 0.5188)))
  expect_that(round(res$AIC, 4), equals(145.6360))
  # Put in a column of NAs
  BRS[, 3] <- NA
  res <- occSS.t(BRS)
  expect_that(round(res$real[, 1], 4), 
      is_equivalent_to(c(0.3579, 0.3017, 0.1471, NA_real_, 0.5434, 0.3969)))
  expect_that(round(as.vector(res$real[-4, 2:3]), 4), 
      is_equivalent_to(c(0.1882, 0.1078, 0.0349, 0.2149, 0.1511, 0.5726, 0.6070,
        0.4511, 0.8380, 0.7086)))
  expect_that(res$real[4, ], 
      is_equivalent_to(rep(NA_real_, 3)))
  expect_that(round(res$AIC, 4), equals(102.4882))
  # All ones:
  tst <- matrix(1, 39, 5)
  res <- occSS.t(tst)
  expect_that(round(as.vector(res$real[,1]), 4), 
      is_equivalent_to(rep(1, 6)))
  expect_that(as.vector(res$real[, 2:3]), 
      is_equivalent_to(rep(NA_real_, 12)))
  expect_that(round(res$AIC, 4), equals(12))
  # All zeros:
  tst <- matrix(0, 39, 5)
  res <- occSS.t(tst)
  expect_that(as.vector(res$real), 
      is_equivalent_to(rep(NA_real_, 18)))
  expect_that(res$AIC, equals(NA_real_))
  # All NAs:
  tst <- matrix(NA, 39, 5)
  res <- occSS.t(tst)
  expect_that(as.vector(res$real), 
      is_equivalent_to(rep(NA_real_, 18)))
  expect_that(res$AIC,
      equals(NA_real_))
}  )
# ..............................................................

test_that("occSS.T gives right answers",  {
  # Data set (Blue Ridge Salamanders)
  require(wiqid)
  data(salamanders)
  BRS <- salamanders
  res <- occSS.T(BRS)

  expect_that(class(res), equals(c("occupancy", "list")))
  expect_that(names(res), equals(c("call", "beta", "real", "AIC")))
  expect_that(is.call(res$call), is_true())
  expect_that(colnames(res$real), equals(c("est", "lowCI", "uppCI")))
  expect_that(rownames(res$real),
       equals(c("psi", "p1",  "p2",  "p3",  "p4",  "p5")))
  # Values returned by PRESENCE:
  expect_that(round(as.vector(t(res$real)), 4), 
      equals(c( 0.5899, 0.3505, 0.7931,
                          0.1865, 0.0881, 0.3523,
                          0.2197, 0.1251, 0.3566,
                          0.2569, 0.1604, 0.3849,
                          0.2981, 0.1811, 0.4493,
                          0.3428, 0.1860, 0.5436)))
  expect_that(round(res$AIC, 4), equals(165.9228))

  # Put in some NAs
  BRS[c(6,167,130,123,89,154,32,120,127,147)] <- NA
  res <- occSS.T(BRS)
  expect_that(round(as.vector(res$real), 4), 
      is_equivalent_to(c(0.5800, 0.1939, 0.2189, 0.2461, 0.2754, 0.3068, 0.3234,
                         0.0879, 0.1183, 0.1432, 0.1531, 0.1493, 0.7997, 0.3752,
                         0.3692, 0.3892, 0.4442, 0.5276)))
  expect_that(round(res$AIC, 4), equals(150.8355))
  # Put in a row of NAs
  BRS[3,] <- NA 
  res <- occSS.T(BRS)
  expect_that(round(as.vector(res$real), 4), 
      equals(c(0.5494, 0.1913, 0.2211, 0.2542, 0.2904, 0.3294, 0.3087, 0.0847,
               0.1184, 0.1481, 0.1623, 0.1615, 0.7689, 0.3767, 0.3751, 0.4005,
               0.4635, 0.5562)))
  expect_that(round(res$AIC, 4), equals(144.8956))
  # Put in a column of NAs
  BRS[, 3] <- NA
  res <- occSS.T(BRS)
  expect_that(round(as.vector(res$real), 4), 
      equals(c( 0.3666, 0.2247, 0.2722, 0.3255, 0.3837, 0.4455, 0.1889, 0.0877,
                0.1307, 0.1699, 0.1919, 0.1965, 0.5899, 0.4662, 0.4820, 0.5321,
                0.6202, 0.7252)))
  expect_that(round(res$AIC, 4), equals(101.3259))
  # All ones:
  tst <- matrix(1, 39, 5)
  res <- occSS.T(tst)
  expect_that(round(as.vector(res$real[1,1]), 4), 
      equals(1))
  expect_that(as.vector(res$real)[-1], 
      is_equivalent_to(rep(NA_real_, 17)))
  expect_that(round(res$AIC, 4), equals(NA_real_))
  # All zeros:
  tst <- matrix(0, 39, 5)
  res <- occSS.T(tst)
  expect_that(as.vector(res$real), 
      equals(rep(NA_real_, 18)))
  expect_that(res$AIC, equals(NA_real_))
  # All NAs:
  tst <- matrix(NA, 39, 5)
  res <- occSS.T(tst)
  expect_that(as.vector(res$real), 
      equals(rep(NA_real_, 18)))
  expect_that(res$AIC, equals(NA_real_))
}  )
# .............................................

test_that("occSSnr gives right answers",  {
  # Data set (Blue Ridge Salamanders)
  require(wiqid)
  data(salamanders)
  BRS <- salamanders
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

