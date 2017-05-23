# File created 2017-05-23

# test_that code for predict methods

# predict is implemented for occSS, occSS0, occSSrn, occSSrn0

# not implemented for occSScovsite, surv*, 

context("Predict method")

test_that("Predict works for occSS",  {
  # Data set (weta)
  require(wiqid)
  data(weta)
  DH <- weta[, 1:5]
  weta.covs <- weta[, 6:11]
  set.seed(123)
  weta.covs$nons <- rnorm(nrow(weta))

  newdata <- expand.grid(Browsed=c(TRUE, FALSE), nons=c(-1,0,1))
  wetaBn <- occSS(DH, psi ~ Browsed + nons, data=weta.covs)
  # On logit scale:
  predBn <- predict(wetaBn, newdata, "psi")
  expect_equal(attr(predBn, "link"), "logit")
  expect_equal(nrow(predBn), nrow(newdata))
  expect_equal(rownames(predBn), rownames(newdata))
  expect_equal(colnames(predBn), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predBn), 4), c(0.5529, 0.6154, -0.6532, 1.7590))
  # With different CI
  predBn80 <- predict(wetaBn, newdata, "psi", ci=0.8)
  expect_equivalent(round(colMeans(predBn80), 4), c(0.5529, 0.6154, -0.2357, 1.3415))
  # On probability scale
  predBn <- predict(wetaBn, newdata, "psi", type="response")
  expect_equal(nrow(predBn), nrow(newdata))
  expect_equal(colnames(predBn), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predBn), 4), c(0.6228, 0.1278, 0.3491, 0.8217))
  
  wetaB <- occSS(DH, psi ~ Browsed, data=weta.covs)
  predB <- predict(wetaB, newdata, "psi", type="response")
  expect_equivalent(round(colMeans(predB), 4), c(0.6202, 0.1139, 0.3752, 0.8016))
  expect_equal(predB[1, ], predB[3,])

  weta0 <- occSS(DH)  # This should call occSS0
  expect_message(pred0 <- predict(weta0, newdata, "psi", type="response"),
    "This is an intercept-only model")
  expect_equal(nrow(pred0), nrow(newdata))
  expect_equivalent(round(colMeans(pred0), 4), c(0.6166, 0.0885, 0.4356, 0.7702 ))
  expect_equal(pred0[1, ], pred0[2,])
  
  newdata <- data.frame(ObsD = c("A", "B", "C"))
  wetaBnO <- occSS(DH, list(psi ~ Browsed + nons, p ~ ObsD), data=weta.covs)
  predBnO <- predict(wetaBnO, newdata, "p", type="response")
  expect_equal(nrow(predBnO), nrow(newdata))
  expect_equal(colnames(predBnO), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predBnO), 4), c(0.3488, 0.0744, 0.2198, 0.5046))

  newdata1 <- data.frame(ObsD = c("B", "C"))
  predBnO1 <- predict(wetaBnO, newdata1, "p", type="response")
  expect_equal(nrow(predBnO1), nrow(newdata1))
  expect_equal(colnames(predBnO1), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predBnO1), 4), c(0.4120, 0.0803, 0.2679, 0.5732))
  
  # Check error messages
  newdata <- data.frame(Browsed=c(TRUE, FALSE))
  expect_error(predict(wetaBn, newdata, "psi"),
    "Needed variable")
  expect_error(predict(wetaBn, newdata, "whatever"),
    "No submodel found for parameter")
  newdata2 <- data.frame(ObsD = c("B", "C", "X"))
  expect_error(predict(wetaBnO, newdata2, "p"),
    "factor ObsD has new levels X")

  # With probit link
  newdata <- expand.grid(Browsed=c(TRUE, FALSE), nons=c(-1,0,1))
  wetaBn <- occSS(DH, psi ~ Browsed + nons, data=weta.covs, link="probit")
  # On probit scale:
  predBn <- predict(wetaBn, newdata, "psi")
  expect_equal(attr(predBn, "link"), "probit")
  expect_equal(nrow(predBn), nrow(newdata))
  expect_equal(colnames(predBn), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predBn), 4), c(0.3375,  0.3712, -0.3900,  1.0650))
  # On probability scale
  predBn <- predict(wetaBn, newdata, "psi", type="response")
  expect_equal(nrow(predBn), nrow(newdata))
  expect_equal(colnames(predBn), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predBn), 4), c(0.6227, 0.1281, 0.3543, 0.8274))
  
})
# ..............................................................

test_that("predict works with occSSrn0",  {
  # Data set (Blue Ridge Salamanders)
  require(wiqid)
  data(salamanders)
  n <- rep(ncol(salamanders), nrow(salamanders))
  y <- rowSums(salamanders)
  
  res <- occSSrn0(y, n)
  newdata <- data.frame(dummy = 1:3)
  rownames(newdata) <- c("A", "B", "C")
  expect_message(predLam <- predict(res, newdata, "lambda"),
    "This is an intercept-only model")
  expect_equal(attr(predLam, "link"), "log")
  expect_equal(nrow(predLam), nrow(newdata))
  expect_equal(rownames(predLam), rownames(newdata))
  expect_equal(colnames(predLam), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predLam), 4), c(0.1332, 0.4494, -0.7476 ,1.0141))

  expect_message(predr <- predict(res, newdata, "r"),
    "This is an intercept-only model")
  expect_equal(attr(predr, "link"), "logit")
  expect_equal(nrow(predr), nrow(newdata))
  expect_equal(rownames(predr), rownames(newdata))
  expect_equal(colnames(predr), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predr), 4), c(-1.7546, 0.5261, -2.7857, -0.7235))

  predLam <- predict(res, newdata, "lambda", type="response")
  expect_equivalent(round(colMeans(predLam), 4), c(1.1425, 0.5135, 0.4735, 2.7568))
  predr <- predict(res, newdata, "r", type="response")
  expect_equivalent(round(colMeans(predr), 4), c(0.1475, 0.0661, 0.0581, 0.3266))

  res2 <- occSSrn0(y, n, link="probit")
  predLam <- predict(res2, newdata, "lambda")
  expect_equal(attr(predLam, "link"), "log")
  expect_equivalent(round(colMeans(predLam), 4), c(0.1332,  0.4495, -0.7477,  1.0142 ))
  predr <- predict(res2, newdata, "r")
  expect_equal(attr(predr, "link"), "probit")
  expect_equivalent(round(colMeans(predr), 4), c(-1.0474,  0.2870, -1.6098, -0.4849))

  predLam <- predict(res2, newdata, "lambda", type="response")
  expect_equivalent(round(colMeans(predLam), 4), c(1.1425, 0.5135, 0.4734, 2.7571))
  predr <- predict(res2, newdata, "r", type="response")
  expect_equivalent(round(colMeans(predr), 4), c(0.1475, 0.0661, 0.0537, 0.3139))
 } )
 
# ..............................................................

test_that("predict works with occSSrnSite",  {
  require(wiqid)
  data(weta)
  DH <- weta[, 1:5]
  y <- rowSums(DH, na.rm=TRUE)
  n <- rowSums(!is.na(DH))
  weta.covs <- weta[, 6:11]
  set.seed(123)
  weta.covs$nons <- rnorm(nrow(weta))

  newdata <- expand.grid(Browsed=c(TRUE, FALSE), nons=c(-1,0,1))
  rownames(newdata) <- LETTERS[1:nrow(newdata)]
  # Full model for lambda
  res <- occSSrnSite(y, n, lambda ~ Browsed + nons, data=weta.covs)
  expect_silent(predLam <- predict(res, newdata, "lambda"))
  expect_equal(attr(predLam, "link"), "log")
  expect_equal(nrow(predLam), nrow(newdata))
  expect_equal(rownames(predLam), rownames(newdata))
  expect_equal(colnames(predLam), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predLam), 4), c(0.1830,  0.4186, -0.6375,  1.0035))
  predLam <- predict(res, newdata, "lambda", type="response")
  expect_equivalent(round(colMeans(predLam), 4), c(1.2527, 0.5216, 0.5544, 2.8358))

  expect_message(predr <- predict(res, newdata, "r"),
    "This is an intercept-only model")
  expect_equal(attr(predr, "link"), "logit")
  expect_equal(nrow(predr), nrow(newdata))
  expect_equal(rownames(predr), rownames(newdata))
  expect_equal(colnames(predr), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predr), 4), c(-1.4015,  0.4381, -2.2602, -0.5428))
  predr <- predict(res, newdata, "r", type="response")
  expect_equivalent(round(colMeans(predr), 4), c(0.1976, 0.0695, 0.0945, 0.3675))

  res2 <- occSSrnSite(y, n, lambda ~ Browsed + nons, data=weta.covs, link="probit")
  predLam <- predict(res2, newdata, "lambda")
  expect_equal(attr(predLam, "link"), "log")
  expect_equivalent(round(colMeans(predLam), 4), c(0.1830,  0.4187, -0.6376,  1.0035))
  predLam <- predict(res2, newdata, "lambda", type="response")
  expect_equivalent(round(colMeans(predLam), 4), c(1.2527, 0.5217, 0.5543, 2.8360))

  predr <- predict(res2, newdata, "r")
  expect_equal(attr(predr, "link"), "probit")
  expect_equivalent(round(colMeans(predr), 4), c(-0.8503,  0.2500, -1.3402, -0.3603))
  predr <- predict(res2, newdata, "r", type="response")
  expect_equivalent(round(colMeans(predr), 4), c(0.1976, 0.0695, 0.0901, 0.3593))

  # Model for lambda and r
  res <- occSSrnSite(y, n, list(lambda ~ nons, r ~ Browsed), data=weta.covs)
  expect_silent(predLam <- predict(res, newdata, "lambda"))
  expect_equal(attr(predLam, "link"), "log")
  expect_equal(nrow(predLam), nrow(newdata))
  expect_equal(rownames(predLam), rownames(newdata))
  expect_equal(colnames(predLam), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predLam), 4), c(0.3014,  0.4099, -0.5020,  1.1047))
  predLam <- predict(res, newdata, "lambda", type="response")
  expect_equivalent(round(colMeans(predLam), 4), c(1.3535, 0.5552, 0.6062, 3.0266))

  expect_silent(predr <- predict(res, newdata, "r"))
  expect_equal(attr(predr, "link"), "logit")
  expect_equal(nrow(predr), nrow(newdata))
  expect_equal(rownames(predr), rownames(newdata))
  expect_equal(colnames(predr), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predr), 4), c(-1.5367,  0.5044, -2.5253, -0.5480))
  predr <- predict(res, newdata, "r", type="response")
  expect_equivalent(round(colMeans(predr), 4), c(0.1804, 0.0729, 0.0775, 0.3675 ))

  res2 <- occSSrnSite(y, n, list(lambda ~ nons, r ~ Browsed), data=weta.covs, link="probit")
  predLam <- predict(res2, newdata, "lambda")
  expect_equal(attr(predLam, "link"), "log")
  expect_equivalent(round(colMeans(predLam), 4), c(0.3014,  0.4100, -0.5021,  1.1049 ))
  predLam <- predict(res2, newdata, "lambda", type="response")
  expect_equivalent(round(colMeans(predLam), 4), c(1.3535, 0.5553, 0.6061, 3.0271))

  predr <- predict(res2, newdata, "r")
  expect_equal(attr(predr, "link"), "probit")
  expect_equivalent(round(colMeans(predr), 4), c(-0.9242,  0.2821, -1.4771, -0.3713 ))
  predr <- predict(res2, newdata, "r", type="response")
  expect_equivalent(round(colMeans(predr), 4), c(0.1804, 0.0729, 0.0728, 0.3563))

  # Intercept only:
  newdata <- expand.grid(Browsed=c(TRUE, FALSE), nons=c(-1,0,1))
  res <- occSSrnSite(y, n)
  rownames(newdata) <- LETTERS[1:nrow(newdata)]
  expect_message(predLam <- predict(res, newdata, "lambda"),
    "This is an intercept-only model")
  expect_equal(attr(predLam, "link"), "log")
  expect_equal(nrow(predLam), nrow(newdata))
  expect_equal(rownames(predLam), rownames(newdata))
  expect_equal(colnames(predLam), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predLam), 4), c(0.1506,  0.3329, -0.5018,  0.8031))
  predLam <- predict(res, newdata, "lambda", type="response")
  expect_equivalent(round(colMeans(predLam), 4), c(1.1626, 0.3870, 0.6055, 2.2324))

  expect_message(predr <- predict(res, newdata, "r"),
    "This is an intercept-only model")
  expect_equal(attr(predr, "link"), "logit")
  expect_equal(nrow(predr), nrow(newdata))
  expect_equal(rownames(predr), rownames(newdata))
  expect_equal(colnames(predr), c("est", "SE", "lowCI", "uppCI"))
  expect_equivalent(round(colMeans(predr), 4), c(-1.3296,  0.4147, -2.1423, -0.5168))
  predr <- predict(res, newdata, "r", type="response")
  expect_equivalent(round(colMeans(predr), 4), c(0.2092, 0.0686, 0.1050, 0.3736))

  res2 <- occSSrn0(y, n, link="probit")
  predLam <- predict(res2, newdata, "lambda")
  expect_equal(attr(predLam, "link"), "log")
  expect_equivalent(round(colMeans(predLam), 4), c(0.1506,  0.3329, -0.5019,  0.8031))
  predLam <- predict(res2, newdata, "lambda", type="response")
  expect_equivalent(round(colMeans(predLam), 4), c(1.1626, 0.3870, 0.6054, 2.2326))

  predr <- predict(res2, newdata, "r")
  expect_equal(attr(predr, "link"), "probit")
  expect_equivalent(round(colMeans(predr), 4), c(-0.8091,  0.2386, -1.2768, -0.3414))
  predr <- predict(res2, newdata, "r", type="response")
  expect_equivalent(round(colMeans(predr), 4), c(0.2092, 0.0686, 0.1008, 0.3664))
} )
 
 