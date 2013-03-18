
# test_that code for rich series of functions

# library(testthat)
# test_file("test-rich.R")
# library(wiqid)

if(FALSE) {
context("Rarifaction")

test_that("richRarify gives the right answers",  {
  # Killarny bird data
  data(KillarneyBirds)
  kil <- KillarneyBirds[, -10]
  kilCntVec <- rowSums(kil)
  res <- richChao1(kilCntVec)
  expect_that(names(res), equals(c("Chao1", "Chao1Low", "Chao1Upp", "Chao1SD")))
 
} )
}

context("Chao1 and ACE")

test_that("richChao1 gives the right answers",  {
  # Killarny bird data
  data(KillarneyBirds)
  kil <- KillarneyBirds[, -10]
  kilCntVec <- rowSums(kil)
  res <- richChao1(kilCntVec)
  expect_that(names(res), equals(c("Chao1", "Chao1Low", "Chao1Upp", "Chao1SD")))
  expect_that(round(res, 2), 
      is_equivalent_to(c(33, 31.28, 45.36, 2.65))) # EstS output
  res <- richChao1(kil)
  expect_that(names(res), equals(c("Chao1", "Chao1Low", "Chao1Upp", "Chao1SD")))
  expect_that(round(res, 2), 
      is_equivalent_to(c(33, 31.28, 45.36, 2.65))) # EstS output

  # Bias-corrected version
  res <- richChao1(kilCntVec, correct=TRUE)
  expect_that(names(res), equals(c("Chao1", "Chao1Low", "Chao1Upp", "Chao1SD")))
  expect_that(round(res, 2), 
      is_equivalent_to(c(32.2, 31.14, 41.37, 1.84)))
        # EstS output is 32.2, 31.13, 41.36, 1.84, but is rounded down.
  # no singletons:
  res <- richChao1(kilCntVec+1)
  expect_that(round(res, 2), 
      is_equivalent_to(c(31, 31, 31, 0)))
        # EstS output is 31, 30.99, 31, 0, but is rounded down.
  # no singletons or doubletons:
  res <- richChao1(kilCntVec+2)
  expect_that(round(res, 2), 
      is_equivalent_to(c(31, 31, 31, 0)))
        # EstS output is 31, 30.99, 31, 0, but is rounded down.
  #TODO check for no doubletons when singletons present.
})


test_that("richACE gives the right answers",  {
  # Killarny bird data
  data(KillarneyBirds)
  kil <- KillarneyBirds[, -10]
  kilCntVec <- rowSums(kil)
  res <- richACE(kilCntVec)
  expect_that(length(res), equals(1))
  expect_that(round(res, 2), 
      is_equivalent_to(33.03)) # EstS output is 33.02
  res <- richACE(kil)
  expect_that(round(res, 2), 
      is_equivalent_to(33.03)) # EstS output is 33.02
  # no singletons:
  res <- richACE(kilCntVec+1)
  expect_that(round(res, 2), 
      is_equivalent_to(31)) # EstS output
  # no singletons or doubletons:
  res <- richACE(kilCntVec+2)
  expect_that(round(res, 2), 
      is_equivalent_to(31)) # EstS output
} )


context("Chao2 and ICE")

test_that("richChao2 gives the right answers",  {
  # Killarny bird data
  data(KillarneyBirds)
  kil <- KillarneyBirds[, -10]
  kilInc <- (kil > 0) * 1
  res <- richChao2(kilInc)
  expect_that(names(res), equals(c("Chao2", "Chao2Low", "Chao2Upp", "Chao2SD")))
  expect_that(round(res, 2), 
      is_equivalent_to(c(34.6, 31.65, 50.94, 3.85))) # EstS output
  res <- richChao2(kil)
  expect_that(round(res, 2), 
      is_equivalent_to(c(34.6, 31.65, 50.94, 3.85))) # EstS output

  # Bias-corrected version
  res <- richChao2(kilInc, correct=TRUE)
  expect_that(round(res, 2), 
      is_equivalent_to(c(33.22, 31.36, 44.88, 2.63)))
        # EstS output is 33.22, 31.35, 44.88, 2.62, but is rounded down.
  
  #TODO check for no uniques, no duplicates.
} )


test_that("richICE gives the right answers",  {
  # Killarny bird data
  data(KillarneyBirds)
  kil <- KillarneyBirds[, -10]
  kilInc <- (kil > 0) * 1
  res <- richICE(kilInc)
  expect_that(length(res), equals(1))
  expect_that(round(res, 2), 
      is_equivalent_to(34.86)) # EstS output
  res <- richICE(kil)
  expect_that(round(res, 2), 
      is_equivalent_to(34.86)) # EstS output

} )





