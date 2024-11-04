context("testing base functions")

test_that("test sort chromosomes", {
  
  chroms <- c("12","7","Y","X")
  expect_sorted <- c("7","12","X","Y")
  
  chroms_sorted <- sortChroms(chroms)
  
  expect_equal( chroms_sorted, expect_sorted )
  
})