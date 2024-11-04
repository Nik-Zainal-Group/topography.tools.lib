context("testing base functions")

test_that("test sortChroms", {
  
  chroms <- c("12","7","Y","X")
  expect_sorted <- c("7","12","X","Y")
  
  chroms_sorted <- sortChroms(chroms)
  
  expect_equal(chroms_sorted,expect_sorted)
  
})

test_that("test sortPositions", {
  
  positions <- data.frame(chr=c("2","2","1","1"),position=c(2,1,2,1))
  expect_sorted <- positions[c(4,3,2,1),]
  
  positions_sorted <- sortPositions(positions)
  
  expect_equal(positions_sorted,expect_sorted)
  
})

test_that("test sortBed", {
  
  bed_table <- data.frame(chr=c("2","2","1","1"),start=c(3,2,2,1),end=c(1,4,3,2))
  expect_sorted <- bed_table[c(4,3,1,2),]
  expect_sorted[3,c("start","end")] <- c(1,3)
  
  bed_table_sorted <- sortBed(bed_table)
  
  expect_equal(bed_table_sorted,expect_sorted)
  
})

test_that("test checkBedRegionsOverlap", {
  
  bed_table1 <- data.frame(chr=c("2","2","1","1"),start=c(3,2,2,1),end=c(1,4,3,2))
  bed_table2 <- data.frame(chr=c("2","2","1","1"),start=c(3,5,2,6),end=c(1,6,3,8))
  
  bed_table_res1 <- checkBedRegionsOverlap(bed_table1)
  bed_table_res2 <- checkBedRegionsOverlap(bed_table2)
  
  expect_equal(bed_table_res1,c("1","2"))
  expect_equal(bed_table_res2,NULL)
  
})
