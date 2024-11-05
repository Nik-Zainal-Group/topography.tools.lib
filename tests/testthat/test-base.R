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

test_that("test breakDownOverlappingBedRegions", {
  
  bed_table <- data.frame(chr=c(1,1),
                          start=c(5000,10000),
                          end=c(12000,15000),
                          signal=c(1,1),
                          text=c(10,20),
                          id=c(1,2),
                          stringsAsFactors = F)
  bed_table_expected <- data.frame(row.names = c("1_5000_9999","1_10000_12000","1_12001_15000"),
                                   id=c("1","1;2","2"),
                                   chr=c("1","1","1"),
                                   start=c(5000,10000,12001),
                                   end=c(9999,12000,15000),
                                   signal=c(1,1,1),
                                   text=c("10","10;20","20"),
                                   stringsAsFactors = F)
  
  bed_table_res <- breakDownOverlappingBedRegions(bed_table = bed_table,aggregateTextColumns = "text")
  
  expect_equal(bed_table_res,bed_table_expected)
  
})

test_that("test assignBedRegionsToNonOverlappingSets", {
  
  bed_table <- data.frame(chr=c(1,1),
                          start=c(5000,10000),
                          end=c(12000,15000),
                          signal=c(1,1),
                          text=c(10,20),
                          id=c(1,2),
                          stringsAsFactors = F)
  sets_expected <- list(`1`=1,`2`=2)
  
  sets_res <- assignBedRegionsToNonOverlappingSets(bed_table = bed_table)
  
  expect_equal(sets_res,sets_expected)
  
})

test_that("test plotBrokenDownBedRegions", {
  
  bed_table <- data.frame(chr=c(1,1),
                          start=c(1,5),
                          end=c(10,15),
                          signal=c(1,1),
                          text=c(10,20),
                          id=c(1,2),
                          stringsAsFactors = F)
  
  filename <- "test-plotBrokenDownBedRegions.pdf"
  
  sets_res <- plotBrokenDownBedRegions(bed_table = bed_table,
                                       filename = filename,
                                       chrom = 1,
                                       pstart = 0,
                                       pend = 20)
  
  expect_true(file.exists(filename))
  
  unlink(filename)
  
})

test_that("test getIMD", {
  
  positions <- data.frame(chr=c("2","2","1","1","1"),
                          position=c(2,1,2,1,5),
                          stringsAsFactors = F)
  expect_IMD <- cbind(positions[c(4,3,5,2,1),],
                      data.frame(leftIMD=c(NA,1,3,NA,1),
                                 rightIMD=c(1,3,NA,1,NA),
                                 aveIMD=c(1,2,3,1,1),
                                 stringsAsFactors = F))
  
  positions_IMD <- getIMD(positions)
  
  expect_equal(positions_IMD,expect_IMD)
  
})
