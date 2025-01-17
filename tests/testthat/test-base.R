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

test_that("test mergeAdjacentBedRegions", {
  
  bed_table <- data.frame(chr=c(1,1,2),
                          start=c(1,11,10),
                          end=c(10,15,15),
                          id=c("A","B","C"),
                          stringsAsFactors = F)
  expect_merge <- data.frame(id=c("A;B","C"),
                             chr=c("1","2"),
                             start=c(1,10),
                             end=c(15,15),
                             stringsAsFactors = F)
  
  bed_table_merged <- mergeAdjacentBedRegions(bed_table)
  
  rownames(expect_merge) <- c(1,nrow(expect_merge))
  rownames(bed_table_merged) <- c(1,nrow(bed_table_merged))
  
  expect_equal(bed_table_merged,expect_merge)
  
})

test_that("test mergeAdjacentBedRegions with signal", {
  
  bed_table <- data.frame(chr=c(1,1,2),
                          start=c(1,11,10),
                          end=c(10,15,15),
                          id=c("A","B","C"),
                          signal=c(2,1,1),
                          stringsAsFactors = F)
  expect_merge <- data.frame(id=c("A;B","C"),
                             chr=c("1","2"),
                             start=c(1,10),
                             end=c(15,15),
                             signal=c(5/3,1),
                             stringsAsFactors = F)
  
  bed_table_merged <- mergeAdjacentBedRegions(bed_table,
                                              aggregateSignalMode="weightedmean")
  
  rownames(expect_merge) <- c(1,nrow(expect_merge))
  rownames(bed_table_merged) <- c(1,nrow(bed_table_merged))
  
  expect_equal(bed_table_merged,expect_merge)
  
})

test_that("test getChromosomesBedTable", {
  
  chrTable <- getChromosomesBedTable(genomev = "hg19")
  
  expect_equal(nrow(chrTable),24)
  
})

test_that("test trimNfromBed", {
  
  bed_table <- data.frame(chr=c(1),
                          start=c(1),
                          end=c(300000),
                          stringsAsFactors = F)
  expected_trimmed <- data.frame(chr=c(1,1),
                                 start=c(10001,227418),
                                 end=c(177417,267719),
                                 stringsAsFactors = F)
  
  bed_table_trimmed <- trimNfromBed(bed_table = bed_table,
                                    genomev = "hg19")
  

  
  expect_equal(bed_table_trimmed,expected_trimmed)
  
})

test_that("test plotBedSignalRegion single signal", {
  
  bed_table <- data.frame(chr=c(1,1),
                          start=c(10000,50000),
                          end=c(100000,150000),
                          signal=c(1,2),
                          text=c(10,20),
                          id=c(1,2),
                          stringsAsFactors = F)
  
  filename <- "test-plotBedSignalRegion.pdf"
  
  sets_res <- plotBedSignalRegion(bed_table = bed_table,
                                  fileout = filename,
                                  pchr = 1,
                                  pstart = 1,
                                  pend = 200000)
  
  expect_true(file.exists(filename))
  
  unlink(filename)
  
})

test_that("test plotBedSignalRegion with two signals", {
  
  bed_table <- data.frame(chr=c(1,1),
                          start=c(10000,50000),
                          end=c(100000,150000),
                          signal=c(1,2),
                          text=c(10,20),
                          id=c(1,2),
                          stringsAsFactors = F)
  bed_table2 <- data.frame(chr=c(1),
                          start=c(90000),
                          end=c(220000),
                          signal=c(5),
                          id=c(1),
                          stringsAsFactors = F)
  
  filename <- "test-plotBedSignalRegion2.pdf"
  
  sets_res <- plotBedSignalRegion(bed_table = bed_table,
                                  bed_table2 = bed_table2,
                                  fileout = filename,
                                  pchr = 1,
                                  pstart = 1,
                                  pend = 200000)
  
  expect_true(file.exists(filename))
  
  unlink(filename)
  
})

test_that("test extendBedRegionsWithNoOverlap", {
  
  bed_table <- data.frame(chr=c(1,1),
                          start=c(3000,15000),
                          end=c(5000,20000),
                          id=c(1,2),
                          stringsAsFactors = F)
  
  extended_expected <- data.frame(chr=c(1,1),
                                  start=c(1,10001),
                                  end=c(10000,25000),
                                  id=c(1,2),
                                  stringsAsFactors = F)
  
  bed_table_extended <- extendBedRegionsWithNoOverlap(bed_table = bed_table,
                                                      extended = 5000)
  
  expect_equal(bed_table_extended,extended_expected)

})

test_that("test getIRD", {
  
  bed_table <- data.frame(chr=c(1,1,1,2),
                          start=c(3000,15000,20001,10000),
                          end=c(5000,20000,25000,20000),
                          id=c(1,2,3,4),
                          stringsAsFactors = F)
  
  res_ird <- getIRD(bed_table = bed_table)
  
  expect_equal(res_ird$aveIRD,c(10000,5000.5,1,NA))
  
})

test_that("test distanceOfPositionToNearestBedRegion", {
  
  positions <- data.frame(chr = c(1,1,1,3,2,2,2,2),
                          position = c(100,200,300,400,500,600,700,800),
                          id=paste0("y",c(1,2,3,4,5,6,7,8)),
                          stringsAsFactors = F)
  bed_table <- data.frame(chr = c(1,1,2,2),
                          start = c(95,350,450,506),
                          end = c(205,405,505,705),
                          id=paste0("x",c(1,2,3,4)),
                          stringsAsFactors = F)
  
  res_dist <- distanceOfPositionToNearestBedRegion(positions = positions,
                                                   bed_table = bed_table)
  
  expect_equal(res_dist$distanceToNearestRegion,c(0,0,-50,0,0,0,95,NA))

})

test_that("test extendBedRegions", {
  
  bed_table <- data.frame(chr=c(1,1,1,2),
                          start=c(3000,15000,20001,10000),
                          end=c(5000,20000,25000,20000),
                          id=c(1,2,3,4),
                          stringsAsFactors = F)
  
  res_ext <- extendBedRegions(bed_table = bed_table,
                              extended = 4000)
  
  expect_equal(c(res_ext$start,res_ext$end),c(1,11000,16001, 6000,9000,24000,29000,24000))
  
})

test_that("test bedpeBreakpointsToPositions", {
  
  sv_bedpe <- data.frame(chrom1=c(1,1,2),
                         start1=c(1000,10000,20000),
                         end1=c(1001,10001,20001),
                         chrom2=c(1,2,2),
                         start2=c(2000,11000,30000),
                         end2=c(2001,11001,30001),
                         svclass=c("tandem-duplication","translocation","deletion"),
                         stringsAsFactors = F)
  
  positions_expected <- data.frame(chr=c("1","1","1","2","2","2"),
                                   position=c(1000,2000,10000,11000,20000,30000),
                                   svclass=c("tandem-duplication","tandem-duplication",
                                             "translocation","translocation",
                                             "deletion","deletion"),
                                   stringsAsFactors = F)
  
  positions <- bedpeBreakpointsToPositions(sv_bedpe = sv_bedpe,
                                           copycolumns = "svclass")
  rownames(positions) <- 1:nrow(positions)
  
  expect_equal(positions,positions_expected)
  
})
