context("testing window data related functions")

test_that("test formatPositionDataToWindow", {
  
  positions <- data.frame(chr=c("2","2","1","1"),
                          position=c(2,1,2,1),
                          stringsAsFactors = F)
  expect_windowdata <- data.frame(row.names = c("2_1_500000","1_1_500000"),
                                  counts = c(2,2),
                                  stringsAsFactors = F)
  
  windowdata <- formatPositionDataToWindow(positions)
  
  expect_equal(windowdata,expect_windowdata)
  
})

test_that("test formatBedDataToWindow", {

  bed_table <- data.frame(chr=c("2","2","1","1"),
                          start=c(1,10,40,60),
                          end=c(5,20,50,120),
                          stringsAsFactors = F)
  expect_windowdata <- data.frame(row.names = c("1_1_500000","2_1_500000"),
                                  signal = c(0.000144,0.000032),
                                  stringsAsFactors = F)

  windowdata <- formatBedDataToWindow(bed_table)
  windowdata <- windowdata[windowdata$signal>0,,drop=F]

  expect_equal(windowdata,expect_windowdata)

})

test_that("test formatBedDataToWindow with signal", {
  
  bed_table <- data.frame(chr=c("2","2","1","1"),
                          start=c(1,10,40,60),
                          end=c(5,20,50,120),
                          signal = c(1,2,3,4),
                          stringsAsFactors = F)
  expect_windowdata <- data.frame(row.names = c("1_1_500000","2_1_500000"),
                                  signal = c(0.000554,0.000054),
                                  stringsAsFactors = F)
  
  windowdata <- formatBedDataToWindow(bed_table)
  windowdata <- windowdata[windowdata$signal>0,,drop=F]
  
  expect_equal(windowdata,expect_windowdata)
  
})

test_that("test formatBedpeBreakpointsDataToWindow", {
  
  sv_bedpe <- data.frame(chrom1=c(1,1,2),
                         start1=c(1000,10000,20000),
                         end1=c(1001,10001,20001),
                         chrom2=c(1,2,2),
                         start2=c(2000,11000,30000),
                         end2=c(2001,11001,30001),
                         svclass=c("tandem-duplication","translocation","deletion"),
                         stringsAsFactors = F)
  expect_windowdata <- data.frame(row.names = c("1_1_500000","2_1_500000"),
                                  counts = c(3,3),
                                  stringsAsFactors = F)
  
  res_windowdata <- formatBedpeBreakpointsDataToWindow(sv_bedpe = sv_bedpe,
                                                       genomev = "hg19")
  
  expect_equal(res_windowdata,expect_windowdata)
  
})

test_that("test sortWindowData", {
  
  windowdata <- data.frame(row.names = c("2_1_500000","1_500001_1000000","1_1_500000"),
                           counts = c(2,3,7),
                           stringsAsFactors = F)
  expect_windowdata_sorted <- windowdata[c(3,2,1),,drop=F]
  
  res_windowdata <- sortWindowData(windowData = windowdata)
  
  expect_equal(res_windowdata,expect_windowdata_sorted)
})

test_that("test mergeWindowData", {
  
  windowdata1 <- data.frame(row.names = c("2_1_500000","1_500001_1000000","1_1_500000"),
                           counts = c(2,3,7),
                           stringsAsFactors = F)
  windowdata2 <- data.frame(row.names = c("2_1_500000","1_500001_1000000","1_1_500000","2_500001_1000000"),
                            signal = c(2,3,7,8),
                            stringsAsFactors = F)
  
  expect_windowdata_merged <-  data.frame(row.names = c("1_1_500000","1_500001_1000000","2_1_500000","2_500001_1000000"),
                                          counts = c(7,3,2,NA),
                                          signal = c(7,3,2,8),
                                          stringsAsFactors = F)
  
  res_windowdata <- mergeWindowData(windowData1 = windowdata1,
                                    windowData2 = windowdata2)
  
  expect_equal(res_windowdata,expect_windowdata_merged)
})

test_that("test getBedFromWindowData", {
  
  windowdata <- data.frame(row.names = c("2_1_500000","1_500001_1000000","1_1_500000"),
                           counts = c(2,3,7),
                           stringsAsFactors = F)
  expect_windowdata_bed <- cbind(data.frame(row.names = rownames(windowdata),
                                            chr = as.character(c(2,1,1)),
                                            start = c(1,500001,1),
                                            end = c(500000,1000000,500000),
                                            stringsAsFactors = F),
                                 windowdata)
  
  res_windowdata <- getBedFromWindowData(windowData = windowdata)
  
  expect_equal(res_windowdata,expect_windowdata_bed)
})
