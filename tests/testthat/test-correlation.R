context("testing correlation functions")


test_that("test correlatePositionsWithBedRegions with classes", {
  positions <- data.frame(chr = c(1,1,1,1,1,1,1,1),
                          position = c(100,200,300,400,500,600,700,800),
                          id=paste0("y",c(1,2,3,4,5,6,7,8)),
                          class=c("C","B","C","A","A","B","B","C"),
                          stringsAsFactors = F)
  bed_table <- data.frame(chr = c(1,1,1,1),
                          start = c(95,206,350,506),
                          end = c(205,405,505,705),
                          id=paste0("x",c(1,2,3,4)),
                          class=c("g2","g1","g1","g2"),
                          stringsAsFactors = F)
  samplingRegions <- data.frame(chr=c(1,1),
                                start=c(1,5001),
                                end=c(5000,10000),
                                stringsAsFactors = F)
  
  expected_summaryOverlaps <- data.frame(row.names = c("positions","regions"),
                                  noverlap=c(7,4),
                                  ntotal=c(8,4),
                                  stringsAsFactors = F)
  
  res_obj <- correlatePositionsWithBedRegions(positions = positions,
                                              bed_table = bed_table,
                                              samplingRegions = samplingRegions,
                                              genomev = NULL,
                                              nsamples = 100)
  
  expect_equal(res_obj$summaryOverlaps,expected_summaryOverlaps)
  expect_true(all(res_obj$pvalue_positionsInRegionClasses[,"noMatch"]>0.5) & res_obj$pvalue_positionsInRegionClasses["B","g2"]<0.5 & res_obj$pvalue_positionsInRegionClasses["A","g1"]<0.5)
  
})


test_that("test correlatePositionsWithBedRegions with random seed", {
  positions <- data.frame(chr = c(1,1,1,1,1,1,1,1),
                          position = c(100,200,300,400,500,600,700,800),
                          id=paste0("y",c(1,2,3,4,5,6,7,8)),
                          class=c("C","B","C","A","A","B","B","C"),
                          stringsAsFactors = F)
  bed_table <- data.frame(chr = c(1,1,1,1),
                          start = c(95,206,350,506),
                          end = c(205,405,505,705),
                          id=paste0("x",c(1,2,3,4)),
                          class=c("g2","g1","g1","g2"),
                          stringsAsFactors = F)
  samplingRegions <- data.frame(chr=c(1,1),
                                start=c(1,5001),
                                end=c(5000,10000),
                                stringsAsFactors = F)
  
  res_obj <- correlatePositionsWithBedRegions(positions = positions,
                                              bed_table = bed_table,
                                              samplingRegions = samplingRegions,
                                              genomev = NULL,randomSeed = 1,
                                              nsamples = 10)
  res_obj2 <- correlatePositionsWithBedRegions(positions = positions,
                                               bed_table = bed_table,
                                               samplingRegions = samplingRegions,
                                               genomev = NULL,randomSeed = 1,
                                               nsamples = 10)
  
  expect_equal(res_obj$sampled_PostionsInAnyRegion,res_obj2$sampled_PostionsInAnyRegion)
  
})