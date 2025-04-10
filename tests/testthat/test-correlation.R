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

test_that("test correlatePositionsWithBedRegions with classes and reuse resampling", {
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
                                              returnResampledPositions = TRUE,
                                              returnResampledBedRegions = TRUE,
                                              genomev = NULL,
                                              nsamples = 10)
  res_obj2 <- correlatePositionsWithBedRegions(positions = positions,
                                               bed_table = bed_table,
                                               samplingRegions = samplingRegions,
                                               resampled_positions_list = res_obj$resampled_positions_list,
                                               resampled_bed_regions_list = res_obj$resampled_bed_regions_list,
                                               genomev = NULL,
                                               nsamples = 10)
  
  expect_equal(res_obj$sampled_PostionsInAnyRegion,res_obj2$sampled_PostionsInAnyRegion)

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

test_that("test correlateBedRegions without classes", {
  bed_table1 <- data.frame(chr = c(1,1,1),
                           start = c(100,220,300),
                           end = c(150,350,420),
                           id=paste0("x",c(1,2,3)),
                           stringsAsFactors = F)
  bed_table2 <- data.frame(chr = c(1,1,1,1),
                           start = c(200,320,450,230),
                           end = c(500,400,550,250),
                           id=paste0("y",c(1,2,3,4)),
                           stringsAsFactors = F)
  samplingRegions <- data.frame(chr=c(1,1),
                                start=c(1,3001),
                                end=c(3000,6000),
                                stringsAsFactors = F)
  
  expected_summaryOverlaps <- data.frame(row.names = c("bed_table1","bed_table2"),
                                         noverlap=c(2,3),
                                         ntotal=c(3,4),
                                         stringsAsFactors = F)
  
  res_obj <- correlateBedRegions(bed_table1 = bed_table1,
                                 bed_table2 = bed_table2,
                                 samplingRegions = samplingRegions,
                                 genomev = NULL,
                                 nsamples = 100)
  
  expect_equal(res_obj$summaryOverlaps,expected_summaryOverlaps)
  expect_true(res_obj$pvalue_Regions1overlappingAnyRegion2<0.5 & res_obj$pvalue_Regions2overlappingAnyRegion1<0.5)
  
})

test_that("test correlateBedRegions with classes", {
  bed_table1 <- data.frame(chr = c(1,1,1),
                           start = c(100,220,300),
                           end = c(150,350,420),
                           id=paste0("x",c(1,2,3)),
                           class=c("A","B","B"),
                           stringsAsFactors = F)
  bed_table2 <- data.frame(chr = c(1,1,1,1),
                           start = c(200,320,450,230),
                           end = c(500,400,550,250),
                           id=paste0("y",c(1,2,3,4)),
                           class=c("M","M","N","N"),
                           stringsAsFactors = F)
  samplingRegions <- data.frame(chr=c(1,1),
                                start=c(1,5001),
                                end=c(5000,10000),
                                stringsAsFactors = F)
  
  expected_summaryOverlaps <- data.frame(row.names = c("bed_table1","bed_table2"),
                                         noverlap=c(2,3),
                                         ntotal=c(3,4),
                                         stringsAsFactors = F)
  
  res_obj <- correlateBedRegions(bed_table1 = bed_table1,
                                 bed_table2 = bed_table2,
                                 samplingRegions = samplingRegions,
                                 genomev = NULL,
                                 nsamples = 100)
  
  expect_equal(res_obj$summaryOverlaps,expected_summaryOverlaps)
  expect_true(all(res_obj$pvalue_regions1overlappingRegion2classes[,"noMatch"]>0.5) & res_obj$pvalue_regions1overlappingRegion2classes["B","M"]<0.5 & res_obj$pvalue_regions1overlappingRegion2classes["B","N"]<0.5)
  
})

test_that("test correlateBedRegions with classes and reuse resampling", {
  bed_table1 <- data.frame(chr = c(1,1,1),
                           start = c(100,220,300),
                           end = c(150,350,420),
                           id=paste0("x",c(1,2,3)),
                           class=c("A","B","B"),
                           stringsAsFactors = F)
  bed_table2 <- data.frame(chr = c(1,1,1,1),
                           start = c(200,320,450,230),
                           end = c(500,400,550,250),
                           id=paste0("y",c(1,2,3,4)),
                           class=c("M","M","N","N"),
                           stringsAsFactors = F)
  samplingRegions <- data.frame(chr=c(1,1),
                                start=c(1,5001),
                                end=c(5000,10000),
                                stringsAsFactors = F)
  
  res_obj <- correlateBedRegions(bed_table1 = bed_table1,
                                 bed_table2 = bed_table2,
                                 samplingRegions = samplingRegions,
                                 returnResampledBedRegions1 = TRUE,
                                 returnResampledBedRegions2 = TRUE,
                                 genomev = NULL,
                                 nsamples = 10)
  res_obj2 <- correlateBedRegions(bed_table1 = bed_table1,
                                  bed_table2 = bed_table2,
                                  samplingRegions = samplingRegions,
                                  resampled_bed_regions1_list = res_obj$resampled_bed_regions1_list,
                                  resampled_bed_regions2_list = res_obj$resampled_bed_regions2_list,
                                  genomev = NULL,
                                  nsamples = 10)
  
  expect_equal(res_obj$sampled_Regions1overlappingAnyRegion2,res_obj2$sampled_Regions1overlappingAnyRegion2)

})



test_that("test multipleCorrelations", {
  bed_table1 <- data.frame(chr = c(1,1,1),
                           start = c(100,220,300),
                           end = c(150,350,420),
                           id=paste0("x",c(1,2,3)),
                           class=c("A","B","B"),
                           stringsAsFactors = F)
  bed_table2 <- data.frame(chr = c(1,1,1,1),
                           start = c(200,320,450,230),
                           end = c(500,400,550,250),
                           id=paste0("y",c(1,2,3,4)),
                           class=c("M","M","N","N"),
                           stringsAsFactors = F)
  bed_table3 <- data.frame(chr = c(1,1,1),
                           start = c(200,1200,1100),
                           end = c(500,1600,1300),
                           id=paste0("z",c(1,2,3)),
                           class=c("E","E","F"),
                           stringsAsFactors = F)
  samplingRegions <- data.frame(chr=c(1,1),
                                start=c(1,2001),
                                end=c(2000,4000),
                                stringsAsFactors = F)
  
  expected_summaryOverlaps <- data.frame(row.names = "bed1",
                                         bed1=3,
                                         bed2=2,
                                         bed3=2,
                                         total=3,
                                         stringsAsFactors = F)
  
  res_obj <- multipleCorrelations(referenceEntities = bed_table1,
                                  compareEntitiesList = list(bed1=bed_table1,bed2=bed_table2,bed3=bed_table3),
                                  referenceEntitiesName = "bed1",
                                  resampleCompareEntities = FALSE,
                                  resampleReferenceEntitiesAllowOverlap = TRUE,
                                  samplingRegions = samplingRegions,
                                  genomev = NULL,
                                  nsamples = 100)
  
  expect_equal(res_obj$counts_refEntitiesWithCompEntities,expected_summaryOverlaps)
  expect_true(!is.null(res_obj$pvalues_refEntitiesWithCompEntities))
  
})



test_that("test multipleCorrelations with positions", {
  bed_table1 <- data.frame(chr = c(1,1,1),
                           start = c(100,220,400),
                           end = c(150,350,420),
                           id=paste0("x",c(1,2,3)),
                           class=c("A","B","B"),
                           stringsAsFactors = F)
  pos_table1 <- data.frame(chr = c(1,1,1),
                           position = c(240,320,410),
                           id=paste0("y",c(1,2,3)),
                           class=c("M","M","N"),
                           stringsAsFactors = F)
  pos_table2 <- data.frame(chr = c(1,1),
                           position = c(120,200),
                           id=paste0("y",c(1,2)),
                           class=c("M","M"),
                           stringsAsFactors = F)
  samplingRegions <- data.frame(chr=c(1,1),
                                start=c(1,2001),
                                end=c(2000,4000),
                                stringsAsFactors = F)
  
  expected_summaryOverlaps <- data.frame(row.names = "bed1",
                                         pos1=2,
                                         pos2=1,
                                         total=3,
                                         stringsAsFactors = F)
  
  res_obj <- multipleCorrelations(referenceEntities = bed_table1,
                                  compareEntitiesList = list(pos1=pos_table1,pos2=pos_table2),
                                  referenceEntitiesName = "bed1",
                                  resampleCompareEntities = FALSE,
                                  resampleReferenceEntitiesAllowOverlap = TRUE,
                                  samplingRegions = samplingRegions,
                                  genomev = NULL,
                                  nsamples = 100)
  
  expect_equal(res_obj$counts_refEntitiesWithCompEntities,expected_summaryOverlaps)
  expect_true(!is.null(res_obj$pvalues_refEntitiesWithCompEntities))
  
})
