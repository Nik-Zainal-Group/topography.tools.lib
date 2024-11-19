context("testing intersection functions")

test_that("test intersectionStats", {
  
  idMap1to2 <- list(x2=c("y1","y2","y4"),x3=c("y1","y2"))
  idMap2to1 <- list(y1=c("x2","x3"),y2=c("x2","x3"),y4="x2")
  idclassmap1 <- data.frame(id=c("x1","x2","x3"),
                            class=c("A","B","B"),
                            stringsAsFactors = F)
  idclassmap2 <- data.frame(id=c("y1","y2","y3","y4"),
                            class=c("X","X","Y","Y"),
                            stringsAsFactors = F)
  res_obj1 <- intersectionStats(idMap = idMap1to2,
                                idclassmap1 = idclassmap1,
                                idclassmap2 = idclassmap2)
  res_obj2 <- intersectionStats(idMap = idMap2to1,
                                idclassmap1 = idclassmap2,
                                idclassmap2 = idclassmap1)
  
  expected_classes1 <- data.frame(row.names = c("A","B"),
                                  X=c(0,2),
                                  Y=c(0,1),
                                  noMatch=c(1,0),
                                  stringsAsFactors = F)
  expected_classes2 <- data.frame(row.names = c("X","Y"),
                                  A=c(0,0),
                                  B=c(2,1),
                                  noMatch=c(0,1),
                                  stringsAsFactors = F)
  
  expect_equal(res_obj1$totalId1matchingAnyId2,2)
  expect_equal(res_obj2$totalId1matchingAnyId2,3)
  expect_equal(res_obj1$countsTable_classes1_in_classes2,expected_classes1)
  expect_equal(res_obj2$countsTable_classes1_in_classes2,expected_classes2)
  
  
})


test_that("test intersectPositionsAndBedRegions_nonOverlapping without classes", {
  positions <- data.frame(chr = c(1,1,1,1,1,1,1,1),
                          position = c(100,200,300,400,500,600,700,800),
                          id=paste0("y",c(1,2,3,4,5,6,7,8)),
                          stringsAsFactors = F)
  bed_table <- data.frame(chr = c(1,1,1,1),
                          start = c(95,206,450,506),
                          end = c(205,405,505,705),
                          id=paste0("x",c(1,2,3,4)),
                          stringsAsFactors = F)
  
  expected_classes1 <- data.frame(row.names = "anyPosition",
                                  anyRegion=7,
                                  noMatch=1,
                                  stringsAsFactors = F)
  expected_classes2 <- data.frame(row.names = "anyRegion",
                                  anyPosition=4,
                                  noMatch=0,
                                  stringsAsFactors = F)
  
  res_obj <- intersectPositionsAndBedRegions_nonOverlapping(positions = positions,
                                                            bed_table = bed_table)
  
  expect_equal(res_obj$totalPostionsInAnyRegion,7)
  expect_equal(res_obj$totalRegionsAtAnyPosition,4)
  expect_equal(res_obj$countsTable_positionsInRegionClasses,expected_classes1)
  expect_equal(res_obj$countsTable_regionsAtPositionClasses,expected_classes2)
  
})


test_that("test intersectPositionsAndBedRegions_nonOverlapping with classes", {
  positions <- data.frame(chr = c(1,1,1,1,1,1,1,1),
                          position = c(100,200,300,400,500,600,700,800),
                          id=paste0("y",c(1,2,3,4,5,6,7,8)),
                          class=c("C","B","C","A","A","B","B","C"),
                          stringsAsFactors = F)
  bed_table <- data.frame(chr = c(1,1,1,1),
                          start = c(95,206,450,506),
                          end = c(205,405,505,705),
                          id=paste0("x",c(1,2,3,4)),
                          class=c("g2","g1","g1","g2"),
                          stringsAsFactors = F)
  
  expected_classes1 <- data.frame(row.names = c("C","B","A"),
                                  g2=c(1,3,0),
                                  g1=c(1,0,2),
                                  noMatch=c(1,0,0),
                                  stringsAsFactors = F)
  expected_classes2 <- data.frame(row.names = c("g2","g1"),
                                  C=c(1,1),
                                  B=c(2,0),
                                  A=c(0,2),
                                  noMatch=c(0,0),
                                  stringsAsFactors = F)
  
  res_obj <- intersectPositionsAndBedRegions_nonOverlapping(positions = positions,
                                                            bed_table = bed_table)
  
  expect_equal(res_obj$totalPostionsInAnyRegion,7)
  expect_equal(res_obj$totalRegionsAtAnyPosition,4)
  expect_equal(res_obj$countsTable_positionsInRegionClasses,expected_classes1)
  expect_equal(res_obj$countsTable_regionsAtPositionClasses,expected_classes2)
  
})


test_that("test intersectPositionsAndBedRegions without classes", {
  positions <- data.frame(chr = c(1,1,1,1,1,1,1,1),
                          position = c(100,200,300,400,500,600,700,800),
                          id=paste0("y",c(1,2,3,4,5,6,7,8)),
                          stringsAsFactors = F)
  bed_table <- data.frame(chr = c(1,1,1,1),
                          start = c(95,206,350,506),
                          end = c(205,405,505,705),
                          id=paste0("x",c(1,2,3,4)),
                          stringsAsFactors = F)
  
  expected_classes1 <- data.frame(row.names = "anyPosition",
                                  anyRegion=7,
                                  noMatch=1,
                                  stringsAsFactors = F)
  expected_classes2 <- data.frame(row.names = "anyRegion",
                                  anyPosition=4,
                                  noMatch=0,
                                  stringsAsFactors = F)
  
  res_obj <- intersectPositionsAndBedRegions(positions = positions,
                                             bed_table = bed_table)
  
  expect_equal(res_obj$totalPostionsInAnyRegion,7)
  expect_equal(res_obj$totalRegionsAtAnyPosition,4)
  expect_equal(res_obj$countsTable_positionsInRegionClasses,expected_classes1)
  expect_equal(res_obj$countsTable_regionsAtPositionClasses,expected_classes2)
  
})

test_that("test intersectPositionsAndBedRegions with classes", {
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
  
  expected_classes1 <- data.frame(row.names = c("C","B","A"),
                                  g2=c(1,3,0),
                                  g1=c(1,0,2),
                                  noMatch=c(1,0,0),
                                  stringsAsFactors = F)
  expected_classes2 <- data.frame(row.names = c("g2","g1"),
                                  C=c(1,1),
                                  B=c(2,0),
                                  A=c(0,2),
                                  noMatch=c(0,0),
                                  stringsAsFactors = F)
  
  res_obj <- intersectPositionsAndBedRegions(positions = positions,
                                             bed_table = bed_table)
  
  expect_equal(res_obj$totalPostionsInAnyRegion,7)
  expect_equal(res_obj$totalRegionsAtAnyPosition,4)
  expect_equal(res_obj$countsTable_positionsInRegionClasses,expected_classes1)
  expect_equal(res_obj$countsTable_regionsAtPositionClasses,expected_classes2)
  
})

test_that("test intersectBed_nonOverlapping without classes", {
  bed_table1 <- data.frame(chr = c(1,1,1),
                           start = c(100,300,500),
                           end = c(200,400,600),
                           id=paste0("x",c(1,2,3)),
                           stringsAsFactors = F)
  bed_table2 <- data.frame(chr = c(1,1,1,1),
                           start = c(250,380,550,650),
                           end = c(320,520,580,700),
                           id=paste0("y",c(1,2,3,4)),
                           stringsAsFactors = F)
  
  expected_classes1 <- data.frame(row.names = "anyRegion",
                                  anyRegion=2,
                                  noMatch=1,
                                  stringsAsFactors = F)
  expected_classes2 <- data.frame(row.names = "anyRegion",
                                  anyRegion=3,
                                  noMatch=1,
                                  stringsAsFactors = F)
  
  res_obj <- intersectBed_nonOverlapping(bed_table1 = bed_table1,
                                         bed_table2 = bed_table2)
  
  expect_equal(res_obj$totalRegions1overlappingAnyRegion2,2)
  expect_equal(res_obj$totalRegions2overlappingAnyRegion1,3)
  expect_equal(res_obj$countsTable_regions1overlappingRegion2classes,expected_classes1)
  expect_equal(res_obj$countsTable_regions2overlappingRegion1classes,expected_classes2)
  
})

test_that("test intersectBed_nonOverlapping with classes", {
  bed_table1 <- data.frame(chr = c(1,1,1),
                           start = c(100,300,500),
                           end = c(200,400,600),
                           id=paste0("x",c(1,2,3)),
                           class=c("A","B","B"),
                           stringsAsFactors = F)
  bed_table2 <- data.frame(chr = c(1,1,1,1),
                           start = c(250,380,550,650),
                           end = c(320,520,580,700),
                           id=paste0("y",c(1,2,3,4)),
                           class=c("M","M","N","N"),
                           stringsAsFactors = F)
  
  expected_classes1 <- data.frame(row.names = c("A","B"),
                                  M=c(0,2),
                                  N=c(0,1),
                                  noMatch=c(1,0),
                                  stringsAsFactors = F)
  expected_classes2 <- data.frame(row.names = c("M","N"),
                                  A=c(0,0),
                                  B=c(2,1),
                                  noMatch=c(0,1),
                                  stringsAsFactors = F)
  
  res_obj <- intersectBed_nonOverlapping(bed_table1 = bed_table1,
                                         bed_table2 = bed_table2)
  
  expect_equal(res_obj$totalRegions1overlappingAnyRegion2,2)
  expect_equal(res_obj$totalRegions2overlappingAnyRegion1,3)
  expect_equal(res_obj$countsTable_regions1overlappingRegion2classes,expected_classes1)
  expect_equal(res_obj$countsTable_regions2overlappingRegion1classes,expected_classes2)
  
})


test_that("test intersectBed without classes", {
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
  
  expected_classes1 <- data.frame(row.names = "anyRegion",
                                  anyRegion=2,
                                  noMatch=1,
                                  stringsAsFactors = F)
  expected_classes2 <- data.frame(row.names = "anyRegion",
                                  anyRegion=3,
                                  noMatch=1,
                                  stringsAsFactors = F)
  
  res_obj <- intersectBed(bed_table1 = bed_table1,
                          bed_table2 = bed_table2)
  
  expect_equal(res_obj$totalRegions1overlappingAnyRegion2,2)
  expect_equal(res_obj$totalRegions2overlappingAnyRegion1,3)
  expect_equal(res_obj$countsTable_regions1overlappingRegion2classes,expected_classes1)
  expect_equal(res_obj$countsTable_regions2overlappingRegion1classes,expected_classes2)
  
})


test_that("test intersectBed with classes", {
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
  
  expected_classes1 <- data.frame(row.names = c("A","B"),
                                  M=c(0,2),
                                  N=c(0,1),
                                  noMatch=c(1,0),
                                  stringsAsFactors = F)
  expected_classes2 <- data.frame(row.names = c("M","N"),
                                  A=c(0,0),
                                  B=c(2,1),
                                  noMatch=c(0,1),
                                  stringsAsFactors = F)
  
  res_obj <- intersectBed(bed_table1 = bed_table1,
                          bed_table2 = bed_table2)
  
  expect_equal(res_obj$totalRegions1overlappingAnyRegion2,2)
  expect_equal(res_obj$totalRegions2overlappingAnyRegion1,3)
  expect_equal(res_obj$countsTable_regions1overlappingRegion2classes,expected_classes1)
  expect_equal(res_obj$countsTable_regions2overlappingRegion1classes,expected_classes2)
  
})
