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
