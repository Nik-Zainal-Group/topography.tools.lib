context("testing sampling functions")

test_that("test randomPositionInRegions", {
  bed_table <- data.frame(chr=c(1,1,2,2),
                          start=c(10001,20001,30001,40001),
                          end=c(20000,30000,40000,50000),
                          stringsAsFactors = F)
  sampled_position <- randomPositionInRegions(samplingRegions = bed_table)
  expect_true(is.integer(sampled_position$position))
})

test_that("test resamplePositions", {
  positions <- data.frame(chr=c(1,1,2,2),
                          position=c(15000,16000,31000,32000),
                          stringsAsFactors = F)
  bed_table <- data.frame(chr=c(1,1,2,2),
                          start=c(10001,20001,30001,40001),
                          end=c(20000,30000,40000,50000),
                          stringsAsFactors = F)
  sampled_position <- resamplePositions(positions = positions,
                                        samplingRegions = bed_table,
                                        genomev = "hg19")
  # checking that the number of sampled positions for each chromosome is preserved
  expect_equal(c(2,2),as.vector(table(sampled_position$chr)))
})

test_that("test resamplePositions from hg19", {
  positions <- data.frame(chr=c(1,1,2,2),
                          position=c(15000,16000,31000,32000),
                          stringsAsFactors = F)
  sampled_position <- resamplePositions(positions = positions,
                                        genomev = "hg19")
  # checking that the number of sampled positions for each chromosome is preserved
  expect_equal(c(2,2),as.vector(table(sampled_position$chr)))
})