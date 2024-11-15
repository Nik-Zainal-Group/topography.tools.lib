context("testing annotation functions")

test_that("test annotateBedWithGenes protein_coding only", {
  bed_table <- data.frame(chr=c(1,8),
                          start=c(10000,128500000),
                          end=c(20000,129000000),
                          id=c("testEmpty","chr8_128mb"),
                          stringsAsFactors = F)
  bed_table_annotated <- annotateBedWithGenes(bed_table = bed_table,
                                              proteinCodingOnly = TRUE,
                                              genomev = "hg19")
  expect_equal(bed_table_annotated$geneNames,c("","MYC"))
})