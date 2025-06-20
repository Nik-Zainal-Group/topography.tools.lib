#file with commands to set up the topography-tools-lib R package

# install.packages("devtools")
# install.packages("roxygen2")
#
# devtools::create("topography.tools.lib")

# usethis::use_package("BSgenome.Hsapiens.UCSC.hg38")
# usethis::use_package("BSgenome.Hsapiens.1000genomes.hs37d5")

devtools::document()
devtools::install()
# devtools::install(dependencies = FALSE)

genetable_hg19 <- read.table("data/genes/Hg19_Gene_List_withIDs.bed",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
genetable_hg38 <- read.table("data/genes/Hg38_Gene_List_withIDs.bed",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
samplingRegions_hg19 <- read.table("data/sampling/mappableAndNotBlacklistedRegions_hg19.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
samplingRegions_hg38 <- read.table("data/sampling/mappableAndNotBlacklistedRegions_hg38.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
timing_regions_hg19 <- read.table("data/timing/timing_table_annotated_hg19_blacklistFiltered.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
timing_regions_hg38 <- read.table("data/timing/timing_table_annotated_hg38_blacklistFiltered.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F)

usethis::use_data(genetable_hg19,
                  genetable_hg38,
                  samplingRegions_hg19,
                  samplingRegions_hg38,
                  timing_regions_hg19,
                  timing_regions_hg38,
                  internal = TRUE,
                  overwrite = TRUE)

#test all
devtools::test()

#some individual tests
devtools::test(pkg = ".",filter = "base")
devtools::test(pkg = ".",filter = "intersections")
devtools::test(pkg = ".",filter = "annotations")
devtools::test(pkg = ".",filter = "sampling")
devtools::test(pkg = ".",filter = "correlation")
devtools::test(pkg = ".",filter = "windowdata")



