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

usethis::use_data(genetable_hg19,
                  genetable_hg38,
                  internal = TRUE,
                  overwrite = TRUE)

#test all
devtools::test()

#some individual tests
devtools::test(pkg = ".",filter = "base")



