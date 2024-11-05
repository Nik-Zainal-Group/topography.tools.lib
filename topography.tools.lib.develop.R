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

#test all
devtools::test()

#some individual tests
devtools::test(pkg = ".",filter = "base")



