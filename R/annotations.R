

#' Annotate bed regions with genes
#'
#' Given a table of bed regions, annotate each regions with the genes that overlap them.
#' 
#' @param bed_table data frame containing bed regions, with required columns chr, start, end, id and optionally class. Value in the id column must be unique.
#' @param proteinCodingOnly if TRUE then annotate only protein coding genes 
#' @param genomev genome version, hg19 or hg38
#' @param verbose set to FALSE to suppress the info messages. Warning and error messages will still be shown
#' @return annotated bed_table
#' @export
annotateBedWithGenes <- function(bed_table,
                                 proteinCodingOnly = TRUE,
                                 genomev = "hg19",
                                 verbose = TRUE){
  if(genomev=="hg19"){
    genetable <- genetable_hg19
  }else if(genomev=="hg38"){
    genetable <- genetable_hg38
  }else{
    message("[error annotatedBedWithGenes] invalid genomev, cannot annotate genes. Use hg19 or hg38.")
    return(NULL)
  }
  
  # check if only protein coding genes were requested
  if(proteinCodingOnly){
    genetable <- genetable[genetable$genetype=="protein_coding",,drop=F]
  }
  # set geneid as id, genetype as class
  genetable$id <- genetable$geneid
  genetable$class <- genetable$genetype
  rownames(genetable) <- genetable$geneid
  
  if(!startsWith(as.character(bed_table$chr[1]),prefix = "chr")) genetable$chr <- substr(genetable$chr,4,5)

  # remove genetable areas of no interest
  genetable <- genetable[genetable$chr %in% unique(bed_table$chr),,drop=F]
  genetablereduced <- NULL
  for(chr in unique(bed_table$chr)){
    minbed <- min(bed_table$start[bed_table$chr==chr])
    maxbed <- max(bed_table$end[bed_table$chr==chr])
    tmpgenetable <- genetable[genetable$chr==chr,,drop=F]
    genetablereduced <- dplyr::bind_rows(genetablereduced,tmpgenetable[!(tmpgenetable$start>maxbed | tmpgenetable$end<minbed),,drop=F])
  }
  
  corr_res <- intersectBed(bed_table1 = bed_table,
                           bed_table2 = genetablereduced,
                           verbose = verbose)
  
  bed_table <- corr_res$annotatedBedRegions1
  colnames(bed_table)[which(colnames(bed_table)=="idAnnotation")] <- "geneIds"
  colnames(bed_table)[which(colnames(bed_table)=="classAnnotation")] <- "geneTypes"
  colnames(bed_table)[which(colnames(bed_table)=="nMatches")] <- "nGenes"
  
  geneNameAnnotation <- sapply(bed_table$geneIds, function(x){
    if(is.na(x)){
      return("")
    }else{
      ids <- strsplit(x,split = ";")[[1]]
      tmpgenenames <- genetablereduced[ids,c("genename"),drop=T]
      return(paste(tmpgenenames,collapse = ";"))
    }
  },USE.NAMES = F)
  
  bed_table <- cbind(bed_table,geneNames=geneNameAnnotation)
  rownames(bed_table) <- 1:nrow(bed_table)
  
  return(bed_table)
}
