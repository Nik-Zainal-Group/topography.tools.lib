
#' @importFrom foreach %dopar%
NULL

#' Sort chromosomes
#'
#' Sort a list of chromosome names.
#' 
#' @param chroms characters vector with chromosome names, chr prefix is optional, so, for example, "1" and "chr1" are both acceptable
#' @param decreasing if TRUE chromosomes are sorted in decreasing order from higher chromosome number
#' @return ordered chromosome names vector
#' @export
sortChroms <- function(chroms,
                       decreasing=FALSE){
  txtchroms <- as.character(chroms)
  n <- sapply(txtchroms,function(x){
    if(startsWith(x,"chr")) x <- substr(x,4,5)
    if(x=="X"){
      return(23)
    }else if (x=="Y"){
      return(24)
    }else if (x=="M"){
      return(25)
    }else{
      return(as.numeric(x))
    }
  },USE.NAMES = F)
  return(txtchroms[order(n,decreasing = decreasing)])
}

#' Sort positions
#'
#' Sort a list of positions. Positions may be genomic coordinates of small variants or polymorphisms
#' 
#' @param positions data frame with at least two required columns: chr, position
#' @return data frame of ordered positions
#' @export
sortPositions <- function(positions){
  if(!all(colnames(positions) %in% c("chr","position"))){
    misscol <- setdiff(c("chr","position"),colnames(positions))
    message("[error sortPositions] missing columns: ",paste(misscol,collapse = ", "))
    return(NULL)
  }
  if(nrow(positions)==0){
    message("[warning sortPositions] no positions to sort")
    return(positions)
  }
  positions$chr <- as.character(positions$chr)
  sortedPositions <- NULL
  chroms <- sortChroms(unique(positions$chr))
  for (chrom in chroms){
    # chrom <- chroms[1]
    tmpTable <- positions[positions$chr==chrom,,drop=F]
    tmpTable <- tmpTable[order(tmpTable$position),,drop=F]
    sortedPositions <- rbind(sortedPositions,tmpTable)
  }
  return(sortedPositions)
}

#' Sort bed
#'
#' Sort a list of bed regions, according to the chromosome (chr column) and start
#' position (start column). If a bed region has end < start then the start and end
#' are swapped so that start <= end in the sorted bed table.
#' 
#' @param bed_table data frame with at least three required columns: chr, start, end
#' @return data frame of ordered bed regions
#' @export
sortBed <- function(bed_table){
  if(!all(colnames(bed_table) %in% c("chr","start","end"))){
    misscol <- setdiff(c("chr","start","end"),colnames(bed_table))
    message("[error sortBed] missing columns: ",paste(misscol,collapse = ", "))
    return(NULL)
  }
  if(nrow(bed_table)==0){
    message("[warning sortBed] no bed regions to sort")
    return(bed_table)
  }
  bed_table$chr <- as.character(bed_table$chr)
  sortedBed <- NULL
  chroms <- sortChroms(unique(bed_table$chr))
  for (chrom in chroms){
    # chrom <- chroms[1]
    tmpTable <- bed_table[bed_table$chr==chrom,,drop=F]
    # check if you need to reorder some start and end
    selinvert <- tmpTable$end < tmpTable$start
    if(sum(selinvert)>0){
      message("[warning sortBed] instances of start > end in ",sum(selinvert)," regions in chromosome ",chrom,
              ". Fixed by swapping start and end so that start comes first.")
      selstarts <- tmpTable$start[selinvert]
      tmpTable$start[selinvert] <- tmpTable$end[selinvert]
      tmpTable$end[selinvert] <- selstarts
    } 
    tmpTable <- tmpTable[order(tmpTable$start),,drop=F]
    sortedBed <- rbind(sortedBed,tmpTable)
  }
  return(sortedBed)
}

#' Check existence of bed regions overlaps
#'
#' Return NULL if all the bed regions in bed_table are disjoint, otherwise return
#' the list of chromosomes where the region overlaps occur.
#' 
#' @param bed_table data frame with at least three required columns: chr, start, end
#' @return NULL if no overlaps are found, or a list of chromosome names where overlaps were found
#' @export
checkBedRegionsOverlap <- function(bed_table){
  bed_table <- sortBed(bed_table)
  chroms <- unique(bed_table$chr)
  
  overlapChroms <- NULL
  for (chrom in chroms) {
    # chrom <- chroms[1]
    tmpTable <- bed_table[bed_table$chr==chrom,,drop=F]
    if(nrow(tmpTable)>1){
      existsOverlap <- any(tmpTable$end[1:(nrow(tmpTable)-1)] >= tmpTable$start[2:nrow(tmpTable)])
      if(existsOverlap) overlapChroms <- c(overlapChroms,chrom)
    }
  }
  if(!is.null(overlapChroms)){
    message("[info checkBedRegionsOverlap] Regions overlap identified in chromosomes ",paste(overlapChroms,collapse = ", "))
  }
  return(overlapChroms)
}
