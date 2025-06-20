
#' Format positions data to window
#'
#' Segment a genome into windows of a fixed size and count how many positions are
#' in each window.
#' 
#' @param positions data frame with at least two required columns: chr, position
#' @param genomev genome version to use, hg19 or hg38
#' @param windowSize window size to use. The last window at the end of each chromosome is
#' likely to be shorter so that the final window does not cover coordinates that are outside the genome.
#' @return data frame of counts of positions in windows across the whole genome.
#' Returns only the windows with non-zero data
#' @export
formatPositionDataToWindow <- function(positions,
                                       genomev="hg19",
                                       windowSize=500000,
                                       verbose=FALSE){
  
  if(genomev=="hg19" & startsWith(as.character(positions$chr[1]),"chr")) positions$chr <- substr(positions$chr,4,5)
  if(genomev=="hg38" & !startsWith(as.character(positions$chr[1]),"chr")) positions$chr <- paste0("chr",positions$chr)
  
  # get the chromosomes
  chromsTable <- getChromosomesBedTable(genomev = genomev)
  chroms <- chromsTable$chr
  chromEnds <- chromsTable$end
  names(chromEnds) <- chromsTable$chr
  
  sigCountsTable <- data.frame(matrix(0,nrow = 0,ncol = 1,
                                      dimnames = list(c(),"counts")),
                               stringsAsFactors = F)
  
  for(i in 1:nrow(positions)){
    if(verbose & (i %% 10000 == 0)) message("[info formatPositionDataToWindow] running position ",i," of ",nrow(positions))
    # i <- 1
    p <- positions[i,"position"]
    wub <- ceiling(p/windowSize)*windowSize
    wlb <- wub - windowSize + 1
    countMutation <- TRUE
    if(wlb > chromEnds[positions[i,"chr"]]){
      message("[warning formatPositionDataToWindow] position ",positions[i,"chr"],":",p," is outside chromosome ",positions[i,"chr"]," and will not be counted.")
      countMutation <- FALSE
    }else if(wub > chromEnds[positions[i,"chr"]]){
      wub <- chromEnds[positions[i,"chr"]]
    }
    if(countMutation){
      id <- paste(positions[i,"chr"],sprintf("%d",wlb),sprintf("%d",wub),sep = "_")
      if(!(id %in% rownames(sigCountsTable))){
        # add new row
        newrow <- data.frame(matrix(0,nrow = 1,ncol = 1,
                                    dimnames = list(id,"counts")),
                             stringsAsFactors = F)
        sigCountsTable <- rbind(sigCountsTable,newrow)
      }
      sigCountsTable[id,"counts"] <- sigCountsTable[id,"counts"] + 1
    }
  }
  
  return(sigCountsTable)
}



#' Format bed data to window
#'
#' Segment a genome into windows of a fixed size and for each window find
#' a set of overlapping bed regions and calculate an aggregated signal.
#' If the bed regions in the bed_table do not have a signal column, then
#' the aggregated signal for a window will simply be the proportion of the window
#' that is covered by the bed regions. This is equivalent to specifying a signal
#' of 1 for all bed regions. If a signal column is specified, then the aggregated
#' signal will be a weighted sum, where the signal of each region is multiplied
#' by the proportion of the window covered by that region.
#' 
#' @param bed_table data frame with at least three required columns: chr, start, end, and optionally signal
#' @param genomev genome version to use, hg19 or hg38
#' @param windowSize window size to use. The last window at the end of each chromosome is
#' likely to be shorter so that the final window does not cover coordinates that are outside the genome.
#' @return data frame of aggregated bed region signal for each window across the whole genome
#' @export
formatBedDataToWindow <- function(bed_table,
                                  genomev="hg19",
                                  windowSize=500000){
  # check which chromosomes have overlap if any
  overlapChroms <- checkBedRegionsOverlap(bed_table)
  if(!is.null(overlapChroms)){
    message("[error formatBedDataToWindow] bed_table regions should not overlap, ",
            "you can break them down using the function breakDownOverlappingBedRegions.")
    return(NULL)
  }
  if(!"signal" %in% colnames(bed_table)){
    # set bed signal to 1 as default to obtain the proportion of window that
    # is covered by any region
    bed_table$signal <- 1
  }
  
  # check chrom notation
  if(genomev=="hg19" & startsWith(as.character(bed_table$chr[1]),"chr")) bed_table$chr <- substr(bed_table$chr,4,5)
  if(genomev=="hg38" & !startsWith(as.character(bed_table$chr[1]),"chr")) bed_table$chr <- paste0("chr",bed_table$chr)
  
  # get the chromosomes
  chromsTable <- getChromosomesBedTable(genomev = genomev)
  chroms <- chromsTable$chr
  # sort the bed, just to be sure it's in order
  bed_table <- sortBed(bed_table)
  
  # let's start
  finalTable <- NULL
  for (chr in chroms){
    message("[info formatBedDataToWindow] formatting chromosome ",chr)
    # chr <- chroms[1]
    chrsize <- chromsTable$end[chromsTable$chr==chr]
    nwindows <- ceiling(chrsize/windowSize)
    tmp_bed_table <- bed_table[bed_table$chr==chr,,drop=F]
    tmp_region_pos <- 1
    for (wi in 1:nwindows){
      # wi <- 1
      wstart <- windowSize*(wi-1) + 1
      wend <- min(windowSize*wi,chrsize)
      windowActualSize <- wend - wstart + 1
      windowID <- paste(chr,sprintf("%d",wstart),sprintf("%d",wend),sep = "_")
      
      # prepare the new row
      newrow <- data.frame(row.names = windowID,
                           signal=0,
                           stringsAsFactors = F)
      slideWindow <- FALSE
      while(tmp_region_pos<=nrow(tmp_bed_table) & !slideWindow){
        # 1. check if the current region starts in a previous window and
        #    ends in the current window (add signal, region+1)
        # 2. check if the current region starts in a previous window and
        #    ends in a future window (add signal, slide)
        # 3. check if the current region starts and ends in this window
        #    (add signal, region+1)
        # 4. check if the current region starts in this window and ends in
        #    a future window (add signal, slide)
        # 5. check if the current region is not in this window at all (slide)
        rstart <- tmp_bed_table$start[tmp_region_pos]
        rend <- tmp_bed_table$end[tmp_region_pos]
        rsignal <- tmp_bed_table$signal[tmp_region_pos]
        # now check
        if(rstart < wstart & rend <= wend){
          # 1.
          regionActualSize <- rend - wstart + 1
          newrow$signal <- newrow$signal + rsignal * regionActualSize/windowActualSize
          tmp_region_pos <- tmp_region_pos + 1
        }else if(rstart < wstart & rend > wend){
          # 2. 
          newrow$signal <- newrow$signal + rsignal
          slideWindow <- TRUE
        }else if(rstart >= wstart & rend <= wend){
          # 3.
          regionActualSize <- rend - rstart + 1
          newrow$signal <- newrow$signal + rsignal * regionActualSize/windowActualSize
          tmp_region_pos <- tmp_region_pos + 1
        }else if(rstart >= wstart & rstart <= wend & rend > wend){
          # 4.
          regionActualSize <- wend - rstart + 1
          newrow$signal <- newrow$signal + rsignal * regionActualSize/windowActualSize
          slideWindow <- TRUE
        }else if (rstart > wend){
          #  5.
          slideWindow <- TRUE
        }else{
          message("[error formatBedDataToWindow] coding error, unexpected window/region status.")
          slideWindow <- TRUE
        }
      }
      
      # add the new row
      finalTable <- rbind(finalTable,
                          newrow)
    }
  }
  return(finalTable)
}

#' Format structural variant breakpoint position data to window
#'
#' Segment a genome into windows of a fixed size and count how many structural
#' variant breakpoints are in each window.
#' 
#' @param sv_bedpe data.frame with required columns: chrom1, start1, end1, chrom2, start2, end2
#' @param genomev genome version to use, hg19 or hg38
#' @param windowSize window size to use. The last window at the end of each chromosome is
#' likely to be shorter so that the final window does not cover coordinates that are outside the genome.
#' @return data frame of counts of positions in windows across the whole genome.
#' Returns only the windows with non-zero data
#' @export
formatBedpeBreakpointsDataToWindow <- function(sv_bedpe,
                                               genomev="hg19",
                                               windowSize=500000){
  # first of all convert the breakpoints into positions
  positions <- bedpeBreakpointsToPositions(sv_bedpe = sv_bedpe)
  return(formatPositionDataToWindow(positions = positions,
                                    genomev = genomev,
                                    windowSize = windowSize))
}


#' Sort window data
#'
#' Sort the rows of a data frame containing window data.
#' 
#' @param windowData data.frame where the row names have the format chr_start_end
#' @return sorted window data
#' @export
sortWindowData <- function(windowData){
  # get location of each row
  locationsTable <- as.data.frame(do.call(rbind,sapply(rownames(windowData),
                                                       function(x) strsplit(x,split = "_"),
                                                       simplify = T,USE.NAMES = F)),
                                  stringsAsFactors = F)
  colnames(locationsTable) <- c("chr","start","end")
  locationsTable$start <- as.numeric(locationsTable$start)
  locationsTable$end <- as.numeric(locationsTable$end)
  locationsTable <- cbind(locationsTable,order=1:nrow(locationsTable))
  locationsTable <- sortBed(bed_table = locationsTable)
  windowData <- windowData[locationsTable$order,,drop=F]
  return(windowData)
}

getWindowSizeFromRowNames <- function(windowData){
  tmpSplit <- strsplit(rownames(windowData)[1],split = "_")[[1]]
  return(as.numeric(tmpSplit[3])-as.numeric(tmpSplit[2])+1)
}

#' Merge window data
#'
#' Combine two data frames containing window data. The function expects that the
#' two data frames do not share column names, while they may share row names, where
#' the row names indicate the windows. The two data frame should use the same window
#' size. If a window is present in one data frame but missing in the second, a 
#' missing value filler will be use, NA as default.
#' 
#' @param windowData1 data.frame where the row names have the format chr_start_end
#' @param windowData2 data.frame where the row names have the format chr_start_end
#' @param missingValueFill filler to use for the missing data, NA as default
#' @return merged window data
#' @export
mergeWindowData <- function(windowData1,
                            windowData2,
                            missingValueFill=NA){
  # some columns are shared by windowData1 and windowData2, for now let's say 
  # we won't allow merge, we can improve later allowing sum/average options
  if(is.null(windowData1)) return(windowData2)
  if(is.null(windowData2)) return(windowData1)
  
  if(length(intersect(colnames(windowData1),colnames(windowData2)))>0){
    message("[error mergeWindowData] merging of windowData tables that share column names is not allowed")
    return(NULL)
  }
  windowData2_private_rownames <- setdiff(rownames(windowData2),rownames(windowData1))
  windowData2_private_colnames <- setdiff(colnames(windowData2),colnames(windowData1))
  
  if(getWindowSizeFromRowNames(windowData1)!=getWindowSizeFromRowNames(windowData2)){
    message("[warning mergeWindowData] it appears that windowData1 and windowData2 use different window sizes")
  }
  
  resultWindowData <- windowData1
  # let's start by adding the new rows to windowData1 (if any)
  if(length(windowData2_private_rownames)>0){
    resultWindowData <- rbind(resultWindowData,as.data.frame(matrix(missingValueFill,
                                                                    nrow = length(windowData2_private_rownames),
                                                                    ncol = ncol(resultWindowData),
                                                                    dimnames = list(windowData2_private_rownames,colnames(resultWindowData))),
                                                             stringsAsFactors = F))
  }
  # now add the new columns
  resultWindowData[,windowData2_private_colnames] <- missingValueFill
  # and fill in with windowData2 values
  resultWindowData[rownames(windowData2),windowData2_private_colnames] <- windowData2[rownames(windowData2),windowData2_private_colnames]
  # and return sorted
  return(sortWindowData(windowData = resultWindowData))
}

#' Get bed table from window data
#'
#' A data frame containing window data has the window coordinates as row names. 
#' This function can be used to extract the window coordinates from the row names
#' and add them to the table as chr, start, end columns.
#' 
#' @param windowData data.frame where the row names have the format chr_start_end
#' @return bed data frame
#' @export
getBedFromWindowData <- function(windowData){
  # get location of each row
  locationsTable <- as.data.frame(do.call(rbind,sapply(rownames(windowData),
                                                       function(x) strsplit(x,split = "_"),
                                                       simplify = T,USE.NAMES = F)),
                                  stringsAsFactors = F)
  colnames(locationsTable) <- c("chr","start","end")
  locationsTable$start <- as.numeric(locationsTable$start)
  locationsTable$end <- as.numeric(locationsTable$end)
  locationsTable <- cbind(locationsTable,windowData)
  return(locationsTable)
}

getBedCoordinatesFromWindowDataRownames <- function(windowData){
  # get location of each row
  locationsTable <- getBedFromWindowData(windowData = windowData)
  locationsTable <- locationsTable[,c("chr","start","end"),drop=F]
  return(locationsTable)
}
