
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
  if(!all(c("chr","position") %in% colnames(positions))){
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
  if(!all(c("chr","start","end") %in% colnames(bed_table))){
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



#' Break down overlapping bed regions
#'
#' If a bed_table contains overlapping bed regions, split the overlapping regions
#' and return a table of non-overlapping bed regions. An id column in the input
#' bed_table data frame is required to keep track of which returned region is derived
#' from which input region or regions. For example, given region id A, start 1, end 10,
#' and region id B, start 5, end 15, the returned table will contain three regions:
#' region id A, start 1, end 4, region id "A;B", start 5, end 10, and region id B,
#' start 11, end 15.
#' 
#' @param bed_table data frame with required columns: chr, start, end, id, and optionally signal
#' @param aggregateSignalMode can be either: mean or sum, relevant only if the signal column is present
#' @param aggregateTextColumns is a list of columns of bed_table that should be aggregated
#' @return updated bed_table with non-overlapping regions
#' @export
breakDownOverlappingBedRegions <- function(bed_table,
                                           aggregateSignalMode="mean",
                                           aggregateTextColumns=NULL){
  requiredcolumns <- c("id","chr","start","end")
  acceptedSignalModes <- c("sum","mean")
  
  if(!all(requiredcolumns %in% colnames(bed_table))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table))
    message("[error breakDownOverlappingBedRegions] missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  # check aggregateSignalMode
  if(!aggregateSignalMode %in% acceptedSignalModes){
    message("[error breakDownOverlappingBedRegions] invalid aggregateSignalMode, please use on of: ",paste(acceptedSignalModes,collapse = ", "))
    return(NULL)
  }
  # check that id is unique
  if(length(unique(bed_table$id)) < nrow(bed_table)){
    message("[error breakDownOverlappingBedRegions] id column in bed_table does not have unique values.")
    return(NULL)
  }
  # check that the requested text columns exist
  if(!is.null(aggregateTextColumns)){
    if(!all(aggregateTextColumns %in% colnames(bed_table))){
      missingcolumns <- setdiff(aggregateTextColumns,colnames(bed_table))
      message("[warning breakDownOverlappingBedRegions] ignoring requested aggregate text columns because they are missing: ",paste(missingcolumns,collapse = ", "))
      aggregateTextColumns <- setdiff(aggregateTextColumns,missingcolumns)
      if(length(aggregateTextColumns)==0) aggregateTextColumns <- NULL
    }
  }
  
  # make sure id is a character
  bed_table$id <- as.character(bed_table$id)
  # make sure it is sorted
  bed_table <- sortBed(bed_table)
  
  # check if we have the signal column
  signalColumn <- "signal" %in% colnames(bed_table)
  returncolums <- requiredcolumns
  if(signalColumn) returncolums <- c(returncolums,"signal")
  
  # check which chromosomes have overlap if any
  overlapChroms <- checkBedRegionsOverlap(bed_table)
  if(is.null(overlapChroms)){
    message("[info breakDownOverlappingBedRegions] no overlapping regions found in bed_table.")
    return(bed_table[,union(returncolums,aggregateTextColumns),drop=F])
  }else{
    # need to break things down one chromosome at a time, though only in the chromosomes where there is overlap
    finalTable <- NULL
    chroms <- unique(bed_table$chr)
    for(chrom in chroms){
      # chrom <- chroms[1]
      tmpTable <- bed_table[bed_table$chr==chrom,,drop=F]
      if(!chrom %in% overlapChroms){
        # there is no overlap in this chromosome
        message("[info breakDownOverlappingBedRegions] no overlapping regions found in chromosome ",chrom)
        rownames(tmpTable) <- paste(tmpTable$chr,sprintf("%d",tmpTable$start),sprintf("%d",tmpTable$end),sep = "_")
        finalTable <- rbind(finalTable,tmpTable[,returncolums,drop=F])
      }else{
        message("[info breakDownOverlappingBedRegions] breaking overlapping regions in chromosome ",chrom)
        # there are overlapping regions in this chromosome, so we need to break things down
        allpositions <- unique(c(tmpTable$start,tmpTable$end))
        allpositions <- allpositions[order(allpositions)]
        # allpositionsTable <- data.frame(position = allpositions,
        #                                 stringsAsFactors = F)
        # allpositionsTable$isstart <- allpositionsTable$position %in% tmpTable$start
        # allpositionsTable$isend <- allpositionsTable$position %in% tmpTable$end
        isstart <- allpositions %in% tmpTable$start
        isend <- allpositions %in% tmpTable$end
        
        newTable <- data.frame(chr=rep(chrom,length(allpositions)-1),
                               start=allpositions[1:(length(allpositions)-1)],
                               end=allpositions[2:length(allpositions)],
                               stringsAsFactors = F)
        # when starting from an end position, we need to add 1
        newTable$start[isend[1:(length(allpositions)-1)]] <- newTable$start[isend[1:(length(allpositions)-1)]] + 1
        # when ending in a start position, we need to subtract 1
        newTable$end[isstart[2:length(allpositions)]] <- newTable$end[isstart[2:length(allpositions)]] - 1
        # need to remove cases in which starts becomes > than end
        newTable <- newTable[newTable$start<=newTable$end,,drop=F]
        
        # if there are cases where start and end happen at the same position,
        # they need to be added as additional regions
        bothstartandend <- allpositions[isstart & isend]
        if(length(bothstartandend)>0){
          newTable <- rbind(newTable,data.frame(chr=rep(chrom,length(bothstartandend)),
                                                start=bothstartandend,
                                                end=bothstartandend,
                                                stringsAsFactors = F))
        }
        # sort? perhaps there is no need
        newTable <- sortBed(newTable)
        # now we need to find all the original regions that overlap each new region
        # discard new segments with no overlap
        # aggregate signal if one or more segments overlap and include their ids
        # brute force approach will check each original region against all new regions
        # I can speed up if I exclude new regions from the check once they end before 
        # the start of an old region, though I need a quick access using ids as row names
        # 
        rownames(newTable) <- paste(newTable$chr,sprintf("%d",newTable$start),sprintf("%d",newTable$end),sep = "_")
        copyOfNewTable <- newTable
        newTable$hasOverlap <- FALSE
        newTable$signalsum <- 0
        newTable$noverlaps <- 0
        newTable$id <- NA
        
        # find the regions that can be matched directly to non overlapping regions
        tmpTable$bedid <- paste(tmpTable$chr,sprintf("%d",tmpTable$start),sprintf("%d",tmpTable$end),sep = "_")
        tmpTable$matchedid <- tmpTable$bedid %in% rownames(newTable)
        countMatchedids <- table(tmpTable$bedid[tmpTable$matchedid])
        uniquelyMatchedIds <- names(countMatchedids)[countMatchedids==1]
        if(length(uniquelyMatchedIds)>0){
          message("[info breakDownOverlappingBedRegions] found ", length(uniquelyMatchedIds)," bed segments (",
                  sprintf("%.2f",length(uniquelyMatchedIds)/nrow(tmpTable)*100),"%) that do not need to be broken down, they will be processed quickly.")
          # we have some uniquely matched ids that we can deal with in one go
          uniqueIDtmpTable <- tmpTable[tmpTable$bedid %in% uniquelyMatchedIds,]
          rownames(uniqueIDtmpTable) <- uniqueIDtmpTable$bedid
          newTable[uniquelyMatchedIds,"hasOverlap"] <- TRUE
          newTable[uniquelyMatchedIds,"noverlaps"] <- 1
          newTable[uniquelyMatchedIds,"id"] <- uniqueIDtmpTable[uniquelyMatchedIds,"id"]
          if(signalColumn) newTable[uniquelyMatchedIds,"signalsum"] <- uniqueIDtmpTable[uniquelyMatchedIds,"signal"]
          # update tmpTable
          tmpTable <- tmpTable[!tmpTable$bedid %in% uniquelyMatchedIds,,drop=F]
          if(nrow(tmpTable)>0){
            # make sure it is sorted
            tmpTable <- sortBed(tmpTable)
          }
        }
        
        if(nrow(tmpTable)>0){
          percmilestone <- 10
          for (i in 1:nrow(tmpTable)){
            # i <- 1
            if(i==nrow(tmpTable)){
              message("[info breakDownOverlappingBedRegions] chromosome ",chrom," progress: 100%")
            }else if(i/nrow(tmpTable)*100>=percmilestone) {
              message("[info breakDownOverlappingBedRegions] chromosome ",chrom," progress: ",percmilestone,"%")
              percmilestone <- percmilestone+10
            }
            if(tmpTable[i,"matchedid"]){
              # non need to look for overlapping new regions, we have one matching exactly
              ids <- tmpTable[i,"bedid"]
              newTable[ids,"hasOverlap"] <- TRUE
              newTable[ids,"noverlaps"] <- newTable[ids,"noverlaps"] + 1
              if(is.na(newTable[ids,"id"])){
                newTable[ids,"id"] <- tmpTable[i,"id"]
              }else{
                newTable[ids,"id"] <- paste(c(newTable[ids,"id"],tmpTable[i,"id"]),collapse = ";")
              }
              if(signalColumn) {
                newTable[ids,"signalsum"] <- newTable[ids,"signalsum"] + tmpTable[i,"signal"]
              }
            }else{
              # need to check for overlapping new regions
              # discard new regions that are lower than start
              copyOfNewTable <- copyOfNewTable[(copyOfNewTable$end >= tmpTable$start[i]),,drop=F]
              overlap <- !(copyOfNewTable$start > tmpTable$end[i] | copyOfNewTable$end < tmpTable$start[i])
              overlapTable <- copyOfNewTable[overlap,,drop=F]
              if(nrow(overlapTable)>0){
                ids <- paste(overlapTable$chr,sprintf("%d",overlapTable$start),sprintf("%d",overlapTable$end),sep = "_")
                newTable[ids,"hasOverlap"] <- TRUE
                newTable[ids,"noverlaps"] <- newTable[ids,"noverlaps"] + 1
                for(idsi in ids){
                  if(is.na(newTable[idsi,"id"])){
                    newTable[idsi,"id"] <- tmpTable[i,"id"]
                  }else{
                    newTable[idsi,"id"] <- paste(c(newTable[idsi,"id"],tmpTable[i,"id"]),collapse = ";")
                  }
                }
                if(signalColumn) {
                  newTable[ids,"signalsum"] <- newTable[ids,"signalsum"] + tmpTable[i,"signal"]
                }
              }
            }
            
          }
        }
        
        # now remove lines if noverlaps is 0
        newTable <- newTable[newTable$hasOverlap,,drop=F]
        if(signalColumn){
          if(aggregateSignalMode=="sum"){
            newTable$signal <- newTable$signalsum
          }else if(aggregateSignalMode=="mean"){
            newTable$signal <- newTable$signalsum/newTable$noverlaps
          }
        }
        finalTable <- rbind(finalTable,newTable[,returncolums,drop=F])
        
      }
    }
    
    if(!is.null(aggregateTextColumns)){
      message("[info breakDownOverlappingBedRegions] aggregating requested text columns...")
      # add the aggregated text columns requested
      aggregatedTextColumns <- sapply(finalTable$id,function(x){
        ids <- strsplit(x,split = ";")[[1]]
        tmpBed <- bed_table[bed_table$id %in% ids,aggregateTextColumns,drop=F]
        apply(tmpBed, 2, function(x) paste(x,collapse = ";"))
      })
      if(length(aggregateTextColumns)>1) {
        aggregatedTextColumns <- t(aggregatedTextColumns)
      }else{
        aggregatedTextColumns <- data.frame(aggregatedTextColumns,
                                            stringsAsFactors = F)
        colnames(aggregatedTextColumns) <- aggregateTextColumns
      }
      finalTable <- cbind(finalTable,aggregatedTextColumns)
    }
    
    message("[info breakDownOverlappingBedRegions] done.")
    return(finalTable)
  }
}


#' Assign bed regions to non-overlapping sets
#'
#' Given a table with bed regions, assign each region to a set, so that each
#' set contains no overlapping bed regions. If the regions in the input bed_table
#' do not overlap, then they will be all assigned to the same set. For example,
#' given region id A, start 1, end 10, and region id B, start 5, end 15, the 
#' function returns the list list(A=1,B=2).
#' 
#' @param bed_table data frame with required columns: chr, start, end, id, and optionally signal
#' @param brokenDown_bed_table This should be the result of breakDownOverlappingBedRegions(bed_table,...),
#' in case this has already been computed. If left NULL (the default), breakDownOverlappingBedRegions(bed_table)
#' is called by this function with default parameters
#' @return list object assigning a set number (integer) to each bed region id
#' @export
assignBedRegionsToNonOverlappingSets <- function(bed_table,
                                                 brokenDown_bed_table=NULL){
  if(is.null(brokenDown_bed_table)){
    message("[info assignBedRegionsToNonOverlappingSets] running breakDownOverlappingBedRegions...")
    brokenDown_bed_table <- breakDownOverlappingBedRegions(bed_table = bed_table)
  }
  
  levelAssignmentList <- list()
  for(i in 1:nrow(brokenDown_bed_table)){
    # i <- 1
    segmentIds <- strsplit(x = brokenDown_bed_table[i,"id"],split = ";")[[1]]
    # consider only ids of segments that are fully in the region
    segmentIds <- segmentIds[segmentIds %in% bed_table$id]
    if(length(segmentIds)>0){
      assignedIds <- segmentIds[segmentIds %in% names(levelAssignmentList)]
      if(length(assignedIds)==0){
        for (j in 1:length(segmentIds)) {
          id <- segmentIds[j]
          levelAssignmentList[[as.character(id)]] <- j
        }
      }else{
        # some segments have already been assigned, need to find the first free position
        # for each unassigned id I have
        unassignedIds <- setdiff(segmentIds,assignedIds)
        if(length(unassignedIds)>0){
          for (j in 1:length(unassignedIds)) {
            # j <- 1
            id <- unassignedIds[j]
            busyLanes <- sapply(assignedIds,function(x) levelAssignmentList[[as.character(x)]],USE.NAMES = F)
            found <- FALSE
            lane <- 1
            while (!found) {
              if(!lane %in% busyLanes){
                found <- TRUE
              }else{
                lane <- lane + 1
              }
            }
            levelAssignmentList[[as.character(id)]] <- lane
            
          }
        }
      }
    }
  }
  return(levelAssignmentList)
}


#' Plot broken down bed regions
#'
#' Given a table with bed regions, plot the bed regions within a given genomic
#' location (chrom:pstart-pend) so that overlapping regions are shown on separate
#' lines. This plot function is useful to visualise overlapping bed regions.
#' This function calls breakDownOverlappingBedRegions (if brokenDown_bed_table is NULL)
#' and assignBedRegionsToNonOverlappingSets, and the results from these function
#' calls are returned in a return object list.
#' 
#' @param bed_table data frame with required columns: chr, start, end, id, and optionally signal
#' @param brokenDown_bed_table This should be the result of breakDownOverlappingBedRegions(bed_table,...),
#' in case this has already been computed. If left NULL (the default), breakDownOverlappingBedRegions(bed_table)
#' is called by this function with default parameters
#' @param filename optional file name to plot to file, use .pdf extension
#' @param chrom chromosome location to plot
#' @param pstart genomic start location to plot
#' @param pend genomic end location to plot
#' @param region_colour plot colour for regions in bed_table
#' @param segments_colour plot colour of non-overlapping segment regions in brokenDown_bed_table
#' @return list object with the results of breakDownOverlappingBedRegions and assignBedRegionsToNonOverlappingSets function calls
#' @export
plotBrokenDownBedRegions <- function(bed_table,
                                     brokenDown_bed_table=NULL,
                                     filename=NULL,
                                     chrom,
                                     pstart,
                                     pend,
                                     region_colour="#0067a5",
                                     segments_colour="#F38400"){
  
  if(is.null(brokenDown_bed_table)){
    message("[info plotBrokenDownBedRegions] running breakDownOverlappingBedRegions...")
    brokenDown_bed_table <- breakDownOverlappingBedRegions(bed_table = bed_table)
  }
  
  # make sure ids are characters
  bed_table$id <- as.character(bed_table$id)
  brokenDown_bed_table$id <- as.character(brokenDown_bed_table$id)
  # plot only segment in a given region of interest (include partial overlap)
  selection <- !(bed_table$start > pend | bed_table$end < pstart) & bed_table$chr==chrom
  if(sum(selection)==0){
    message("[warning plotBrokenDownBedRegions] nothing to plot in the requested region.")
    return(NULL)
  }
  bed_table <- bed_table[selection,,drop=F]
  # if there is something in the bed_table then there must be something in the brokenDown_bed_table
  selection <- !(brokenDown_bed_table$start > pend | brokenDown_bed_table$end < pstart) & brokenDown_bed_table$chr==chrom
  brokenDown_bed_table <- brokenDown_bed_table[selection,,drop=F]
  levelAssignmentList <- assignBedRegionsToNonOverlappingSets(bed_table = bed_table,
                                                              brokenDown_bed_table = brokenDown_bed_table)
  
  # OK now I should be able to print
  ymax <- 2 + max(unlist(levelAssignmentList))
  # xmin <- min(bed_table$start)
  # xmax <- max(bed_table$end)
  xmin <- pstart
  xmax <- pend
  halftick <- 0.15
  
  
  if(!is.null(filename)){
    cairo_pdf(filename = filename,width = 9,height = 5)
  }
  
  par(mar=c(5,4,4,3))
  plot(NA,bty="n",
       ylim=c(0,ymax),
       xlim=c(xmin,xmax),
       ylab="",
       yaxt = "n",
       xlab="position")
  abline(v = unique(c(bed_table$start,bed_table$end)),lty=3)
  for(i in 1:nrow(bed_table)){
    # i <- 1
    ypos <- ymax - levelAssignmentList[[as.character(bed_table$id[i])]]
    lines(x=c(bed_table$start[i],bed_table$end[i]),
          y=rep(ypos,2),
          col=region_colour,
          lwd=2)
    lines(x=rep(bed_table$start[i],2),
          y=c(ypos+halftick,ypos-halftick),
          col=region_colour,
          lwd=2)
    lines(x=rep(bed_table$end[i],2),
          y=c(ypos+halftick,ypos-halftick),
          col=region_colour,
          lwd=2)
  }
  for(i in 1:nrow(brokenDown_bed_table)){
    # i <- 1
    ypos <- 1
    lines(x=c(brokenDown_bed_table$start[i],brokenDown_bed_table$end[i]),
          y=rep(ypos,2),
          col=segments_colour,
          lwd=2)
    lines(x=rep(brokenDown_bed_table$start[i],2),
          y=c(ypos+halftick,ypos-halftick),
          col=segments_colour,
          lwd=2)
    lines(x=rep(brokenDown_bed_table$end[i],2),
          y=c(ypos+halftick,ypos-halftick),
          col=segments_colour,
          lwd=2)
  }
  
  # close the file
  if(!is.null(filename)) dev.off()
  
  # return the broken down table and the assignment level
  returnObj <- list()
  returnObj$levelAssignmentList <- levelAssignmentList
  returnObj$brokenDown_bed_table <- brokenDown_bed_table
  return(returnObj)
}

#' Compute inter-mutational distance of a set of positions
#'
#' This function returns both the left and right IMD of each given position
#' as well as the average IMD.
#' 
#' @param positions data frame with required columns: chr and position 
#' @return left, right and average inter-mutational distance
#' @export
getIMD <- function(positions){
  # some checks
  requiredcolumns_pos <- c("chr","position")
  if(!all(requiredcolumns_pos %in% colnames(positions))){
    missingcolumns <- setdiff(requiredcolumns_pos,colnames(positions))
    message("[error getIMD] missing required columns in positions table: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  
  # compute the IMD, for each chromosome separately, sort first
  positions <- sortPositions(positions = positions)
  chroms <- unique(positions$chr)
  new_positions <- NULL
  for(chrom in chroms){
    # chrom <- chroms[1]
    chr_positions <- positions[positions$chr==chrom,,drop=F]
    if(nrow(chr_positions)>0){
      # get ready
      chr_positions$leftIMD <- NA
      chr_positions$rightIMD <- NA
      chr_positions$aveIMD <- NA
      # there is at least one row, so if it is only one row IMD=NA
      # and just add to final table, otherwise we calculate the IMD and fill the table
      if(nrow(chr_positions)>1){
        IMD <- chr_positions$position[2:nrow(chr_positions)] - chr_positions$position[1:(nrow(chr_positions)-1)]
        chr_positions[2:nrow(chr_positions),"leftIMD"] <- IMD
        chr_positions[1:(nrow(chr_positions)-1),"rightIMD"] <- IMD
        chr_positions[,"aveIMD"] <- apply(chr_positions[,c("leftIMD","rightIMD")],1,mean,na.rm=T)
      }
      new_positions <- rbind(new_positions,chr_positions)
    }
  }
  return(new_positions)
}



#' Merge adjacent bed regions
#'
#' Given a table with bed regions, merge regions that are adjacent and return
#' the updated table. The bed regions should not be overlapping. If there are
#' overlapping bed regions in the bed_table, consider running breakDownOverlappingBedRegions
#' to obtain non-overlapping segments and then run this function to merge the
#' adjacent segments. For example, given region id A, start 1, end 10, and region id B,
#' start 11, end 15, the function returns the merged region with id "A;B", start 1, end 15.
#' 
#' 
#' @param bed_table data frame with required columns: chr, start, end, id
#' @return updated bed_table with merged adjacent regions
#' @export
mergeAdjacentBedRegions <- function(bed_table){
  # check required columns
  requiredcolumns <- c("id","chr","start","end")
  if(!all(requiredcolumns %in% colnames(bed_table))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table))
    message("[error mergeAdjacentBedRegions] missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  # check overlaps
  res_check <- checkBedRegionsOverlap(bed_table=bed_table)
  if(!is.null(res_check)){
    message("[error mergeAdjacentBedRegions] bed_table regions should not overlap, ",
            "you can break them down using the function breakDownOverlappingBedRegions.")
    return(NULL)
  }
  # now sort
  bed_table <- sortBed(bed_table = bed_table)
  # now merge
  newTable <- NULL
  chroms <- unique(bed_table$chr)
  for(chrom in chroms){
    # chrom <- chroms[1]
    chrTable <- bed_table[bed_table$chr==chrom,requiredcolumns,drop=F]
    if(nrow(chrTable)==1){
      # only one row so nothing to merge with
      newTable <- rbind(newTable,chrTable)
    }else{
      # sure more than one
      regdist <- chrTable$start[2:(nrow(chrTable))] - chrTable$end[1:(nrow(chrTable)-1)]
      mergepos <- which(regdist==1)
      if(length(mergepos)==0){
        # nothing to merge in this chrom
        newTable <- rbind(newTable,chrTable)
      }else{
        # find where to merge and where to copy rows
        # group the adjacent windows in case of multiple consecutive merges
        mergeGroups <- list()
        if(length(mergepos)==1) {
          mergeGroups[["1"]] <- mergepos
        }else{
          currentMergeGroup <- 1
          for(i in 1:length(mergepos)){
            mergeGroups[[as.character(currentMergeGroup)]] <- c(mergeGroups[[as.character(currentMergeGroup)]],mergepos[i])
            if(i<length(mergepos)){
              if(mergepos[i+1]-mergepos[i]>1) currentMergeGroup <- currentMergeGroup + 1
            }
          }
        }
        # check for rows to copy before the merges
        if(min(mergeGroups[["1"]])>1){
          # yes we copy from 1 to the first merge
          newTable <- rbind(newTable,chrTable[1:(min(mergeGroups[["1"]])-1),,drop=F])
        }
        # now get to merge
        for (ni in 1:length(mergeGroups)){
          # ni <- 1
          n <- names(mergeGroups)[ni]
          mergeRows <- min(mergeGroups[[n]]):(max(mergeGroups[[n]])+1)
          newTable <- rbind(newTable,data.frame(id=paste(chrTable[mergeRows,"id"],collapse = ";"),
                                                chr=chrom,
                                                start=min(chrTable[mergeRows,"start"]),
                                                end=max(chrTable[mergeRows,"end"]),
                                                stringsAsFactors = F))
          # now I should check if there are rows to copy after the mergeGroups or in between
          if(ni==length(mergeGroups)){
            # ok we are at the end
            startingrow <- max(mergeGroups[[n]])+2
            if(startingrow<=nrow(chrTable)){
              # and we got something to add
              newTable <- rbind(newTable,chrTable[startingrow:nrow(chrTable),,drop=F])
            }
          }else{
            # there is another merge later
            startingrow <- max(mergeGroups[[n]])+2
            endingrow <- min(mergeGroups[[names(mergeGroups)[ni+1]]])-1
            if(startingrow<=endingrow){
              newTable <- rbind(newTable,chrTable[startingrow:endingrow,,drop=F])
            }
          }
        }
      }
    }
    
  }
  return(newTable)
}

#' Get chromosomes bed table
#'
#' Return a bed table with the list of chromosomes.
#' 
#' 
#' @param genomev hg19 or hg38
#' @return chromosomes bed table
#' @export
getChromosomesBedTable <- function(genomev){
  # select reference genome
  if(genomev=="hg19"){
    expected_chroms <- paste0(c(seq(1:22),"X","Y"))
    genomeSeq <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  }else if(genomev=="hg38"){
    expected_chroms <- paste0("chr",c(seq(1:22),"X","Y"))
    genomeSeq <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }
  
  # get chrom lengths info
  chromsTable <- as.data.frame(GenomeInfoDb::seqinfo(genomeSeq))
  chromsTable <- chromsTable[expected_chroms,]
  
  # set up table
  regions_table <- data.frame(chr=rownames(chromsTable),
                              start=rep(1,nrow(chromsTable)),
                              end=chromsTable$seqlengths,
                              stringsAsFactors = F)
  return(regions_table)
}


#' Remove N from bed table
#'
#' Given a table of bed regions and a reference genome, trim and split the bed
#' bed regions to remove reference genome N positions.
#' 
#' 
#' @param bed_table data frame with required columns: chr, start, end
#' @param genomev hg19 or hg38
#' @return updated bed_table
#' @export
trimNfromBed <- function(bed_table,
                         genomev,
                         verbose=FALSE){
  # check required columns
  requiredcolumns <- c("chr","start","end")
  if(!all(requiredcolumns %in% colnames(bed_table))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table))
    message("[error trimNfromBed] missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  
  # check for chr prefix
  chrPrefix <- FALSE
  if(startsWith(x = as.character(bed_table$chr[1]),prefix = "chr")) chrPrefix <- TRUE
  
  # select reference genome
  if(genomev=="hg19"){
    expected_chroms <- paste0(c(seq(1:22),"X","Y"))
    genomeSeq <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
    if(chrPrefix) bed_table$chr <- substr(bed_table$chr,4,5)
  }else if(genomev=="hg38"){
    expected_chroms <- paste0("chr",c(seq(1:22),"X","Y"))
    genomeSeq <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    if(!chrPrefix) bed_table$chr <- paste0("chr",bed_table)
  }
  
  # great, now let's check for N and remove them
  regions_table_final <- NULL
  for (i in 1:nrow(bed_table)) {
    # i <- 1
    if((i %% 50 == 0) & verbose){
      message("[info trimNfromBed] processing row ",i," of ",nrow(bed_table))
    }
    currentSeq <- as.character(BSgenome::getSeq(genomeSeq, as.character(bed_table$chr[i]), start=bed_table$start[i], bed_table$end[i]))
    
    # does it contain N?
    hasN <- grepl(pattern = "N",x = currentSeq,fixed = TRUE)
    if(!hasN){
      # just add row as it is
      regions_table_final <- rbind(regions_table_final,bed_table[i,,drop=F])
    }else{
      if(verbose) message("[info trimNfromBed] Found N in row ",i,": splitting")
      # I need all the positions where N is, then use distance from next to find the regions
      # whenever I have two non-consecutive N then I have a region
      isN <- strsplit(currentSeq,split = "")[[1]]=="N"
      if(!all(isN)){
        positions <- bed_table$start[i]:bed_table$end[i]
        positions <- positions[isN]
        if(positions[1]!=bed_table$start[i]) positions <- c(bed_table$start[i]-1,positions)
        if(positions[length(positions)]!=bed_table$end[i]) positions <- c(positions,bed_table$end[i]+1)
        pdist <- positions[2:length(positions)]-positions[1:(length(positions)-1)]
        segmentsPos <- which(pdist>1)
        newrows <- data.frame(chr=rep(bed_table$chr[i],length(segmentsPos)),
                              start=positions[segmentsPos]+1,
                              end=positions[segmentsPos+1]-1,
                              stringsAsFactors = F)
        if(verbose) message("[info trimNfromBed] -> row ",i," split into ",nrow(newrows))
        regions_table_final <- rbind(regions_table_final,newrows)
      }else{
        if(verbose) message("[info trimNfromBed] -> row ",i," is 100% N")
      }
    }
  }
  
  return(regions_table_final)
}


#' Remove N from bed table
#'
#' Given a table of bed regions and a reference genome, trim and split the bed
#' bed regions to remove reference genome N positions.
#' 
#' 
#' @param bed_table data frame with required columns: chr, start, end, signal
#' @param bed_table2 data frame with required columns: chr, start, end, signal
#' @param fileout name of the output file for the plot, use .pdf extension
#' @param pchr chromosome name of the region to plot
#' @param pstart left boundary position of the genomic region to plot
#' @param pend right boundary position of the genomic region to plot
#' @param signalColour colour of the signal line
#' @param signalColour2 colour of the second signal line from bed_table2
#' @param highlightRegions highlightRegions is a list of bed tables (chr, start, end)
#' @param highlightRegionsColours highlightRegionsColours is a list of colours.
#' The names of the list need to match the names of the highlightRegions list
#' @param highlightPositions highlightPositions is a list of position tables (chr, position, optional text)
#' @param highlightPositionsColours highlightPositionsColours is a list of colours.
#' The names of the list need to match the names of the highlightPositions list
#' @param plotGenes if TRUE the coding genes in the regions will be plotted on top
#' @param genomev hg19 or hg38
#' @param main title of the plot
#' @param lwd line width of the plot
#' @param cexlabels scaling parameter for the labels
#' @param ylabel ylabel for the bed_table signal
#' @param ylabel2 ylabel for the bed_table2 signal
#' @export
plotBedSignalRegion <- function(bed_table,
                                bed_table2=NULL,
                                fileout=NULL,
                                pchr,
                                pstart,
                                pend,
                                signalColour="black",
                                signalColour2="grey",
                                highlightRegions=NULL,
                                highlightRegionsColours=NULL,
                                highlightPositions=NULL,
                                highlightPositionsColours=NULL,
                                plotGenes=TRUE,
                                genomev="hg19",
                                main="",
                                lwd=1.5,
                                cexlabels=1,
                                ylabel="signal",
                                ylabel2="signal"){
  
  # check overlaps
  res_check <- checkBedRegionsOverlap(bed_table=bed_table)
  if(!is.null(res_check)){
    message("[warning plotBedSignalRegion] some bed_table regions overlap. Signal in ",
            "overlapping segments will be summed. If you prefer to resolve overlapping segments ",
            "yourself, you can use the function breakDownOverlappingBedRegions.")
  }
  if(!is.null(bed_table2)){
    res_check <- checkBedRegionsOverlap(bed_table=bed_table2)
    if(!is.null(res_check)){
      message("[warning plotBedSignalRegion] some bed_table2 regions overlap. Signal in ",
              "overlapping segments will be summed. If you prefer to resolve overlapping segments ",
              "yourself, you can use the function breakDownOverlappingBedRegions.")
    }
  }
  
  # set plot parameters in inch
  mbottom <- 0.9
  mtop <- 1
  mleft <- 1
  mright <- 0.5
  if(!is.null(highlightRegions)){
    mright <- 0.1+(max(strwidth(names(highlightRegions),units = "inch",cex = 1,ps = par(ps=12))))
  }
  if(!is.null(bed_table2)){
    mright <- max(mright,1)
  }
  datawidth <- 6
  dataheight <- 1.2
  highlightRegionHeightInch <- 0.2
  basesPerInch <- (pend-pstart)/datawidth
  
  # convert to character just in case
  pchr <- as.character(pchr)
  
  # select the data that is overlapping the region of interest
  regionBed <- bed_table[bed_table$chr==pchr,,drop=F]
  if(nrow(regionBed)>0){
    selection <- !(regionBed$start > pend | regionBed$end < pstart)
    regionBed <- regionBed[selection,,drop=F]
  }
  # same for bed_table2 if any
  regionBed2 <- NULL
  if(!is.null(bed_table2)){
    regionBed2 <- bed_table2[bed_table2$chr==pchr,,drop=F]
    if(nrow(regionBed2)>0){
      selection <- !(regionBed2$start > pend | regionBed2$end < pstart)
      regionBed2 <- regionBed2[selection,,drop=F]
    }
  }
  
  # before we plot, let's check that we have data for all positions in the region
  # we can do that by adding one bed region from start to end with zero signal
  # and break down the overlap
  regionBed <- regionBed[,c("chr", "start", "end", "signal"),drop=F]
  plotBed <- data.frame(chr=pchr,
                        start=pstart,
                        end=pend,
                        signal=0,
                        stringsAsFactors = F)
  regionBed <- rbind(regionBed,plotBed)
  regionBed$id <- 1:nrow(regionBed)
  res_bd <- breakDownOverlappingBedRegions(bed_table = regionBed,
                                           aggregateSignalMode = "sum")
  # remove segments that are outside the plot region
  plotregionId <- as.character(nrow(regionBed))
  selectRows <- sapply(res_bd$id,function(id){
    ids <- strsplit(id,split = ";")[[1]]
    return(plotregionId %in% ids)
  },USE.NAMES = F)
  res_bd <- res_bd[selectRows,,drop=F]
  
  # same for bed_table2
  res_bd2 <- NULL
  if(!is.null(regionBed2)){
    regionBed2 <- regionBed2[,c("chr", "start", "end", "signal"),drop=F]
    plotBed <- data.frame(chr=pchr,
                          start=pstart,
                          end=pend,
                          signal=0,
                          stringsAsFactors = F)
    regionBed2 <- rbind(regionBed2,plotBed)
    regionBed2$id <- 1:nrow(regionBed2)
    res_bd2 <- breakDownOverlappingBedRegions(bed_table = regionBed2,
                                              aggregateSignalMode = "sum")
    # remove segments that are outside the plot region
    plotregionId <- as.character(nrow(regionBed2))
    selectRows <- sapply(res_bd2$id,function(id){
      ids <- strsplit(id,split = ";")[[1]]
      return(plotregionId %in% ids)
    },USE.NAMES = F)
    res_bd2 <- res_bd2[selectRows,,drop=F]
  }
  
  # get Genes
  genetable <- NULL
  levelAssignmentList <- NULL
  if(plotGenes){
    # read the gene files
    if(genomev=="hg19"){
      genetable <- genetable_hg19
    }else if(genomev=="hg38"){
      genetable <- genetable_hg38
    }else{
      message("[error plotBedSignalRegion] invalid genomev, cannot annotate genes. Use hg19 or hg38.")
      return(NULL)
    }
    # only protein coding
    genetable <- genetable[genetable$genetype=="protein_coding",,drop=F]
    if(!startsWith(as.character(bed_table$chr[1]),prefix = "chr")) genetable$chr <- substr(genetable$chr,4,5)
    # select only the relevant part of the table
    genetable <- genetable[genetable$chr==pchr,,drop=F]
    genetable <- genetable[!(genetable$start > pend | genetable$end < pstart),,drop=F]
    if(nrow(genetable)>0){
      genetable$id <- 1:nrow(genetable)
      tmpgenetable <- genetable
      # when assigning levels I need to consider if the text with the gene name also overlaps
      geneNamesBasesSize <- basesPerInch*strwidth(tmpgenetable$genename,units = "inch",cex = 0.5,ps = par(ps=12))
      for(j in 1:nrow(tmpgenetable)) {
        genemiddlepoint <- (max(tmpgenetable$start[j],pstart)+min(tmpgenetable$end[j],pend))/2
        tmpgenetable$start[j] <- round(min(tmpgenetable$start[j],genemiddlepoint-geneNamesBasesSize[j]/2))
        tmpgenetable$end[j] <- round(max(tmpgenetable$end[j],genemiddlepoint+geneNamesBasesSize[j]/2))
      }
      genes_res <- breakDownOverlappingBedRegions(bed_table = tmpgenetable)
      levelAssignmentList <- assignBedRegionsToNonOverlappingSets(bed_table = tmpgenetable,
                                                                  brokenDown_bed_table = genes_res)
    }
  }
  
  
  # infer more parameters for plotting
  signalMin <- min(0,min(res_bd$signal))
  signalMax <- max(1,max(res_bd$signal))
  ydatagap <- 0.05*(signalMax-signalMin)
  ydatagapTop <- 0.15*(signalMax-signalMin)
  ylimData <- c(signalMin-ydatagap,signalMax+ydatagapTop)
  ySize <- ylimData[2]-ylimData[1]
  highlightRegionHeightData <- ySize*highlightRegionHeightInch/dataheight
  nhighlightregions <- 0
  if(!is.null(highlightRegions)) {
    nhighlightregions <- length(highlightRegions)
  }
  nhighlightpositions <- 0
  if(!is.null(highlightPositions)) {
    nhighlightpositions <- length(highlightPositions)
  }
  # check if we need to plot genes
  nGeneLayers <- 0
  if(!is.null(levelAssignmentList)) {
    nGeneLayers <- 1 + max(unlist(levelAssignmentList))
  }
  genesGap <- 0.6*highlightRegionHeightData
  # determine actual ylim
  ylim <- c(ylimData[1],ylimData[2]+(nhighlightregions+nGeneLayers)*highlightRegionHeightData)
  # determine position of text for highlight regions
  xtextgap <- (pend-pstart)*0.01
  xtextpos <- (pend+xtextgap)/1e6
  
  # determine plot size in inch
  pwidth <- mleft + datawidth + mright
  pheight <- mbottom + dataheight + (nhighlightregions+nGeneLayers)*highlightRegionHeightInch + mtop
  
  #  ok now plot
  if(!is.null(fileout)) cairo_pdf(filename = fileout,width = pwidth,height = pheight)
  par(mai=c(mbottom,mleft,mtop,mright),mgp=c(2.5,0.9,0))
  xlab <- ifelse(startsWith(as.character(pchr),prefix = "chr"),substr(pchr,4,8),pchr)
  plot(NA,
       xlim=c(pstart,pend)/1e6,
       bty="n",
       ylim=ylim,
       main=main,
       ylab="",
       las=1,
       xaxs="i",
       yaxs="i",
       yaxt="n",
       xlab=paste0("chromosome ",xlab," (Mb)"))
  # ylabel
  xylabelpos <- (pstart - 0.1*(pend-pstart))/1e6
  text(x = xylabelpos,
       y = (ylimData[1]+ylimData[2])/2,
       labels=ylabel,srt=90,adj=0.5,xpd=T,col=signalColour)
  lines(x=c(pstart,pend)/1e6,y=c(ylimData[1],ylimData[1]),col="black",xpd=T)
  lines(x=c(pstart,pstart)/1e6,y=c(ylimData[1],ylimData[2]),col="black",xpd=T)
  # ok figure out the axis from 0 to something just below ylimData[2]
  yaxisgaps <- c(1,2,5)
  gapi <- 1
  gapscale <- 1
  gapfound <- FALSE
  gaptarget <- ylimData[2]/2
  yaxisfinalgap <- NULL
  while (!gapfound) {
    currentgap <- yaxisgaps[gapi]/gapscale
    # check if we are below gaptarget
    if(currentgap<=gaptarget){
      # if the next one up is greater than target we are done
      if(gapi==3){
        nextgap <- yaxisgaps[1]/gapscale*10
      }else{
        nextgap <- yaxisgaps[gapi+1]/gapscale
      }
      if(nextgap>gaptarget){
        # found it
        yaxisfinalgap <- currentgap
        gapfound <- TRUE
      }else{
        # need to go higher
        if(gapi==3){
          gapi <- 1
          gapscale <- gapscale/10
        }else{
          gapi <- gapi + 1
        }
      }
    }else{
      # need to go lower
      if(gapi==1){
        gapi <- 3
        gapscale <- gapscale*10
      }else{
        gapi <- gapi - 1
      }
    }
  }
  
  if(!is.null(yaxisfinalgap)) {
    axis(side = 2,
         at=seq(0,ylimData[2],yaxisfinalgap),
         las=2,
         col=signalColour,
         col.ticks=signalColour,
         col.axis=signalColour)
  }
  if(nrow(res_bd)>0){
    # if there is at least one segment to draw, draw a line
    currentSegment <- c(res_bd[1,"start"],res_bd[1,"end"])
    lines(as.numeric(currentSegment)/1e6,
          rep(res_bd$signal[1],2),
          lwd=lwd,
          col=signalColour)
    if(nrow(res_bd)>1){
      for(i in 2:nrow(res_bd)){
        # i <- 2
        previousSegment <- currentSegment
        currentSegment <- c(res_bd[i,"start"],res_bd[i,"end"])
        lines(rep(as.numeric(previousSegment[2]),2)/1e6,
              c(res_bd$signal[i-1],res_bd$signal[i]),
              lwd=lwd,
              col=signalColour)
        lines(as.numeric(currentSegment)/1e6,
              rep(res_bd$signal[i],2),
              lwd=lwd,
              col=signalColour)
      }
    }
  }
  
  # now plot the highlight regions if any
  if(nhighlightregions>0){
    for(i in 1:nhighlightregions){
      # i <- 1
      n <- names(highlightRegions)[i]
      # highlight region range
      rbottom <- ylimData[2]+highlightRegionHeightData*(i-1)
      rtop <- ylimData[2]+highlightRegionHeightData*i
      
      # draw some divisory dotted line
      abline(h=rbottom,lty=3,col="darkgrey")
      
      # write label on the right
      text(x = xtextpos,y=(rbottom+rtop)/2,
           labels=n,adj=0,xpd=TRUE,cex=cexlabels)
      
      # I need to find out whether we have any region to plot 
      regionBed <- highlightRegions[[n]][highlightRegions[[n]]$chr==pchr,,drop=F]
      if(nrow(regionBed)>0){
        selection <- !(regionBed$start > pend | regionBed$end < pstart)
        regionBed <- regionBed[selection,,drop=F]
      }
      if(nrow(regionBed)>0){
        for (j in 1:nrow(regionBed)){
          # j <- 1
          rect(xleft = regionBed[j,"start"]/1e6,
               ybottom = rbottom,
               xright = regionBed[j,"end"]/1e6,
               ytop = rtop,
               col = highlightRegionsColours[[n]],
               border = NA)
        }
      }
      
    }
    # more dotted lines if genes are plotted on top
    if(!is.null(levelAssignmentList)) abline(h=rtop,lty=3,col="darkgrey")
  }
  
  # now plot genes if any
  if(!is.null(levelAssignmentList)) {
    highlightregionsTop <- ylimData[2]+highlightRegionHeightData*nhighlightregions
    for(i in 1:nrow(genetable)){
      # i <- 1
      xleft <- max(genetable$start[i],pstart)
      xright <- min(genetable$end[i],pend)
      textxpos <- (xleft+xright)/2
      ypos <- highlightregionsTop+genesGap+(levelAssignmentList[[as.character(i)]]-1)*highlightRegionHeightData
      lines(x=c(xleft,xright)/1e6,y=c(ypos,ypos),lwd=3)
      text(genetable$genename[i],cex = 0.5,y = ypos+0.5*highlightRegionHeightData,x=textxpos/1e6,adj=0.5,xpd=T)
    }
  }
  
  # now plot the highlight positions if any
  if(nhighlightpositions>0){
    for(i in 1:nhighlightpositions){
      # i <- 1
      n <- names(highlightPositions)[i]
      
      # I need to find out whether we have any position to plot 
      regionPositions <- highlightPositions[[n]][highlightPositions[[n]]$chr==pchr,,drop=F]
      if(nrow(regionPositions)>0){
        selection <- regionPositions$position <= pend & regionPositions$position >= pstart
        regionPositions <- regionPositions[selection,,drop=F]
      }
      if(nrow(regionPositions)>0){
        for (j in 1:nrow(regionPositions)){
          # j <- 1
          abline(v=regionPositions[j,"position"]/1e6,
                 col=highlightPositionsColours[[n]],
                 lwd=lwd)
          if("text" %in% colnames(regionPositions)){
            text(x=(regionPositions[j,"position"]-2*xtextgap)/1e6,
                 y=sum(ylimData)/2,labels=regionPositions[j,"text"],
                 col=highlightPositionsColours[[n]],
                 srt=90,adj=0.5)
          }
        }
      }
    }
  }
  
  # I can only add the second signal after everything else has been plotted
  if(!is.null(res_bd2)){
    signalMin <- min(0,min(res_bd2$signal))
    signalMax <- max(1,max(res_bd2$signal))
    ydatagap <- 0.05*(signalMax-signalMin)
    ydatagapTop <- 0.15*(signalMax-signalMin)
    ylimData <- c(signalMin-ydatagap,signalMax+ydatagapTop)
    ySize <- ylimData[2]-ylimData[1]
    ylim <- c(ylimData[1],ylimData[2])
    par(fig=c(0,1,mbottom/pheight,(mbottom+dataheight)/pheight),
        new=TRUE,
        mai=c(0,mleft,0,mright),mgp=c(2.5,0.9,0))
    plot(NA,
         xlim=c(pstart,pend)/1e6,
         bty="n",
         ylim=ylim,
         main="",
         ylab="",
         las=1,
         xaxs="i",
         yaxs="i",
         yaxt="n",
         xaxt="n",
         xlab="")
    # ylabel
    xylabelpos <- (pend + 0.1*(pend-pstart))/1e6
    text(x = xylabelpos,
         y = (ylimData[1]+ylimData[2])/2,
         labels=ylabel2,srt=90,adj=0.5,xpd=T,col=signalColour2)
    lines(x=c(pend,pend)/1e6,y=c(ylimData[1],ylimData[2]),col="black",xpd=T)
    # ok figure out the axis from 0 to something just below ylimData[2]
    yaxisgaps <- c(1,2,5)
    gapi <- 1
    gapscale <- 1
    gapfound <- FALSE
    gaptarget <- ylimData[2]/2
    yaxisfinalgap <- NULL
    while (!gapfound) {
      currentgap <- yaxisgaps[gapi]/gapscale
      # check if we are below gaptarget
      if(currentgap<=gaptarget){
        # if the next one up is greater than target we are done
        if(gapi==3){
          nextgap <- yaxisgaps[1]/gapscale*10
        }else{
          nextgap <- yaxisgaps[gapi+1]/gapscale
        }
        if(nextgap>gaptarget){
          # found it
          yaxisfinalgap <- currentgap
          gapfound <- TRUE
        }else{
          # need to go higher
          if(gapi==3){
            gapi <- 1
            gapscale <- gapscale/10
          }else{
            gapi <- gapi + 1
          }
        }
      }else{
        # need to go lower
        if(gapi==1){
          gapi <- 3
          gapscale <- gapscale*10
        }else{
          gapi <- gapi - 1
        }
      }
    }
    
    if(!is.null(yaxisfinalgap)) {
      axis(side = 4,
           at=seq(0,ylimData[2],yaxisfinalgap),
           las=2,
           col=signalColour2,
           col.ticks=signalColour2,
           col.axis=signalColour2)
    }
    if(nrow(res_bd2)>0){
      # if there is at least one segment to draw, draw a line
      currentSegment <- c(res_bd2[1,"start"],res_bd2[1,"end"])
      lines(as.numeric(currentSegment)/1e6,
            rep(res_bd2$signal[1],2),
            lwd=lwd,
            col=signalColour2)
      if(nrow(res_bd2)>1){
        for(i in 2:nrow(res_bd2)){
          # i <- 2
          previousSegment <- currentSegment
          currentSegment <- c(res_bd2[i,"start"],res_bd2[i,"end"])
          lines(rep(as.numeric(previousSegment[2]),2)/1e6,
                c(res_bd2$signal[i-1],res_bd2$signal[i]),
                lwd=lwd,
                col=signalColour2)
          lines(as.numeric(currentSegment)/1e6,
                rep(res_bd2$signal[i],2),
                lwd=lwd,
                col=signalColour2)
        }
      }
    }
  }
  if(!is.null(fileout)) dev.off()
  
}




# positions need columns: chr, position, id, optionally class
# regions need columns: chr, start, end, id, optionally class
# class can be omitted, and a single class will be added
# regions cannot overlap as each position will be assigned to at most one region
# this function will automatically sort the positions and regions
# for faster assignment of positions to regions
# positions <- data.frame(chr = c(1,1,1,1,1),position = c(100,200,300,310,500),id=paste0("p",c(1,2,3,4,5)),class=paste0("pr",c(1,2,2,2,1)),stringsAsFactors = F)
# bed_table <- data.frame(chr = c(1,1,1,1),start = c(95,195,295,350),end = c(110,210,320,400),id=paste0("r",c(1,2,3,4)),class=c("cr1","cr1","cr2","cr2"),stringsAsFactors = F)

# positions <- data.frame(chr = c(1,1,1,1,1,1,1,1),position = c(100,200,300,400,500,600,700,800),
#                         id=paste0("y",c(1,2,3,4,5,6,7,8)),class=c("C","B","C","A","A","B","B","C"),stringsAsFactors = F)
# bed_table <- data.frame(chr = c(1,1,1,1),start = c(95,206,450,506),end = c(205,405,505,705),
#                         id=paste0("x",c(1,2,3,4)),class=c("g2","g1","g1","g2"),stringsAsFactors = F)

intersectPositionsAndBedRegions_nonOverlapping <- function(positions,
                                                           bed_table){
  # check column requirements
  # check required columns
  requiredcolumns <- c("chr","position","id")
  if(!all(requiredcolumns %in% colnames(positions))){
    missingcolumns <- setdiff(requiredcolumns,colnames(positions))
    message("[error intersectPositionsAndBedRegions_nonOverlapping] positions table missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  requiredcolumns <- c("chr","start","end","id")
  if(!all(requiredcolumns %in% colnames(bed_table))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table))
    message("[error intersectPositionsAndBedRegions_nonOverlapping] bed_table missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  
  # check which chromosomes have overlap if any
  overlapChroms <- checkBedRegionsOverlap(bed_table)
  if(!is.null(overlapChroms)){
    message("[error intersectPositionsAndBedRegions_nonOverlapping] regions should not overlap, ",
            "you can break them down using the function breakDownOverlappingBedRegions, ",
            "or you can run intersectPositionsAndBedRegions.")
    return(NULL)
  }
  # need them to be ordered
  positions <- sortPositions(positions)
  bed_table <- sortBed(bed_table)
  chroms <- sortChroms(union(positions$chr,bed_table$chr))
  # add single class if missing
  if(!("class" %in% colnames(positions))){
    positions$class <- "anyPosition"
  }
  pclasses <- as.character(unique(positions$class))
  if(!("class" %in% colnames(bed_table))){
    bed_table$class <- "anyRegion"
  }
  rclasses <- as.character(unique(bed_table$class))
  # I can obtain a count summary and also add a new column to positions with
  # the annotation of which class each position belongs to
  # countsTable <- rep(0,length(classes)+1)
  # names(countsTable) <- c(classes,"notInRegions")
  # countsTable <- data.frame(matrix(0,nrow = length(chroms),
  #                                  ncol = length(rclasses)+1,
  #                                  dimnames = list(chroms,c(rclasses,"notInRegions"))),
  #                           stringsAsFactors = F,
  #                           check.names = F)
  # classRegionSizes <- rep(0,length(rclasses)+1)
  # names(classRegionSizes) <- c(rclasses,"notInRegions")
  # classAnnotation <- NULL
  
  # now if I want to annotated both positions and bed regions, I need to be
  # able to add bed regions classes and ids to the positions table and
  # positions ids and classes to the bed table. In the case of the bed table,
  # more than one id of positions and more than one class are possible, so
  # lists are a good idea. Use ids for quick access.
  # perhaps I should do the classes at the end from the ids in one go.
  
  # index
  positions$id <- as.character(positions$id)
  bed_table$id <- as.character(bed_table$id)
  rownames(positions) <- positions$id
  rownames(bed_table) <- bed_table$id
  
  # annotation lists
  positions_RegionIDs <- list()
  regions_PositionIDs <- list()
  
  for (chrom in chroms) {
    # chrom <- chroms[1]
    positions_chrom <- positions[positions$chr==chrom,,drop=F]
    regions_chrom <- bed_table[bed_table$chr==chrom,,drop=F]
    # if there are no positions there is nothing to do in this chromosome
    if(nrow(positions_chrom)>0){
      # if there are no regions in this chromosome, all the positions are
      # classified as notInRegions
      # if(nrow(regions_chrom)==0){
      #   countsTableOld[chrom,"notInRegions"] <- nrow(positions_chrom)
      #   classAnnotation <- c(classAnnotation,rep("notInRegions",nrow(positions_chrom)))
      # }else{
      if(nrow(regions_chrom)>0){
        # in this chromosome we have both regions and mutations, ordered by position
        # now we start the core algorithm
        currentPosition <- 1
        currentRegion <- 1
        while (currentPosition<=nrow(positions_chrom) | currentRegion<=nrow(regions_chrom)) {
          # first of all, check that we have not finished with the positions or the regions
          if(currentPosition>nrow(positions_chrom)){
            # no more positions to check, so all remaining regions must be empty
            # just skip to the end
            currentRegion <- nrow(regions_chrom) + 1
          }else if(currentRegion>nrow(regions_chrom)){
            # no more regions to check, so if there are any other positions then they 
            # should be classified as notInRegions
            if(currentPosition<=nrow(positions_chrom)){
              nPositionsRemaining <- nrow(positions_chrom) - currentPosition + 1
              # add them
              # countsTableOld[chrom,"notInRegions"] <- countsTableOld[chrom,"notInRegions"] + nPositionsRemaining
              # classAnnotation <- c(classAnnotation,rep("notInRegions",nPositionsRemaining))
              # skip to the end
              currentPosition <- nrow(positions_chrom) + 1
            }
          }else{
            # we still have positions and regions left to check
            if(positions_chrom[currentPosition,"position"]<regions_chrom[currentRegion,"start"]){
              # the current position happens before the start of the current region, so it is notInRegions
              # add it
              # countsTableOld[chrom,"notInRegions"] <- countsTableOld[chrom,"notInRegions"] + 1
              # classAnnotation <- c(classAnnotation,"notInRegions")
              # move one position forward
              currentPosition <- currentPosition + 1
            }else if(positions_chrom[currentPosition,"position"]<=regions_chrom[currentRegion,"end"]){
              # if you get here, the current position is after start, but within the current region
              # add it to the class of the region
              # regionClass <- as.character(regions_chrom[currentRegion,"class"])
              # countsTableOld[chrom,regionClass] <- countsTableOld[chrom,regionClass] + 1
              # classAnnotation <- c(classAnnotation,regionClass)
              # annotate the position and the region
              positions_RegionIDs[[positions_chrom[currentPosition,"id"]]] <- c(positions_RegionIDs[[positions_chrom[currentPosition,"id"]]],regions_chrom[currentRegion,"id"])
              regions_PositionIDs[[regions_chrom[currentRegion,"id"]]] <- c(regions_PositionIDs[[regions_chrom[currentRegion,"id"]]],positions_chrom[currentPosition,"id"])
              # move one position forward
              currentPosition <- currentPosition + 1
            }else{
              # if you get here, then the current position must be after the end of the current region
              # which means we are done with this region and we move forward
              currentRegion <- currentRegion + 1
            }
          }
        }
      }
      
    }
  }
  # Annotation of positions and regions
  
  # 1. how many positions are associated with each region? must be 0 or n 
  # 2. how many positions for each class are associated with each region? either one or multiple classes or a noPositions class
  # 3. how many positions for each class are associated with each region class or with the noMatch class
  # 4. how many positions are associated with each region class or with the noMatch class (aggregate of 3.) 
  # 5. how many regions are associated with each position? must be 0 or 1 because of no overlapping regions
  # 6. how many regions for each class are associated with each position? either one class or a noMatch class
  # 7. how many regions for each class are associated with each position class or with the noPositions class
  # 8. how many regions are associated with each position class or with the noPositions class (aggregate of 7.)
  
  res_stats <- intersectionStatsComplete(idMap1to2 = positions_RegionIDs,
                                         idMap2to1 = regions_PositionIDs,
                                         idclassmap1 = positions,
                                         idclassmap2 = bed_table)
  
  # annotate
  # res_stats <- intersectionStats(idMap = positions_RegionIDs,
  #                                idclassmap1 = positions,
  #                                idclassmap2 = bed_table)
  # positions <- res_stats$idclassmap1
  # countsTable_positionsInRegionClasses <- res_stats$countsTable_classes1_in_classes2
  # countsTable_positionsInRegionClasses_total <- res_stats$countsTable_total1_in_classes2
  # countsTable_positionsInEachRegion <- res_stats$countsTable_classes1_in_id2
  # countsTable_positionsInEachRegion_total <- res_stats$countsTable_total1_in_id2
  # totalPostionsInAnyRegion <- res_stats$totalId1matchingAnyId2
  
  # if(length(positions_RegionIDs)>0){
  #   idlist <- unlist(positions_RegionIDs)
  #   positions[names(idlist),"regionIdAnnotation"] <- idlist
  #   positions[names(idlist),"regionClassAnnotation"] <- bed_table[idlist,"class"]
  #   positions[is.na(positions$regionClassAnnotation),"regionClassAnnotation"] <- "notInRegions"
  # }else{
  #   positions$regionIdAnnotation <- NA
  #   positions$regionClassAnnotation <- rep("notInRegions",nrow(positions))
  # }
  # # set up the counts table
  # countsTable_positions <- data.frame(matrix(0,nrow = length(pclasses),
  #                                            ncol = length(rclasses)+1,
  #                                            dimnames = list(pclasses,c(rclasses,"notInRegions"))),
  #                                     stringsAsFactors = F,
  #                                     check.names = F)
  # # Populate the countsTable
  # for(pi in pclasses){
  #   # pi <- pclasses[1]
  #   countsTable_tmp <- data.frame(as.list(table(positions$regionClassAnnotation[positions$class==pi])),
  #                                 check.names = F,stringsAsFactors = F,row.names = pi)
  #   countsTable_positions[rownames(countsTable_tmp),colnames(countsTable_tmp)] <- countsTable_tmp
  # }
  # countsTable_positions_total <- apply(countsTable_positions,2,sum)
  
  # add positions annotation to regions
  # 1. how many positions are associated with each region? must be 0 or n 
  # 2. how many positions for each class are associated with each region? either one or multiple classes or a noPositions class
  # 3. how many regions for each class are associated with each position class or with the noPositions class
  # 4. how many regions are associated with each position class or with the noPositions class (aggregate of 3.)
  
  # res_stats <- intersectionStats(idMap = regions_PositionIDs,
  #                                idclassmap1 = bed_table,
  #                                idclassmap2 = positions)
  # bed_table <- res_stats$idclassmap1
  # countsTable_regionsAtPositionClasses <- res_stats$countsTable_classes
  # countsTable_regionsAtPositionClasses_total <- res_stats$countsTable_classes_total
  # countsTable_positionsInEachRegion <- res_stats$countsTable_idclassmap1
  # countsTable_positionsInEachRegion_total <- res_stats$countsTable_idclassmap1_total
  # totalRegionsAtAnyPosition <- res_stats$totalId1matchingAnyId2
  
  # if(length(regions_PositionIDs)>0){
  #   # idlist <- unlist(regions_PositionIDs)
  #   res_annotColumns <- idlist_to_columns(idList = regions_PositionIDs,idmap = positions)
  #   bed_table[row.names(res_annotColumns),"positionIdAnnotation"] <- res_annotColumns$aggreatedIds
  #   bed_table[row.names(res_annotColumns),"positionClassAnnotation"] <- res_annotColumns$aggreatedMappedIds
  #   bed_table[is.na(bed_table$positionClassAnnotation),"positionClassAnnotation"] <- "noPositions"
  #   
  # }else{
  #   bed_table$positionIdAnnotation <- NA
  #   bed_table$positionClassAnnotation <- rep("noPositions",nrow(bed_table))
  # 
  # }
  # # set up the counts table
  # countsTable_regions <- data.frame(matrix(0,nrow = length(rclasses),
  #                                          ncol = length(pclasses)+1,
  #                                          dimnames = list(rclasses,c(pclasses,"noPositions"))),
  #                                   stringsAsFactors = F,
  #                                   check.names = F)
  # countsTable_regions_total <- data.frame(matrix(0,nrow = 1,
  #                                                ncol = length(pclasses)+1,
  #                                                dimnames = list("1",c(pclasses,"noPositions"))),
  #                                         stringsAsFactors = F,
  #                                         check.names = F)
  # # Populate the countsTable
  # for(ri in rclasses){
  #   # ri <- rclasses[1]
  #   select_class_withPositions <- bed_table$class==ri & !bed_table$positionClassAnnotation=="noPositions"
  #   select_class_noPositions <- bed_table$class==ri & bed_table$positionClassAnnotation=="noPositions"
  #   rids <- bed_table$id[select_class_withPositions]
  #   pids <- unique(unlist(regions_PositionIDs[rids]))
  #   tpc <- table(positions[pids,"class"])
  #   countsTable_regions[ri,names(tpc)] <- as.vector(tpc)
  #   # check regions with no positions for this class
  #   countsTable_regions[ri,"noPositions"] <- sum(select_class_noPositions)
  # }
  # # countsTable_regions_total <- apply(countsTable_regions,2,sum)
  # select_class_withPositions <- !bed_table$positionClassAnnotation=="noPositions"
  # select_class_noPositions <- bed_table$positionClassAnnotation=="noPositions"
  # rids <- bed_table$id[select_class_withPositions]
  # pids <- unique(unlist(regions_PositionIDs[rids]))
  # tpc <- table(positions[pids,"class"])
  # countsTable_regions_total["1",names(tpc)] <- as.vector(tpc)
  # # check regions with no positions for this class
  # countsTable_regions_total["1","noPositions"] <- sum(select_class_noPositions)
  
  # collect results in the return object
  returnObj <- list()
  
  returnObj$annotatedPositions <- res_stats$idclassmap1_updated
  returnObj$annotatedBedRegions <- res_stats$idclassmap2_updated
  
  returnObj$totalPostionsInAnyRegion <- totalPostionsInAnyRegion
  returnObj$totalRegionsAtAnyPosition <- totalRegionsAtAnyPosition
  
  returnObj$positions_RegionIDs <- positions_RegionIDs
  returnObj$regions_PositionIDs <- regions_PositionIDs
  
  returnObj$countsTable_positionsInEachRegion <- res_stats$countsTable_classes1_in_id2
  returnObj$countsTable_positionsInEachRegion_total <- res_stats$countsTable_total1_in_id2
  returnObj$countsTable_positionsInRegionClasses <- res_stats$countsTable_classes1_in_classes2
  returnObj$countsTable_positionsInRegionClasses_total <- res_stats$countsTable_total1_in_classes2
  
  returnObj$countsTable_regionsAtEachPosition <- res_stats$countsTable_classes2_in_id1
  returnObj$countsTable_regionsAtEachPosition_total <- res_stats$countsTable_total2_in_id1
  returnObj$countsTable_regionsAtPositionClasses <- res_stats$countsTable_classes2_in_classes1
  returnObj$countsTable_regionsAtPositionClasses_total <- res_stats$countsTable_total2_in_classes1
  
  return(returnObj)
}

# idlist_to_columns <- function(idList,
#                               idmap){
#   aggreatedIds <- unlist(lapply(idList,function(x) paste(x,collapse = ";")))
#   aggreatedMappedIds <- unlist(lapply(idList,function(x) {
#     mapped <- idmap[x,"class"]
#     tmapped <- table(mapped)
#     returnstring <- NULL
#     for(i in 1:length(tmapped)) returnstring <- c(returnstring,paste(tmapped[i],names(tmapped)[i],sep = ":"))
#     return(paste(returnstring,collapse = ";"))
#   }))
#   return(data.frame(aggreatedIds,aggreatedMappedIds,stringsAsFactors = F))
# }

# idlist_to_classCounts <- function(idList,
#                                   idmap){
#   aggreatedIds <- unlist(lapply(idList,function(x) paste(x,collapse = ";")))
#   aggreatedMappedIds <- unlist(lapply(idList,function(x) {
#     mapped <- idmap[x,"class"]
#     tmapped <- table(mapped)
#     
#     return(data.frame(matrix(as.vector(tmapped),nrow = 1,ncol = length(tmapped),dimnames = list("1",names(tmapped)))))
#   }))
#   return(data.frame(aggreatedIds,aggreatedMappedIds,stringsAsFactors = F))
# }



#' Intersection statistics
#'
#' This is a low level function, the engine that calculates the overlap between
#' any two groups of genomics entities. Entities (e.g. positions and regions),
#' have ids and classes, so that ids are unique for a given group of entities,
#' while classes may be repeated, so that subgroups of entities can belong to
#' the same class. The intersectionStats function assumes that, whatever the
#' entities, the intersection of two groups of entities has been computed and
#' is represented by an idMap, which maps the ids from the first group to the
#' ids of the second group. The mapping is bidirectional, i.e. if an entity
#' in the first group intersects an entity in the second group then the converse
#' is also true. To copmute the statistics for the classes, the mapping from id
#' to class for each group needs to be provided (idclassmap1 and idclassmap2).
#' To obtain all statistics, one needs to run intersectionStats twice, swapping
#' the two groups, so that in the second run idMap maps from the ids of the second
#' group to the ids of the first group and idclassmap1 and idclassmap2 are swapped.
#' 
#' 
#' @param idMap list object where the names are ids of the first group, and each id maps to a vector of ids of the second group
#' @param idclassmap1 data frame with required columns id and class, id must be unique
#' @param idclassmap2 data frame with required columns id and class, id must be unique
#' @return object with details intersection statistics
#' @export
intersectionStats <- function(idMap,
                              idclassmap1,
                              idclassmap2){
  
  # index
  idclassmap1$id <- as.character(idclassmap1$id)
  idclassmap2$id <- as.character(idclassmap2$id)
  rownames(idclassmap1) <- idclassmap1$id
  rownames(idclassmap2) <- idclassmap2$id

  # get classes
  classes1 <- unique(idclassmap1$class)
  classes2 <- unique(idclassmap2$class)
  
  # set up the counts tables
  countsTable_classes1_in_classes2 <- data.frame(matrix(0,nrow = length(classes1),
                                                        ncol = length(classes2)+1,
                                                        dimnames = list(classes1,c(classes2,"noMatch"))),
                                                 stringsAsFactors = F,
                                                 check.names = F)
  countsTable_total1_in_classes2 <- data.frame(matrix(0,nrow = 1,
                                                      ncol = length(classes2)+1,
                                                      dimnames = list("1",c(classes2,"noMatch"))),
                                               stringsAsFactors = F,
                                               check.names = F)
  countsTable_classes1_in_id2 <- data.frame(matrix(0,nrow = length(classes1),
                                                   ncol = nrow(idclassmap2)+1,
                                                   dimnames = list(classes1,c(rownames(idclassmap2),"noMatch"))),
                                            stringsAsFactors = F,
                                            check.names = F)
  # annotate idclassmap1
  if(length(idMap)>0){
    # res_annotColumns <- idlist_to_columns(idMap = idMap,idmap = idclassmap2)
    aggreatedIds <- unlist(lapply(idMap,function(x) paste(x,collapse = ";")))
    # aggreatedMappedIds <- unlist(lapply(idMap,function(x) {
    #   mapped <- idclassmap2[x,"class"]
    #   tmapped <- table(mapped)
    #   returnstring <- NULL
    #   for(i in 1:length(tmapped)) returnstring <- c(returnstring,paste(tmapped[i],names(tmapped)[i],sep = ":"))
    #   return(paste(returnstring,collapse = ";"))
    # }))
    aggreatedMappedIds <- NULL
    for(id in names(idMap)){
      # id <- names(idMap)[1]
      x <- idMap[[id]]
      mapped <- idclassmap2[x,"class"]
      tmapped <- table(mapped)
      returnstring <- NULL
      for(i in 1:length(tmapped)) returnstring <- c(returnstring,paste(tmapped[i],names(tmapped)[i],sep = ":"))
      aggreatedMappedIds <- c(aggreatedMappedIds,paste(returnstring,collapse = ";"))
      # countsTable_classes2_in_id1[names(tmapped),id] <- as.vector(tmapped)
    }
    res_annotColumns <- data.frame(aggreatedIds,aggreatedMappedIds,stringsAsFactors = F)
    
    idclassmap1[row.names(res_annotColumns),"idAnnotation"] <- res_annotColumns$aggreatedIds
    idclassmap1[row.names(res_annotColumns),"classAnnotation"] <- res_annotColumns$aggreatedMappedIds
    idclassmap1[is.na(idclassmap1$idAnnotation),"classAnnotation"] <- "noMatch"
    
  }else{
    idclassmap1$idAnnotation <- NA
    idclassmap1$classAnnotation <- rep("noMatch",nrow(idclassmap1))
    
  }
  
  # Populate the countsTable
  for(ci in classes1){
    # ci <- classes1[1]
    select_class_withMatch <- idclassmap1$class==ci & !idclassmap1$classAnnotation=="noMatch"
    select_class_noMatch <- idclassmap1$class==ci & idclassmap1$classAnnotation=="noMatch"
    if(sum(select_class_withMatch)>0){
      rids <- idclassmap1$id[select_class_withMatch]
      
      # pids <- unique(unlist(idMap[rids]))
      # tpc <- table(idclassmap2[pids,"class"])
      
      res_ids <- unlist(idMap[rids])
      res_classes <- unlist(lapply(idMap[rids], function(x){
        unique(idclassmap2[x,"class"])
      }))
      
      tpi <- table(res_ids)
      tpc <- table(res_classes)
      countsTable_classes1_in_id2[ci,names(tpi)] <- as.vector(tpi)
      countsTable_classes1_in_classes2[ci,names(tpc)] <- as.vector(tpc)
    }
    if(sum(select_class_noMatch)>0){
      # check regions with no positions for this class
      countsTable_classes1_in_classes2[ci,"noMatch"] <- sum(select_class_noMatch)
      countsTable_classes1_in_id2[ci,"noMatch"] <- sum(select_class_noMatch)
    }
  }
  # countsTable_regions_total <- apply(countsTable_regions,2,sum)
  select_class_withMatch <- !idclassmap1$classAnnotation=="noMatch"
  select_class_noMatch <- idclassmap1$classAnnotation=="noMatch"
  if(sum(select_class_withMatch)>0){
    rids <- idclassmap1$id[select_class_withMatch]
    
    # pids <- unique(unlist(idMap[rids]))
    # tpc <- table(idclassmap2[pids,"class"])
    
    res_classes <- unlist(lapply(idMap[rids], function(x){
      unique(idclassmap2[x,"class"])
    }))
    tpc <- table(res_classes)

    countsTable_total1_in_classes2["1",names(tpc)] <- as.vector(tpc)
  }
  if(sum(select_class_noMatch)>0){
    # check regions with no positions for this class
    countsTable_total1_in_classes2["1","noMatch"] <- sum(select_class_noMatch)
  }
  
  # return the results
  returnObj <- list()
  returnObj$idclassmap1_updated <- idclassmap1
  returnObj$countsTable_classes1_in_classes2 <- countsTable_classes1_in_classes2
  returnObj$countsTable_total1_in_classes2 <- countsTable_total1_in_classes2
  returnObj$countsTable_classes1_in_id2 <- countsTable_classes1_in_id2
  returnObj$countsTable_total1_in_id2 <- apply(countsTable_classes1_in_id2,2,sum)
  returnObj$totalId1matchingAnyId2 <- sum(!is.na(idclassmap1$idAnnotation))
  return(returnObj)
}

# running both directions and annotating at once
intersectionStatsComplete <- function(idMap1to2,
                                      idMap2to1,
                                      idclassmap1,
                                      idclassmap2){
  # running both directions
  res1 <- intersectionStats(idMap = idMap1to2,
                            idclassmap1 = idclassmap1,
                            idclassmap2 = idclassmap2)
  res2 <- intersectionStats(idMap = idMap2to1,
                            idclassmap1 = idclassmap2,
                            idclassmap2 = idclassmap1)
  
  # combine the objects and return
  returnObj <- res1
  returnObj$idclassmap2_updated <- res2$idclassmap1_updated
  returnObj$countsTable_classes2_in_classes1 <- res2$countsTable_classes1_in_classes2
  returnObj$countsTable_total2_in_classes1 <- res2$countsTable_total1_in_classes2
  returnObj$countsTable_classes2_in_id1 <- res2$countsTable_classes1_in_id2
  returnObj$countsTable_total2_in_id1 <- res2$countsTable_total1_in_id2
  returnObj$totalId2matchingAnyId1 <- res2$totalId1matchingAnyId2
  return(returnObj)
}

