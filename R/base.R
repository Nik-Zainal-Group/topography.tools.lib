
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


