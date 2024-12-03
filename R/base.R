
#' @importFrom foreach %dopar%
NULL

#' @importFrom doRNG %dorng%
NULL

getEntitiesType <- function(entities){
  colnamesPositions <- c("chr","position")
  colnamesBedRegions <- c("chr","start","end")
  isPositions <- all(colnamesPositions %in% colnames(entities))
  isBedRegions <- all(colnamesBedRegions %in% colnames(entities))
  if(isPositions & isBedRegions){
    return("ambiguous")
  }else if(isPositions){
    return("positions")
  }else if(isBedRegions){
    return("bedRegions")
  }else{
    return("unknown")
  }
}

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
            assignedIds <- c(assignedIds,as.character(id))
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
#' @param bed_table data frame with required columns: chr, start, end, id, and optionally signal
#' @param aggregateSignalMode when merging bed regions that have a signal column,
#' options for aggregating signal are sum, mean, or weightedmean (weighted w.r.t size)
#' @return updated bed_table with merged adjacent regions
#' @export
mergeAdjacentBedRegions <- function(bed_table,
                                    aggregateSignalMode="sum"){
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
  
  acceptedSignalModes <- c("sum","mean","weightedmean")
  # check aggregateSignalMode
  if(!aggregateSignalMode %in% acceptedSignalModes){
    message("[error mergeAdjacentBedRegions] invalid aggregateSignalMode, please use on of: ",paste(acceptedSignalModes,collapse = ", "))
    return(NULL)
  }
  
  # add signal to the required columns so we know we need to keep it
  if("signal" %in% colnames(bed_table)) requiredcolumns <- c(requiredcolumns,"signal")
  
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
          newrow <- data.frame(id=paste(chrTable[mergeRows,"id"],collapse = ";"),
                               chr=chrom,
                               start=min(chrTable[mergeRows,"start"]),
                               end=max(chrTable[mergeRows,"end"]),
                               stringsAsFactors = F)
          
          # now add signal id necessary
          if("signal" %in% requiredcolumns){
            if(aggregateSignalMode=="sum"){
              newrow$signal <- sum(chrTable[mergeRows,"signal"])
            }else if(aggregateSignalMode=="mean"){
              newrow$signal <- mean(chrTable[mergeRows,"signal"])
            }else if(aggregateSignalMode=="weightedmean"){
              sizes <- chrTable[mergeRows,"end"] - chrTable[mergeRows,"start"] + 1
              newrow$signal <- as.vector(chrTable[mergeRows,"signal"] %*% sizes)/sum(sizes)
            }
          }
          
          # add to new table
          newTable <- rbind(newTable,newrow)
          
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


#' Plot signal of bed regions
#'
#' Given a table of bed regions with the signal column, plot the signal values
#' across a given region. Signal for segments not included in the given bed regions
#' will be set to zero. A second bed table can be specified, and its signal will be
#' plotted along with the signal of the first bed table.
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
#' @param proteinCodingOnly if TRUE then plot only protein coding genes 
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
                                proteinCodingOnly = TRUE,
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
    if(proteinCodingOnly) genetable <- genetable[genetable$genetype=="protein_coding",,drop=F]
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


#' Extend bed regions preserving no overlap
#'
#' Given a table of non-overlapping bed regions, extend the regions
#' in both directions by a given length. If two regions are next to each
#' other and are at a distance that is less than 2 times the extend length,
#' then they will be both extended by half such distance so that they become
#' adjacent regions. When extending at the beginning of a chromosome, 1 will be
#' the minimum position.
#' 
#' 
#' @param bed_table data frame with required columns: chr, start, end
#' @param extended length in number of bases by which each region should be
#' extended in both directions
#' @return updated bed_table
#' @export
extendBedRegionsWithNoOverlap <- function(bed_table,
                                          extended){
  # check no overlap before extension
  overlapChroms <- checkBedRegionsOverlap(bed_table = bed_table)
  if(!is.null(overlapChroms)){
    message("[error extendBedRegionsWithNoOverlap] regions in bed_table should not overlap, ",
            "you can break them down using the function breakDownOverlappingBedRegions.")
    return(NULL)
  }
  
  tmp_bed_table <- NULL
  if(extended>0){
    # extend each chrom separately so that extensions don't overlap
    tmp_chroms <- unique(bed_table$chr)
    for(tmp_chrom in tmp_chroms){
      tmp_bed_table_chrom <- bed_table[bed_table$chr==tmp_chrom,,drop=F]
      if(nrow(tmp_bed_table_chrom)>0){
        for(i in 1:nrow(tmp_bed_table_chrom)){
          if(i==1) {
            tmp_bed_table_chrom$start[i] <- max(tmp_bed_table_chrom$start[i]-extended,1)
          }
          if(i+1<=nrow(tmp_bed_table_chrom)){
            # we need to check that the end of this extended region does not
            # overlap with the start of the next region
            extended_diff <-  (tmp_bed_table_chrom$start[i+1] - extended) - (tmp_bed_table_chrom$end[i] + extended)
            if(extended_diff>0){
              tmp_bed_table_chrom$end[i] <- tmp_bed_table_chrom$end[i] + extended
              tmp_bed_table_chrom$start[i+1] <- tmp_bed_table_chrom$start[i+1] - extended
            }else{
              # find the place where to break the overlap, halfway between the region boundaries should do
              midpoint <- floor((tmp_bed_table_chrom$end[i] + tmp_bed_table_chrom$start[i+1])/2)
              tmp_bed_table_chrom$end[i] <- midpoint
              tmp_bed_table_chrom$start[i+1] <- midpoint + 1
            }
          }
          if(i==nrow(tmp_bed_table_chrom)){
            tmp_bed_table_chrom$end[i] <- tmp_bed_table_chrom$end[i]+extended
          }
        }
        tmp_bed_table <- rbind(tmp_bed_table,tmp_bed_table_chrom)
      }
    }
  }else{
    tmp_bed_table <- bed_table
  }
  return(tmp_bed_table)
}


#' Extend bed regions 
#'
#' Given a table of bed regions, extend the regions in both directions by a given
#' length. When extending at the beginning of a chromosome, 1 will be
#' the minimum position. This function will extend bed regions allowing bed regions
#' to overlap. Use function extendBedRegionsWithNoOverlap for extending non-overlapping
#' regions and preserving the non-overlap property.
#' 
#' 
#' @param bed_table data frame with required columns: chr, start, end
#' @param extended length in number of bases by which each region should be
#' extended in both directions
#' @return updated bed_table
#' @export
extendBedRegions <- function(bed_table,
                             extended){
  bed_table$end <- bed_table$end + extended
  tmpnewstart <- bed_table$start - extended
  bed_table$start <- sapply(tmpnewstart,function(x) {
    return(ifelse(x>0,x,1))
  })
  return(bed_table)
}

#' Compute inter-region distance of a set of bed regions
#'
#' This function returns both the left and right IRD of each given region
#' as well as the average IRD. Regions must be non-overlapping
#' 
#' @param bed_table data frame with required columns: chr, start, end 
#' @return left, right and average inter-region distance
#' @export
getIRD <- function(bed_table){
  # some checks
  requiredcolumns_pos <- c("chr","start","end")
  if(!all(requiredcolumns_pos %in% colnames(bed_table))){
    missingcolumns <- setdiff(requiredcolumns_pos,colnames(bed_table))
    message("[error getIRD] missing required columns in bed_table: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  
  # I need to check that there are no overlaps
  overlapChroms <- checkBedRegionsOverlap(bed_table = bed_table)
  if(!is.null(overlapChroms)){
    message("[error getIRD] regions in bed_table should not overlap, ",
            "you can break them down using the function breakDownOverlappingBedRegions.")
    return(NULL)
  }
  
  # compute the IRD, for each chromosome separately, sort first
  bed_table <- sortBed(bed_table = bed_table)
  chroms <- unique(bed_table$chr)
  new_bed_table <- NULL
  for(chrom in chroms){
    # chrom <- chroms[1]
    chr_table <- bed_table[bed_table$chr==chrom,,drop=F]
    if(nrow(chr_table)>0){
      # get ready
      chr_table$leftIRD <- NA
      chr_table$rightIRD <- NA
      chr_table$aveIRD <- NA
      # there is at least one row, so if it is only one row IRD=NA
      # and just add to final table, otherwise we calculate the IRD and fill the table
      if(nrow(chr_table)>1){
        IRD <- chr_table$start[2:nrow(chr_table)] - chr_table$end[1:(nrow(chr_table)-1)]
        chr_table[2:nrow(chr_table),"leftIRD"] <- IRD
        chr_table[1:(nrow(chr_table)-1),"rightIRD"] <- IRD
        chr_table[,"aveIRD"] <- apply(chr_table[,c("leftIRD","rightIRD")],1,mean,na.rm=T)
      }
      new_bed_table <- rbind(new_bed_table,chr_table)
    }
  }
  return(new_bed_table)
}



#' Compute distance of positions to the nearest bed region
#'
#' Given a set of positions and a set of bed regions, annotate the distance of
#' each position to the nearest bed region. Distance is 0 if the position is in
#' a region, and NA if the position is on a chromosome where there are no regions.
#' Regions need to be non-overlapping.
#' 
#' @param positions data frame with required columns: chr, position, id 
#' @param bed_table data frame with required columns: chr, start, end, id 
#' @return annotated positions with distance to nearest bed region
#' @export
distanceOfPositionToNearestBedRegion <- function(positions,
                                                 bed_table){
  
  # let's check which positions are inside a bed_table, so have distance 0
  res_int <- intersectPositionsAndBedRegions_nonOverlapping(positions = positions,
                                                            bed_table = bed_table,
                                                            computeStats = TRUE)
  annotatedPositions <- res_int$annotatedPositions
  annotatedPositions$classAnnotation <- NULL
  colnames(annotatedPositions)[which(colnames(annotatedPositions)=="idAnnotation")] <- "nearestRegion"
  annotatedPositions$distanceToNearestRegion[!is.na(annotatedPositions$nearestRegion)] <- 0
  
  # now let's find the distance of the mutations not in the regions
  leftoverPositions <- annotatedPositions[is.na(annotatedPositions$nearestRegion),,drop=F]
  
  if(nrow(leftoverPositions)>0){
    # index the positions and the regions
    rownames(positions) <- positions$id
    rownames(bed_table) <- bed_table$id
    
    # find the nearest region to each position by extending the regions
    # need to know how much to extend, can be tricky, so just extend max chrom 
    # length, so about 250 mil
    bed_table_extended <- extendBedRegionsWithNoOverlap(bed_table = bed_table,
                                                        extended = 250000000)
    res_int_ext <- intersectPositionsAndBedRegions_nonOverlapping(positions = leftoverPositions,
                                                                  bed_table = bed_table_extended,
                                                                  computeStats = TRUE)
    # if some positions still are not assigned, then there are no regions in the chromosome
    annotatedLeftover <- res_int_ext$annotatedPositions
    annotatedLeftover <- annotatedLeftover[!is.na(annotatedLeftover$idAnnotation),,drop=F]
    # now we can update
    if(nrow(annotatedLeftover)>0){
      annotatedPositions[annotatedLeftover$id,"nearestRegion"] <- annotatedLeftover$idAnnotation
      annotatedPositions[annotatedLeftover$id,"distanceToNearestRegion"] <- apply(abs(annotatedLeftover$position - bed_table[annotatedLeftover$idAnnotation,c("start","end")]),1,min)
      # recover sign
      idInvertSign <- annotatedLeftover$id[(annotatedLeftover$position - bed_table[annotatedLeftover$idAnnotation,"start"]) < 0]
      if(length(idInvertSign)>0){
        annotatedPositions[idInvertSign,"distanceToNearestRegion"] <- - annotatedPositions[idInvertSign,"distanceToNearestRegion"]
      }
    }
  }

  # return annotated positions
  return(annotatedPositions)
}
