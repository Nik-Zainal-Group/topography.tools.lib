

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
#' is also true. The mapping from id to class for each group needs to be provided
#' (idclassmap1 and idclassmap2), so that the statistics for the classes can be computed.
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
    aggreatedIds <- unlist(lapply(idMap,function(x) paste(x,collapse = ";")))
    aggreatedMappedIds <- NULL
    for(id in names(idMap)){
      # id <- names(idMap)[1]
      x <- idMap[[id]]
      mapped <- idclassmap2[x,"class"]
      tmapped <- table(mapped)
      returnstring <- NULL
      for(i in 1:length(tmapped)) returnstring <- c(returnstring,paste(tmapped[i],names(tmapped)[i],sep = ":"))
      aggreatedMappedIds <- c(aggreatedMappedIds,paste(returnstring,collapse = ";"))
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

#' Intersection statistics complete
#'
#' This is a low level function, the engine that calculates the overlap between
#' any two groups of genomics entities. Entities (e.g. positions and regions),
#' have ids and classes, so that ids are unique for a given group of entities,
#' while classes may be repeated, so that subgroups of entities can belong to
#' the same class. The intersectionStats function assumes that, whatever the
#' entities, the intersection of two groups of entities has been computed and
#' is represented by an idMap1to2, which maps the ids from the first group to the
#' ids of the second group. The mapping is bidirectional, i.e. if an entity
#' in the first group intersects an entity in the second group then the converse
#' is also true. The converse idMap2to1 can also be provided, or it will be inferred
#' automatically if left NULL.
#' The mapping from id to class for each group needs to be provided
#' (idclassmap1 and idclassmap2), so that the statistics for the classes can be computed.
#' In practice, intersectionStatsComplete calls intersectionStats twice, so that
#' statistics of intersection can be calculated from both direction, i.e. counting elements of 
#' of the first group that overlap the second group, or counting elements of the 
#' second group that overlap the first group.
#' 
#' 
#' @param idMap1to2 list object where the names are ids of the first group, and each id maps to a vector of ids of the second group
#' @param idMap2to1 list object where the names are ids of the second group, and each id maps to a vector of ids of the first group.
#' If NULL, this will be computed automatically from idMap1to2 
#' @param idclassmap1 data frame with required columns id and class, id must be unique
#' @param idclassmap2 data frame with required columns id and class, id must be unique
#' @return object with details intersection statistics
#' @export
intersectionStatsComplete <- function(idMap1to2,
                                      idMap2to1=NULL,
                                      idclassmap1,
                                      idclassmap2){
  
  if(is.null(idMap2to1)){
    idMap2to1 <- reverseIdMap(idMap = idMap1to2)
  }
  
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

reverseIdMap <- function(idMap){
  reverseMap <- list()
  for (id1 in names(idMap)) {
    # id1 <- names(idMap)[1]
    for(id2 in idMap[[id1]]){
      # id2 <- idMap[[id1]][1]
      id2 <- as.character(id2)
      reverseMap[[id2]] <- c(reverseMap[[id2]],id1)
    }
  }
  return(reverseMap)
}

mergeIdMaps <- function(idMap1,
                        idMap2){
  if(length(idMap1)==0){
    return(idMap2)
  }else if(length(idMap2)==0){
    return(idMap1)
  }else{
    for(id in names(idMap2)){
      idMap1[[id]] <- unique(c(idMap1[[id]],idMap2[[id]]))
    }
    return(idMap1)
  }
}

#' Intersect positions with non-overlapping bed regions
#'
#' Given a table of positions and a table of bed regions, find the positions that
#' are contained in the regions, and conversely the regions that overlap given positions.
#' If positions and/or bed regions have classes, then find how many positions for each 
#' class of positions are contained in each region or class of regions, and find how many regions for
#' each class of regions contain each position or each class of positions.
#' This function is restricted to only non-overlapping bed regions. For overlapping
#' bed regions, use the more general function intersectPositionsAndBedRegions.
#' 
#' @param positions data frame containing positions, with required columns chr, position, id and optionally class. Value in the id column must be unique
#' @param bed_table data frame containing bed regions, with required columns chr, start, end, id and optionally class. Value in the id column must be unique. Regions cannot overlap.
#' @param computeStats if FALSE, intersect stats will not be calculated and only the id maps mapping
#' position ids to region ids and viceversa will be returned. This is meant to save compute time when
#' intersectPositionsAndBedRegions_nonOverlapping is run inside intersectPositionsAndBedRegions.
#' @return object with details intersection statistics
#' @export
intersectPositionsAndBedRegions_nonOverlapping <- function(positions,
                                                           bed_table,
                                                           computeStats=TRUE){
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
  # pclasses <- as.character(unique(positions$class))
  if(!("class" %in% colnames(bed_table))){
    bed_table$class <- "anyRegion"
  }
  # rclasses <- as.character(unique(bed_table$class))
  
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
  idMapPositionsToRegions <- list()
  idMapRegionsToPositions <- list()
  
  for (chrom in chroms) {
    # chrom <- chroms[1]
    positions_chrom <- positions[positions$chr==chrom,,drop=F]
    regions_chrom <- bed_table[bed_table$chr==chrom,,drop=F]
    # if there are no positions there is nothing to do in this chromosome
    if(nrow(positions_chrom)>0){
      # if there are no regions in this chromosome, all the positions are
      # classified as noMatch
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
            # should be classified as noMatch
            if(currentPosition<=nrow(positions_chrom)){
              nPositionsRemaining <- nrow(positions_chrom) - currentPosition + 1
              # skip to the end
              currentPosition <- nrow(positions_chrom) + 1
            }
          }else{
            # we still have positions and regions left to check
            if(positions_chrom[currentPosition,"position"]<regions_chrom[currentRegion,"start"]){
              # the current position happens before the start of the current region, so it is noMatch
              # move one position forward
              currentPosition <- currentPosition + 1
            }else if(positions_chrom[currentPosition,"position"]<=regions_chrom[currentRegion,"end"]){
              # if you get here, the current position is after start, but within the current region
              # add it to the class of the region
              # annotate the position and the region
              idMapPositionsToRegions[[positions_chrom[currentPosition,"id"]]] <- c(idMapPositionsToRegions[[positions_chrom[currentPosition,"id"]]],regions_chrom[currentRegion,"id"])
              idMapRegionsToPositions[[regions_chrom[currentRegion,"id"]]] <- c(idMapRegionsToPositions[[regions_chrom[currentRegion,"id"]]],positions_chrom[currentPosition,"id"])
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
  
  if(computeStats){
    res_stats <- intersectionStatsComplete(idMap1to2 = idMapPositionsToRegions,
                                           idMap2to1 = idMapRegionsToPositions,
                                           idclassmap1 = positions,
                                           idclassmap2 = bed_table)
  }
  

  
  # collect results in the return object
  returnObj <- list()
  
  returnObj$idMapPositionsToRegions <- idMapPositionsToRegions
  returnObj$idMapRegionsToPositions <- idMapRegionsToPositions

  if(computeStats){
    returnObj$annotatedPositions <- res_stats$idclassmap1_updated
    returnObj$annotatedBedRegions <- res_stats$idclassmap2_updated
    
    returnObj$totalPostionsInAnyRegion <- res_stats$totalId1matchingAnyId2
    returnObj$totalRegionsAtAnyPosition <- res_stats$totalId2matchingAnyId1
    
    returnObj$countsTable_positionsInEachRegion <- res_stats$countsTable_classes1_in_id2
    returnObj$countsTable_positionsInEachRegion_total <- res_stats$countsTable_total1_in_id2
    returnObj$countsTable_positionsInRegionClasses <- res_stats$countsTable_classes1_in_classes2
    returnObj$countsTable_positionsInRegionClasses_total <- res_stats$countsTable_total1_in_classes2
    
    returnObj$countsTable_regionsAtEachPosition <- res_stats$countsTable_classes2_in_id1
    returnObj$countsTable_regionsAtEachPosition_total <- res_stats$countsTable_total2_in_id1
    returnObj$countsTable_regionsAtPositionClasses <- res_stats$countsTable_classes2_in_classes1
    returnObj$countsTable_regionsAtPositionClasses_total <- res_stats$countsTable_total2_in_classes1
  }
  return(returnObj)
}

#' Intersect positions with bed regions
#'
#' Given a table of positions and a table of bed regions, find the positions that
#' are contained in the regions, and conversely the regions that overlap given positions.
#' If positions and/or bed regions have classes, then find how many positions for each 
#' class of positions are contained in each region or class of regions, and find how many regions for
#' each class of regions contain each position or each class of positions.
#' This function allows for overlapping bed regions. In practice, if there are 
#' overlapping bed regions, the regions in bed_table will be assigned to a minimum number of
#' sets such that each set contains non-overlapping regions. The function intersectPositionsAndBedRegions_nonOverlapping
#' will then be used on the non-overlapping sets separately, and the results merged.
#' 
#' @param positions data frame containing positions, with required columns chr, position, id and optionally class. Value in the id column must be unique
#' @param bed_table data frame containing bed regions, with required columns chr, start, end, id and optionally class. Value in the id column must be unique
#' @return object with details intersection statistics
#' @export
intersectPositionsAndBedRegions <- function(positions,
                                            bed_table){
  # check column requirements
  # check required columns
  requiredcolumns <- c("chr","position","id")
  if(!all(requiredcolumns %in% colnames(positions))){
    missingcolumns <- setdiff(requiredcolumns,colnames(positions))
    message("[error intersectPositionsAndBedRegions] positions table missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  requiredcolumns <- c("chr","start","end","id")
  if(!all(requiredcolumns %in% colnames(bed_table))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table))
    message("[error intersectPositionsAndBedRegions] bed_table missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  
  # add single class if missing
  if(!("class" %in% colnames(positions))){
    positions$class <- "anyPosition"
  }
  if(!("class" %in% colnames(bed_table))){
    bed_table$class <- "anyRegion"
  }
  
  # check which chromosomes have overlap if any
  overlapChroms <- checkBedRegionsOverlap(bed_table)
  if(is.null(overlapChroms)){
    # bed_regions are non-overlapping, we can just use the non-overlapping function
    message("[info intersectPositionsAndBedRegions] bed_table regions are not overlapping: running intersectPositionsAndBedRegions_nonOverlapping")
    return(intersectPositionsAndBedRegions_nonOverlapping(positions = positions,
                                                          bed_table = bed_table))
  }else{
    # bed_regions are overlapping, we need to assign the regions to non-overlapping sets
    message("[info intersectPositionsAndBedRegions] bed_table regions are overlapping: assigning regions to non-overlapping sets and running separately")
    assignedSets <- assignBedRegionsToNonOverlappingSets(bed_table = bed_table)
    assignedSets <- reverseIdMap(assignedSets)
    # index
    rownames(bed_table) <- bed_table$id
    # initialise id maps
    idMapPositionsToRegions <- list()
    idMapRegionsToPositions <- list()
    # run the intersection for each set
    for(i in 1:length(assignedSets)){
      # i <- 1
      message("[info intersectPositionsAndBedRegions] running intersection with non-overlapping set ",i," of ",length(assignedSets))
      si <- names(assignedSets)[i]
      ids <- assignedSets[[si]]
      tmpres <- intersectPositionsAndBedRegions_nonOverlapping(positions = positions,
                                                               bed_table = bed_table[ids,,drop=F],
                                                               computeStats = FALSE)
      idMapPositionsToRegions <- mergeIdMaps(idMap1 = idMapPositionsToRegions,
                                             idMap2 = tmpres$idMapPositionsToRegions)
      idMapRegionsToPositions <- mergeIdMaps(idMap1 = idMapRegionsToPositions,
                                             idMap2 = tmpres$idMapRegionsToPositions)
    }
    
    message("[info intersectPositionsAndBedRegions] calculating intersect stats")
    # now just get the stats
    res_stats <- intersectionStatsComplete(idMap1to2 = idMapPositionsToRegions,
                                           idMap2to1 = idMapRegionsToPositions,
                                           idclassmap1 = positions,
                                           idclassmap2 = bed_table)
    
    # collect results in the return object
    returnObj <- list()
    
    returnObj$idMapPositionsToRegions <- idMapPositionsToRegions
    returnObj$idMapRegionsToPositions <- idMapRegionsToPositions

    returnObj$annotatedPositions <- res_stats$idclassmap1_updated
    returnObj$annotatedBedRegions <- res_stats$idclassmap2_updated
    
    returnObj$totalPostionsInAnyRegion <- res_stats$totalId1matchingAnyId2
    returnObj$totalRegionsAtAnyPosition <- res_stats$totalId2matchingAnyId1
    
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
}



#' Intersect two sets of non-overlapping bed regions
#'
#' Given two tables of bed regions, find the regions from the first table that
#' overlap the regions in the second table, and conversely the regions in the second
#' table that overlap the regions in the first table.
#' If the regions have classes, then find how many regions in the first table for each 
#' class of regions are contained in each region or class of regions in the second table,
#' and find how many regions in the second table for each class of regions contain each
#' region or class of regions in the second table.
#' This function is restricted to only non-overlapping bed regions. For overlapping
#' bed regions, use the more general function intersectBed.
#' 
#' @param bed_table1 data frame containing bed regions, with required columns chr, start, end, id and optionally class. Value in the id column must be unique. Regions cannot overlap.
#' @param bed_table2 data frame containing bed regions, with required columns chr, start, end, id and optionally class. Value in the id column must be unique. Regions cannot overlap.
#' @param computeStats if FALSE, intersect stats will not be calculated and only the id maps mapping
#' position ids to region ids and viceversa will be returned. This is meant to save compute time when
#' intersectBed_nonOverlapping is run inside intersectBed.
#' @return object with details intersection statistics
#' @export
intersectBed_nonOverlapping <- function(bed_table1,
                                        bed_table2,
                                        computeStats = TRUE){
  
  requiredcolumns <- c("chr","start","end","id")
  if(!all(requiredcolumns %in% colnames(bed_table1))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table))
    message("[error intersectBed_nonOverlapping] bed_table1 missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  if(!all(requiredcolumns %in% colnames(bed_table2))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table))
    message("[error intersectBed_nonOverlapping] bed_table2 missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }  
  
  # make sure id is character
  bed_table1$id <- as.character(bed_table1$id)
  bed_table2$id <- as.character(bed_table2$id)
  
  if (nrow(bed_table1)==0 & nrow(bed_table2)==0){
    message("[info intersectBed_nonOverlapping] None of the input bed_table contains regions.")
    com_table <- data.frame(chr=character(),
                            start=numeric(),
                            end=numeric(),
                            bedtable=character(),
                            segment1=numeric(),
                            segment2=numeric(),
                            stringsAsFactors = F)
  }else if (nrow(bed_table1)==0){
    message("[info intersectBed_nonOverlapping] Only bed_table2 contains regions, so they are all private.")
    com_table <- data.frame(chr=bed_table2$chr,
                            start=bed_table2$start,
                            end=bed_table2$end,
                            bedtable="2",
                            segment1=rep(NA,nrow(bed_table2)),
                            segment2=bed_table2$id,
                            stringsAsFactors = F)
  }else if (nrow(bed_table2)==0){
    message("[info intersectBed_nonOverlapping] Only bed_table1 contains regions, so they are all private.")
    com_table <- data.frame(chr=bed_table1$chr,
                            start=bed_table1$start,
                            end=bed_table1$end,
                            bedtable="1",
                            segment1=bed_table1$id,
                            segment2=rep(NA,nrow(bed_table1)),
                            stringsAsFactors = F)
  }else{
    # I need to check that there are no overlaps
    overlapChroms1 <- checkBedRegionsOverlap(bed_table = bed_table1)
    if(!is.null(overlapChroms1)){
      message("[error intersectBed_nonOverlapping] regions in bed_table1 should not overlap, ",
              "you can break them down using the function breakDownOverlappingBedRegions, ",
              "or you can run intersectBed instead.")
      return(NULL)
    }
    overlapChroms2 <- checkBedRegionsOverlap(bed_table = bed_table2)
    if(!is.null(overlapChroms2)){
      message("[error intersectBed_nonOverlapping] regions in bed_table2 should not overlap, ",
              "you can break them down using the function breakDownOverlappingBedRegions, ",
              "or you can run intersectBed instead.")
      return(NULL)
    }
    
    bed_table1$segmentid <- bed_table1$id
    bed_table2$segmentid <- bed_table2$id
    
    # then sort
    bed_table1 <- sortBed(bed_table = bed_table1)
    bed_table2 <- sortBed(bed_table = bed_table2)
    
    # initialise
    com_table <- NULL
    # need to check overlaps, one chromosome at a time
    chroms1 <- unique(bed_table1$chr)
    chroms2 <- unique(bed_table2$chr)
    chroms_all <- union(chroms1,chroms2)
    chroms_all <- sortChroms(chroms = chroms_all)
    for(chrom in chroms_all){
      # chrom <- "1"
      chrom_regions1 <- bed_table1[bed_table1$chr==chrom,,drop=F]
      chrom_regions2 <- bed_table2[bed_table2$chr==chrom,,drop=F]
      if (nrow(chrom_regions1)==0 & nrow(chrom_regions2)==0){
        # nothing to add
        message("[info intersectBed_nonOverlapping] No regions in chromosome ",chrom)
      }else if (nrow(chrom_regions1)==0){
        message("[info intersectBed_nonOverlapping] Only bed_table2 contains regions in chromosome ",chrom)
        com_table <- rbind(com_table,
                           data.frame(chr=chrom_regions2$chr,
                                      start=chrom_regions2$start,
                                      end=chrom_regions2$end,
                                      bedtable="2",
                                      segment1=rep(NA,nrow(chrom_regions2)),
                                      segment2=chrom_regions2$segmentid,
                                      stringsAsFactors = F))
      }else if (nrow(chrom_regions2)==0){
        message("[info intersectBed_nonOverlapping] Only bed_table1 contains regions in chromosome ",chrom)
        com_table <- rbind(com_table,
                           data.frame(chr=chrom_regions1$chr,
                                      start=chrom_regions1$start,
                                      end=chrom_regions1$end,
                                      bedtable="1",
                                      segment1=chrom_regions1$segmentid,
                                      segment2=rep(NA,nrow(chrom_regions1)),
                                      stringsAsFactors = F))
      }else{
        
        chrom_com_table <- NULL
        
        iter1 <- 1
        iter2 <- 1
        currentrow1 <- chrom_regions1[iter1,,drop=F]
        currentrow2 <- chrom_regions2[iter2,,drop=F]
        while (iter1<=nrow(chrom_regions1) & iter2<=nrow(chrom_regions2)) {
          overlapping <- ! (currentrow1$end < currentrow2$start | currentrow1$start > currentrow2$end)
          if(overlapping){
            startdiff <- currentrow1$start - currentrow2$start
            enddiff <- currentrow1$end - currentrow2$end
            sharedid1 <- currentrow1$segmentid
            sharedid2 <- currentrow2$segmentid
            if(startdiff!=0){
              # need to make a private segment
              if(startdiff>0){
                # seg2 starts first
                sharedstart <- currentrow1$start
                newrow <- data.frame(chr=chrom,
                                     start=currentrow2$start,
                                     end=sharedstart-1,
                                     bedtable="2",
                                     segment1=NA,
                                     segment2=currentrow2$segmentid,
                                     stringsAsFactors = F)
              }else{
                # seg1 starts first
                sharedstart <- currentrow2$start
                newrow <- data.frame(chr=chrom,
                                     start=currentrow1$start,
                                     end=sharedstart-1,
                                     bedtable="1",
                                     segment1=currentrow1$segmentid,
                                     segment2=NA,
                                     stringsAsFactors = F)
              }
              chrom_com_table <- rbind(chrom_com_table,newrow)
            }else{
              sharedstart <- currentrow1$start
            }
            
            leftoverSegment <- NULL
            if(enddiff!=0){
              # need to make a private segment
              if(enddiff>0){
                # seg2 ends first
                sharedend <- currentrow2$end
                currentrow1 <- data.frame(chr=chrom,
                                          start=sharedend+1,
                                          end=currentrow1$end,
                                          segmentid=currentrow1$segmentid,
                                          stringsAsFactors = F)
                iter2 <- iter2+1
                if(iter2<=nrow(chrom_regions2)) {
                  currentrow2 <- chrom_regions2[iter2,,drop=F]
                }else{
                  # there are no more rows in table 2, but there is a leftover
                  # segment for the current row in table 1
                  leftoverSegment <- data.frame(chr=chrom,
                                                start=currentrow1$start,
                                                end=currentrow1$end,
                                                bedtable="1",
                                                segment1=currentrow1$segmentid,
                                                segment2=NA,
                                                stringsAsFactors = F)
                  iter1 <- iter1+1
                }
              }else{
                # seg1 ends first
                sharedend <- currentrow1$end
                currentrow2 <- data.frame(chr=chrom,
                                          start=sharedend+1,
                                          end=currentrow2$end,
                                          segmentid=currentrow2$segmentid,
                                          stringsAsFactors = F)
                iter1 <- iter1+1
                if(iter1<=nrow(chrom_regions1)) {
                  currentrow1 <- chrom_regions1[iter1,,drop=F]
                }else{
                  # there are no more rows in table 1, but there is a leftover
                  # segment for the current row in table 2
                  leftoverSegment <- data.frame(chr=chrom,
                                                start=currentrow2$start,
                                                end=currentrow2$end,
                                                bedtable="2",
                                                segment1=NA,
                                                segment2=currentrow2$segmentid,
                                                stringsAsFactors = F)
                  iter2 <- iter2+1
                }
              }
            }else{
              sharedend <- currentrow1$end
              iter1 <- iter1+1
              if(iter1<=nrow(chrom_regions1)) currentrow1 <- chrom_regions1[iter1,,drop=F]
              iter2 <- iter2+1
              if(iter2<=nrow(chrom_regions2)) currentrow2 <- chrom_regions2[iter2,,drop=F]
            }
            newrow <- data.frame(chr=chrom,
                                 start=sharedstart,
                                 end=sharedend,
                                 bedtable="shared",
                                 segment1=sharedid1,
                                 segment2=sharedid2,
                                 stringsAsFactors = F)
            chrom_com_table <- rbind(chrom_com_table,newrow)
            if(!is.null(leftoverSegment)) chrom_com_table <- rbind(chrom_com_table,leftoverSegment)
          }else{
            # there is no overlapping between these two segments
            # add the one that comes first
            startdiff <- currentrow1$start - currentrow2$start
            if(startdiff>0){
              # add seg2
              newrow <- data.frame(chr=chrom,
                                   start=currentrow2$start,
                                   end=currentrow2$end,
                                   bedtable="2",
                                   segment1=NA,
                                   segment2=currentrow2$segmentid,
                                   stringsAsFactors = F)
              chrom_com_table <- rbind(chrom_com_table,newrow)
              iter2 <- iter2+1
              if(iter2<=nrow(chrom_regions2)) currentrow2 <- chrom_regions2[iter2,,drop=F]
            }else{
              # add seg1
              newrow <- data.frame(chr=chrom,
                                   start=currentrow1$start,
                                   end=currentrow1$end,
                                   bedtable="1",
                                   segment1=currentrow1$segmentid,
                                   segment2=NA,
                                   stringsAsFactors = F)
              chrom_com_table <- rbind(chrom_com_table,newrow)
              iter1 <- iter1+1
              if(iter1<=nrow(chrom_regions1)) currentrow1 <- chrom_regions1[iter1,,drop=F]
            }
          }
          
          # now it may be that either one or both chrom regions are exhausted
          # if iter1>nrow(chrom_regions1) & iter2>nrow(chrom_regions2) there is nothing to do here
          # if iter1<=nrow(chrom_regions1) & iter2<=nrow(chrom_regions2) there is nothing to do here
          # if iter1>nrow(chrom_regions1) & iter2<=nrow(chrom_regions2) then add all remaining chrom_regions2
          # if iter1<=nrow(chrom_regions1) & iter2>nrow(chrom_regions2) then add all remaining chrom_regions1
          if(iter1>nrow(chrom_regions1) & iter2<=nrow(chrom_regions2)){
            newrows <- data.frame(chr=chrom_regions2$chr[iter2:nrow(chrom_regions2)],
                                  start=chrom_regions2$start[iter2:nrow(chrom_regions2)],
                                  end=chrom_regions2$end[iter2:nrow(chrom_regions2)],
                                  bedtable="2",
                                  segment1=rep(NA,length(iter2:nrow(chrom_regions2))),
                                  segment2=chrom_regions2$segmentid[iter2:nrow(chrom_regions2)],
                                  stringsAsFactors = F)
            chrom_com_table <- rbind(chrom_com_table,newrows)
          }else if(iter1<=nrow(chrom_regions1) & iter2>nrow(chrom_regions2)){
            newrows <- data.frame(chr=chrom_regions1$chr[iter1:nrow(chrom_regions1)],
                                  start=chrom_regions1$start[iter1:nrow(chrom_regions1)],
                                  end=chrom_regions1$end[iter1:nrow(chrom_regions1)],
                                  bedtable="1",
                                  segment1=chrom_regions1$segmentid[iter1:nrow(chrom_regions1)],
                                  segment2=rep(NA,length(iter1:nrow(chrom_regions1))),
                                  stringsAsFactors = F)
            chrom_com_table <- rbind(chrom_com_table,newrows)
          }
          
        }
        com_table <- rbind(com_table,chrom_com_table)
      }
    }
  }
  
  # add single class if missing
  if(!("class" %in% colnames(bed_table1))){
    bed_table1$class <- "anyRegion"
  }
  if(!("class" %in% colnames(bed_table2))){
    bed_table2$class <- "anyRegion"
  }
  
  # annotation lists
  idMapRegions1ToRegions2 <- list()
  idMapRegions2ToRegions1 <- list()
  
  sharedOnly <- com_table[com_table$bedtable=="shared",,drop=F]
  if(nrow(sharedOnly)>0){
    for(i in 1:nrow(sharedOnly)){
      idMapRegions1ToRegions2[[sharedOnly$segment1[i]]] <- c(idMapRegions1ToRegions2[[sharedOnly$segment1[i]]],sharedOnly$segment2[i])
      idMapRegions2ToRegions1[[sharedOnly$segment2[i]]] <- c(idMapRegions2ToRegions1[[sharedOnly$segment2[i]]],sharedOnly$segment1[i])
    }
  }
  
  if(computeStats){
    if(nrow(bed_table1)>0 & nrow(bed_table2)>0){
      res_stats <- intersectionStatsComplete(idMap1to2 = idMapRegions1ToRegions2,
                                             idMap2to1 = idMapRegions2ToRegions1,
                                             idclassmap1 = bed_table1,
                                             idclassmap2 = bed_table2)
    }else{
      message("[warning intersectBed_nonOverlapping] cannot calculate intersection stats because bed_table1 and/or bed_table2 are empty")
    }

  }
  
  # return object
  returnObj <- list()
  returnObj$intersectionTable <- com_table
  
  returnObj$idMapRegions1ToRegions2 <- idMapRegions1ToRegions2
  returnObj$idMapRegions2ToRegions1 <- idMapRegions2ToRegions1
  
  if(computeStats & nrow(bed_table1)>0 & nrow(bed_table2)>0){
    returnObj$annotatedBedRegions1 <- res_stats$idclassmap1_updated
    returnObj$annotatedBedRegions2 <- res_stats$idclassmap2_updated
    
    returnObj$totalRegions1overlappingAnyRegion2 <- res_stats$totalId1matchingAnyId2
    returnObj$totalRegions2overlappingAnyRegion1 <- res_stats$totalId2matchingAnyId1
    
    returnObj$countsTable_regions1overlappingEachRegion2 <- res_stats$countsTable_classes1_in_id2
    returnObj$countsTable_regions1overlappingEachRegion2_total <- res_stats$countsTable_total1_in_id2
    returnObj$countsTable_regions1overlappingRegion2classes <- res_stats$countsTable_classes1_in_classes2
    returnObj$countsTable_regions1overlappingRegion2classes_total <- res_stats$countsTable_total1_in_classes2
    
    returnObj$countsTable_regions2overlappingEachRegion1 <- res_stats$countsTable_classes2_in_id1
    returnObj$countsTable_regions2overlappingEachRegion1_total <- res_stats$countsTable_total2_in_id1
    returnObj$countsTable_regions2overlappingRegion1classes <- res_stats$countsTable_classes2_in_classes1
    returnObj$countsTable_regions2overlappingRegion1classes_total <- res_stats$countsTable_total2_in_classes1
  }
  
  return(returnObj)
}



#' Intersect two sets of bed regions
#'
#' Given two tables of bed regions, find the regions from the first table that
#' overlap the regions in the second table, and conversely the regions in the second
#' table that overlap the regions in the first table.
#' If the regions have classes, then find how many regions in the first table for each 
#' class of regions are contained in each region or class of regions in the second table,
#' and find how many regions in the second table for each class of regions contain each
#' region or class of regions in the second table.
#' 
#' @param bed_table1 data frame containing bed regions, with required columns chr, start, end, id and optionally class. Value in the id column must be unique.
#' @param bed_table2 data frame containing bed regions, with required columns chr, start, end, id and optionally class. Value in the id column must be unique.
#' @return object with details intersection statistics
#' @export
intersectBed <- function(bed_table1,
                         bed_table2){
  
  requiredcolumns <- c("chr","start","end","id")
  if(!all(requiredcolumns %in% colnames(bed_table1))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table))
    message("[error intersectBed] bed_table1 missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  if(!all(requiredcolumns %in% colnames(bed_table2))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table))
    message("[error intersectBed] bed_table2 missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  
  # annotation lists
  idMapRegions1ToRegions2 <- list()
  idMapRegions2ToRegions1 <- list()
  
  if(nrow(bed_table1)>0 & nrow(bed_table2)>0){
    # make sure id is character
    bed_table1$id <- as.character(bed_table1$id)
    bed_table2$id <- as.character(bed_table2$id)
    
    # index
    rownames(bed_table1) <- bed_table1$id
    rownames(bed_table2) <- bed_table2$id
    
    # add single class if missing
    if(!("class" %in% colnames(bed_table1))){
      bed_table1$class <- "anyRegion"
    }
    if(!("class" %in% colnames(bed_table2))){
      bed_table2$class <- "anyRegion"
    }
    
    message("[info intersectBed] Assigning bed_table1 regions to non-overlapping sets")
    assignedSets1 <- assignBedRegionsToNonOverlappingSets(bed_table = bed_table1)
    assignedSets1 <- reverseIdMap(assignedSets1)
    
    message("[info intersectBed] Assigning bed_table2 regions to non-overlapping sets")
    assignedSets2 <- assignBedRegionsToNonOverlappingSets(bed_table = bed_table2)
    assignedSets2 <- reverseIdMap(assignedSets2)
    
    # run the intersection for each set
    for(i in 1:length(assignedSets1)){
      # i <- 1
      si <- names(assignedSets1)[i]
      ids1 <- assignedSets1[[si]]
      for(j in 1:length(assignedSets2)){
        # j <- 1
        message("[info intersectBed] running intersection between non-overlapping set ",i,
                " of ",length(assignedSets1)," from bed_table1 and non-overlapping set ",j,
                " of ",length(assignedSets2)," from bed_table2")
        sj <- names(assignedSets2)[j]
        ids2 <- assignedSets2[[sj]]
        tmpres <- intersectBed_nonOverlapping(bed_table1 = bed_table1[ids1,,drop=F],
                                              bed_table2 = bed_table2[ids2,,drop=F],
                                              computeStats = FALSE)
        idMapRegions1ToRegions2 <- mergeIdMaps(idMap1 = idMapRegions1ToRegions2,
                                               idMap2 = tmpres$idMapRegions1ToRegions2)
        idMapRegions2ToRegions1 <- mergeIdMaps(idMap1 = idMapRegions2ToRegions1,
                                               idMap2 = tmpres$idMapRegions2ToRegions1)
      }
    }
    
    res_stats <- intersectionStatsComplete(idMap1to2 = idMapRegions1ToRegions2,
                                           idMap2to1 = idMapRegions2ToRegions1,
                                           idclassmap1 = bed_table1,
                                           idclassmap2 = bed_table2)
  }else{
    message("[warning intersectBed] cannot calculate intersection stats because bed_table1 and/or bed_table2 are empty")
  }
  
  # return object
  returnObj <- list()
  
  returnObj$idMapRegions1ToRegions2 <- idMapRegions1ToRegions2
  returnObj$idMapRegions2ToRegions1 <- idMapRegions2ToRegions1
  
  if(nrow(bed_table1)>0 & nrow(bed_table2)>0){
    returnObj$annotatedBedRegions1 <- res_stats$idclassmap1_updated
    returnObj$annotatedBedRegions2 <- res_stats$idclassmap2_updated
    
    returnObj$totalRegions1overlappingAnyRegion2 <- res_stats$totalId1matchingAnyId2
    returnObj$totalRegions2overlappingAnyRegion1 <- res_stats$totalId2matchingAnyId1
    
    returnObj$countsTable_regions1overlappingEachRegion2 <- res_stats$countsTable_classes1_in_id2
    returnObj$countsTable_regions1overlappingEachRegion2_total <- res_stats$countsTable_total1_in_id2
    returnObj$countsTable_regions1overlappingRegion2classes <- res_stats$countsTable_classes1_in_classes2
    returnObj$countsTable_regions1overlappingRegion2classes_total <- res_stats$countsTable_total1_in_classes2
    
    returnObj$countsTable_regions2overlappingEachRegion1 <- res_stats$countsTable_classes2_in_id1
    returnObj$countsTable_regions2overlappingEachRegion1_total <- res_stats$countsTable_total2_in_id1
    returnObj$countsTable_regions2overlappingRegion1classes <- res_stats$countsTable_classes2_in_classes1
    returnObj$countsTable_regions2overlappingRegion1classes_total <- res_stats$countsTable_total2_in_classes1
  }

  return(returnObj)
}
