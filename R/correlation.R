
#' Correlate positions with bed regions
#'
#' Given a set of positions and a set of regions, calculate the overlap between
#' the positions and the regions (see the function intersectPositionsAndBedRegions).
#' If nsamples is greater than zero, then positions and/or regions are resampled 
#' randomly so that a NULL distribution of overlaps can be obtained. Notice that
#' the the bigger nsamples is, the more accurate the p-value for rejecting the NULL
#' hypothesis will be. A returned p-value of 0 just means that the p-value reached
#' the limit of decimal numbers that can be obtained and should be interpreted as
#' "< 1/nsamples".
#' 
#' @param positions data frame with required columns: chr, position, id
#' @param bed_regions data frame with required columns: chr, start, end, id
#' @param nsamples number of resampling used to determine the NULL distribution
#' @param altHypothesis greatherthan (the default) or lowerthan can be used 
#' @param resamplePositionsFlag if TRUE then positions will be resampled to determine the NULL distribution. Note that at least one between resamplePositionsFlag and resampleBedRegionsFlag should be TRUE
#' @param resampleBedRegionsFlag if TRUE then bed_table will be resampled to determine the NULL distribution. Note that at least one between resamplePositionsFlag and resampleBedRegionsFlag should be TRUE
#' @param resampleBedRegionsAllowOverlap allow overlapping bed regions when resampling bed_table
#' @param genomev hg19 or hg38
#' @param samplingRegions supply your own sampling regions for resampling positions and/or bed_table.
#' @param randomSeed set a random seed for the resampling 
#' @param nparallel how many parallel processes to use when running the resampling 
#' @return data frame of ordered positions
#' @export
correlatePositionsWithBedRegions <- function(positions,
                                             bed_table,
                                             nsamples=0,
                                             altHypothesis="greaterthan",
                                             resamplePositionsFlag=TRUE,
                                             resampleBedRegionsFlag=TRUE,
                                             resampleBedRegionsAllowOverlap=TRUE,
                                             genomev="hg19",
                                             samplingRegions=NULL,
                                             randomSeed=NULL,
                                             nparallel=1){
  
  # check required columns
  requiredcolumns <- c("chr","position","id")
  if(!all(requiredcolumns %in% colnames(positions))){
    missingcolumns <- setdiff(requiredcolumns,colnames(positions))
    message("[error correlatePositionsWithBedRegions] positions table missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  requiredcolumns <- c("chr","start","end","id")
  if(!all(requiredcolumns %in% colnames(bed_table))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table))
    message("[error correlatePositionsWithBedRegions] bed_table missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  
  # add single class if missing
  if(!("class" %in% colnames(positions))){
    positions$class <- "anyPosition"
  }
  if(!("class" %in% colnames(bed_table))){
    bed_table$class <- "anyRegion"
  }
  
  # check if both genomev and samplingRegions have been specified and give a warning
  if(!is.null(genomev) & !is.null(samplingRegions)){
    message("[warning correlatePositionsWithBedRegions] both genomev and samplingRegions have been specified,",
            " genomev will be ignore and the samplingRegions table will be used to sample the positions.")
    genomev <- NULL
  }
  
  # sort
  positions <- sortPositions(positions)
  bed_table <- sortBed(bed_table)
  
  # find overlaps
  res_assign <- intersectPositionsAndBedRegions(positions = positions,
                                                bed_table  = bed_table)
  
  # prepare info
  totalPostionsInAnyRegion <- res_assign$totalPostionsInAnyRegion
  totalRegionsAtAnyPosition <- res_assign$totalRegionsAtAnyPosition
  countsTable_positionsInRegionClasses <- res_assign$countsTable_positionsInRegionClasses
  countsTable_regionsAtPositionClasses <- res_assign$countsTable_regionsAtPositionClasses
  
  # update with annotations
  annotatedPositions <- res_assign$annotatedPositions
  annotatedBedRegions <- res_assign$annotatedBedRegions
  
  # summary table
  summaryTable <- data.frame(row.names = c("positions","regions"),
                             noverlap=c(totalPostionsInAnyRegion,totalRegionsAtAnyPosition),
                             ntotal=c(nrow(positions),nrow(bed_table)),
                             stringsAsFactors = F)
  
  # checks before simulations
  if(nsamples>0 & !resamplePositionsFlag & !resampleBedRegionsFlag){
    message("[warning correlatePositionsWithBedRegions] nsamples>0 but both ",
            "resamplePositionsFlag and resampleBedRegionsFlag are set to FALSE. ",
            "no resampling will be performed.")
    nsamples <- 0
  }
  
  # use simulations to find the significance
  # initialise
  sampled_PostionsInAnyRegion <- integer(nsamples)
  sampled_RegionsAtAnyPosition <- integer(nsamples)
  sampled_positionsInRegionClasses <- array(dim = c(nrow(countsTable_positionsInRegionClasses),ncol(countsTable_positionsInRegionClasses),nsamples),
                                            dimnames = list(rownames(countsTable_positionsInRegionClasses),colnames(countsTable_positionsInRegionClasses),1:nsamples))
  sampled_regionsAtPositionClasses <- array(dim = c(nrow(countsTable_regionsAtPositionClasses),ncol(countsTable_regionsAtPositionClasses),nsamples),
                                            dimnames = list(rownames(countsTable_regionsAtPositionClasses),colnames(countsTable_regionsAtPositionClasses),1:nsamples))
  if(nsamples>0){
    # set RNGkind to avoid warning
    RNGkind("L'Ecuyer-CMRG")
    doParallel::registerDoParallel(nparallel)
    if(!is.null(randomSeed)){
      # set random seed now so no need to set it later
      doRNG::registerDoRNG(randomSeed)
    }
    
    res_list <- foreach::foreach(i=1:nsamples) %dorng% {
      message("[info correlatePositionsWithBedRegions] resampling ",i," of ",nsamples)
      if(resamplePositionsFlag){
        resampled_positions <- resamplePositions(positions = positions,
                                                 genomev = genomev,
                                                 randomSeed = NULL,
                                                 samplingRegions = samplingRegions)
      }else{
        resampled_positions <- positions
      }

      # I should resample the original regions and then extend
      if(resampleBedRegionsFlag){
        resampled_bed_table <- resampleBedRegions(bed_table = bed_table,
                                                  genomev = genomev,
                                                  randomSeed = NULL,
                                                  samplingRegions = samplingRegions,
                                                  allowRegionsOverlap = resampleBedRegionsAllowOverlap)

      }else{
        resampled_bed_table <- bed_table
      }
      
      # now get the stats
      res_sample_assign <- intersectPositionsAndBedRegions(positions = resampled_positions,
                                                           bed_table = resampled_bed_table)
      # combine and return
      returnObj <- list()
      returnObj$sampled_PostionsInAnyRegion <- res_sample_assign$totalPostionsInAnyRegion
      returnObj$sampled_RegionsAtAnyPosition <- res_sample_assign$totalRegionsAtAnyPosition
      returnObj$sampled_positionsInRegionClasses <- as.matrix(res_sample_assign$countsTable_positionsInRegionClasses[rownames(countsTable_positionsInRegionClasses),colnames(countsTable_positionsInRegionClasses),drop=F])
      returnObj$sampled_regionsAtPositionClasses <- as.matrix(res_sample_assign$countsTable_regionsAtPositionClasses[rownames(countsTable_regionsAtPositionClasses),colnames(countsTable_regionsAtPositionClasses),drop=F])
      return(returnObj)
    }
    
    # reorganise
    for(i in 1:length(res_list)){
      # i <- 1
      sampled_PostionsInAnyRegion[i] <- res_list[[i]]$sampled_PostionsInAnyRegion
      sampled_RegionsAtAnyPosition[i] <- res_list[[i]]$sampled_RegionsAtAnyPosition
      sampled_positionsInRegionClasses[,,i] <- res_list[[i]]$sampled_positionsInRegionClasses
      sampled_regionsAtPositionClasses[,,i] <- res_list[[i]]$sampled_regionsAtPositionClasses
      
    }
  }
  
  # return object
  returnObj <- list()
  returnObj$annotatedPositions <- annotatedPositions
  returnObj$annotatedBedRegions <- annotatedBedRegions
  returnObj$summaryOverlaps <- summaryTable
  returnObj$nsamples <- nsamples
  returnObj$countsTable_positionsInRegionClasses <- countsTable_positionsInRegionClasses
  returnObj$countsTable_regionsAtPositionClasses <- countsTable_regionsAtPositionClasses
  
  if(nsamples>0){
    message("[info correlatePositionsWithBedRegions] calculating p-values... ")
    
    returnObj$sampled_PostionsInAnyRegion <- sampled_PostionsInAnyRegion
    returnObj$sampled_RegionsAtAnyPosition <- sampled_RegionsAtAnyPosition
    
    if(altHypothesis=="greaterthan"){
      returnObj$pvalue_PostionsInAnyRegion <- sum(totalPostionsInAnyRegion<=sampled_PostionsInAnyRegion)/nsamples
      returnObj$pvalue_RegionsAtAnyPosition <- sum(totalRegionsAtAnyPosition<=sampled_RegionsAtAnyPosition)/nsamples
    }else if(altHypothesis=="lowerthan"){
      returnObj$pvalue_PostionsInAnyRegion <- sum(totalPostionsInAnyRegion>=sampled_PostionsInAnyRegion)/nsamples
      returnObj$pvalue_RegionsAtAnyPosition <- sum(totalRegionsAtAnyPosition>=sampled_RegionsAtAnyPosition)/nsamples
    }else{
      message("[info correlatePositionsWithBedRegions] unknown alternative hypothesis ",altHypothesis,", please use greaterthan or lowerthan.")
    }
    
    # calculate all p-values the easy way
    returnObj$sampled_positionsInRegionClasses <- sampled_positionsInRegionClasses
    returnObj$sampled_regionsAtPositionClasses <- sampled_regionsAtPositionClasses
    
    pvalue_positionsInRegionClasses <- matrix(nrow = nrow(countsTable_positionsInRegionClasses),ncol = ncol(countsTable_positionsInRegionClasses),
                                              dimnames = list(rownames(countsTable_positionsInRegionClasses),colnames(countsTable_positionsInRegionClasses)))
    pvalue_regionsAtPositionClasses <- matrix(nrow = nrow(countsTable_regionsAtPositionClasses),ncol = ncol(countsTable_regionsAtPositionClasses),
                                              dimnames = list(rownames(countsTable_regionsAtPositionClasses),colnames(countsTable_regionsAtPositionClasses)))
    if(altHypothesis=="greaterthan"){
      for(i in 1:nrow(countsTable_positionsInRegionClasses)){
        # i <- 1
        for(j in 1:ncol(countsTable_positionsInRegionClasses)){
          # j <- 1
          pvalue_positionsInRegionClasses[i,j] <- sum(countsTable_positionsInRegionClasses[i,j]<=sampled_positionsInRegionClasses[i,j,])/nsamples
        }
      }
      for(i in 1:nrow(countsTable_regionsAtPositionClasses)){
        # i <- 1
        for(j in 1:ncol(countsTable_regionsAtPositionClasses)){
          # j <- 1
          pvalue_regionsAtPositionClasses[i,j] <- sum(countsTable_regionsAtPositionClasses[i,j]<=sampled_regionsAtPositionClasses[i,j,])/nsamples
        }
      }
      returnObj$pvalue_positionsInRegionClasses <- pvalue_positionsInRegionClasses
      returnObj$pvalue_regionsAtPositionClasses <- pvalue_regionsAtPositionClasses
    }else if(altHypothesis=="lowerthan"){
      for(i in 1:nrow(countsTable_positionsInRegionClasses)){
        # i <- 1
        for(j in 1:ncol(countsTable_positionsInRegionClasses)){
          # j <- 1
          pvalue_positionsInRegionClasses[i,j] <- sum(countsTable_positionsInRegionClasses[i,j]>=sampled_positionsInRegionClasses[i,j,])/nsamples
        }
      }
      for(i in 1:nrow(countsTable_regionsAtPositionClasses)){
        # i <- 1
        for(j in 1:ncol(countsTable_regionsAtPositionClasses)){
          # j <- 1
          pvalue_regionsAtPositionClasses[i,j] <- sum(countsTable_regionsAtPositionClasses[i,j]>=sampled_regionsAtPositionClasses[i,j,])/nsamples
        }
      }
      returnObj$pvalue_positionsInRegionClasses <- pvalue_positionsInRegionClasses
      returnObj$pvalue_regionsAtPositionClasses <- pvalue_regionsAtPositionClasses
    }else{
      message("[info correlatePositionsWithBedRegions] unknown alternative hypothesis ",altHypothesis,", please use greaterthan or lowerthan.")
    }
    
    # calculate also the mean/median/sd of counts
    returnObj$median_PostionsInAnyRegion <- median(sampled_PostionsInAnyRegion)
    returnObj$mean_PostionsInAnyRegion <- mean(sampled_PostionsInAnyRegion)
    returnObj$sd_PostionsInAnyRegion <- sd(sampled_PostionsInAnyRegion)
    returnObj$median_RegionsAtAnyPosition <- median(sampled_RegionsAtAnyPosition)
    returnObj$mean_RegionsAtAnyPosition <- mean(sampled_RegionsAtAnyPosition)
    returnObj$sd_RegionsAtAnyPosition <- sd(sampled_RegionsAtAnyPosition)
    
    median_positionsInRegionClasses <- matrix(nrow = nrow(countsTable_positionsInRegionClasses),ncol = ncol(countsTable_positionsInRegionClasses),
                                              dimnames = list(rownames(countsTable_positionsInRegionClasses),colnames(countsTable_positionsInRegionClasses)))
    mean_positionsInRegionClasses <- matrix(nrow = nrow(countsTable_positionsInRegionClasses),ncol = ncol(countsTable_positionsInRegionClasses),
                                            dimnames = list(rownames(countsTable_positionsInRegionClasses),colnames(countsTable_positionsInRegionClasses)))
    sd_positionsInRegionClasses <- matrix(nrow = nrow(countsTable_positionsInRegionClasses),ncol = ncol(countsTable_positionsInRegionClasses),
                                          dimnames = list(rownames(countsTable_positionsInRegionClasses),colnames(countsTable_positionsInRegionClasses)))
    median_regionsAtPositionClasses <- matrix(nrow = nrow(countsTable_regionsAtPositionClasses),ncol = ncol(countsTable_regionsAtPositionClasses),
                                              dimnames = list(rownames(countsTable_regionsAtPositionClasses),colnames(countsTable_regionsAtPositionClasses)))
    mean_regionsAtPositionClasses <- matrix(nrow = nrow(countsTable_regionsAtPositionClasses),ncol = ncol(countsTable_regionsAtPositionClasses),
                                            dimnames = list(rownames(countsTable_regionsAtPositionClasses),colnames(countsTable_regionsAtPositionClasses)))
    sd_regionsAtPositionClasses <- matrix(nrow = nrow(countsTable_regionsAtPositionClasses),ncol = ncol(countsTable_regionsAtPositionClasses),
                                          dimnames = list(rownames(countsTable_regionsAtPositionClasses),colnames(countsTable_regionsAtPositionClasses)))
    for(i in 1:nrow(countsTable_positionsInRegionClasses)){
      # i <- 1
      for(j in 1:ncol(countsTable_positionsInRegionClasses)){
        # j <- 1
        median_positionsInRegionClasses[i,j] <- median(sampled_positionsInRegionClasses[i,j,])
        mean_positionsInRegionClasses[i,j] <- mean(sampled_positionsInRegionClasses[i,j,])
        sd_positionsInRegionClasses[i,j] <- sd(sampled_positionsInRegionClasses[i,j,])
      }
    }
    for(i in 1:nrow(countsTable_regionsAtPositionClasses)){
      # i <- 1
      for(j in 1:ncol(countsTable_regionsAtPositionClasses)){
        # j <- 1
        median_regionsAtPositionClasses[i,j] <- median(sampled_regionsAtPositionClasses[i,j,])
        mean_regionsAtPositionClasses[i,j] <- mean(sampled_regionsAtPositionClasses[i,j,])
        sd_regionsAtPositionClasses[i,j] <- sd(sampled_regionsAtPositionClasses[i,j,])
      }
    }
    returnObj$median_positionsInRegionClasses <- median_positionsInRegionClasses
    returnObj$mean_positionsInRegionClasses <- mean_positionsInRegionClasses
    returnObj$sd_positionsInRegionClasses <- sd_positionsInRegionClasses
    returnObj$median_regionsAtPositionClasses <- median_regionsAtPositionClasses
    returnObj$mean_regionsAtPositionClasses <- mean_regionsAtPositionClasses
    returnObj$sd_regionsAtPositionClasses <- sd_regionsAtPositionClasses
  }
  return(returnObj)
}






#' Correlate two sets of bed regions
#'
#' Given two sets of regions, calculate the overlap between the regions in the first
#' set and the regions in the second set (see the function intersectBed).
#' If nsamples is greater than zero, then the regions are resampled 
#' randomly so that a NULL distribution of overlaps can be obtained. Notice that
#' the the bigger nsamples is, the more accurate the p-value for rejecting the NULL
#' hypothesis will be. A returned p-value of 0 just means that the p-value reached
#' the limit of decimal numbers that can be obtained and should be interpreted as
#' "< 1/nsamples".
#' 
#' @param bed_table1 data frame with required columns: chr, start, end, id
#' @param bed_table2 data frame with required columns: chr, start, end, id
#' @param nsamples number of resampling used to determine the NULL distribution
#' @param altHypothesis greatherthan (the default) or lowerthan can be used 
#' @param resampleBedRegions1Flag if TRUE then bed_table1 will be resampled to determine the NULL distribution. Note that at least one between resampleBedRegions1Flag and resampleBedRegions2Flag should be TRUE
#' @param resampleBedRegions2Flag if TRUE then bed_table2 will be resampled to determine the NULL distribution. Note that at least one between resampleBedRegions1Flag and resampleBedRegions2Flag should be TRUE
#' @param resampleBedRegions1AllowOverlap allow overlapping bed regions when resampling bed_table1
#' @param resampleBedRegions2AllowOverlap allow overlapping bed regions when resampling bed_table2
#' @param genomev hg19 or hg38
#' @param samplingRegions supply your own sampling regions for resampling positions and/or bed_table.
#' @param randomSeed set a random seed for the resampling 
#' @param nparallel how many parallel processes to use when running the resampling 
#' @return data frame of ordered positions
#' @export
correlateBedRegions <- function(bed_table1,
                                bed_table2,
                                nsamples=0,
                                altHypothesis="greaterthan",
                                resampleBedRegions1Flag=TRUE,
                                resampleBedRegions2Flag=TRUE,
                                resampleBedRegions1AllowOverlap=TRUE,
                                resampleBedRegions2AllowOverlap=TRUE,
                                genomev="hg19",
                                samplingRegions=NULL,
                                randomSeed=NULL,
                                nparallel=1){
  
  # check required columns
  requiredcolumns <- c("chr","start","end","id")
  if(!all(requiredcolumns %in% colnames(bed_table1))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table1))
    message("[error correlateBedRegions] bed_table1 missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  if(!all(requiredcolumns %in% colnames(bed_table2))){
    missingcolumns <- setdiff(requiredcolumns,colnames(bed_table2))
    message("[error correlateBedRegions] bed_table2 missing required columns: ",paste(missingcolumns,collapse = ", "))
    return(NULL)
  }
  
  # add single class if missing
  if(!("class" %in% colnames(bed_table1))){
    bed_table1$class <- "anyRegion"
  }
  if(!("class" %in% colnames(bed_table2))){
    bed_table2$class <- "anyRegion"
  }
  
  # check if both genomev and samplingRegions have been specified and give a warning
  if(!is.null(genomev) & !is.null(samplingRegions)){
    message("[warning correlatePositionsWithBedRegions] both genomev and samplingRegions have been specified,",
            " genomev will be ignore and the samplingRegions table will be used to sample the positions.")
    genomev <- NULL
  }
  
  # sort
  bed_table1 <- sortBed(bed_table1)
  bed_table2 <- sortBed(bed_table2)
  
  # find overlaps
  res_assign <- intersectBed(bed_table1 = bed_table1,
                             bed_table2 = bed_table2)
  
  # prepare info
  totalRegions1overlappingAnyRegion2 <- res_assign$totalRegions1overlappingAnyRegion2
  totalRegions2overlappingAnyRegion1 <- res_assign$totalRegions2overlappingAnyRegion1
  countsTable_regions1overlappingRegion2classes <- res_assign$countsTable_regions1overlappingRegion2classes
  countsTable_regions2overlappingRegion1classes <- res_assign$countsTable_regions2overlappingRegion1classes
  
  # update with annotations
  annotatedBedRegions1 <- res_assign$annotatedBedRegions1
  annotatedBedRegions2 <- res_assign$annotatedBedRegions2
  
  # summary table
  summaryTable <- data.frame(row.names = c("bed_table1","bed_table2"),
                             noverlap=c(totalRegions1overlappingAnyRegion2,totalRegions2overlappingAnyRegion1),
                             ntotal=c(nrow(bed_table1),nrow(bed_table2)),
                             stringsAsFactors = F)
  
  # checks before simulations
  if(nsamples>0 & !resampleBedRegions1Flag & !resampleBedRegions2Flag){
    message("[warning correlateBedRegions] nsamples>0 but both ",
            "resampleBedRegions1Flag and resampleBedRegions2Flag are set to FALSE. ",
            "no resampling will be performed.")
    nsamples <- 0
  }
  
  # use simulations to find the significance
  # initialise
  sampled_Regions1overlappingAnyRegion2 <- integer(nsamples)
  sampled_Regions2overlappingAnyRegion1 <- integer(nsamples)
  sampled_regions1overlappingRegion2classes <- array(dim = c(nrow(countsTable_regions1overlappingRegion2classes),ncol(countsTable_regions1overlappingRegion2classes),nsamples),
                                                     dimnames = list(rownames(countsTable_regions1overlappingRegion2classes),colnames(countsTable_regions1overlappingRegion2classes),1:nsamples))
  sampled_regions2overlappingRegion1classes <-  array(dim = c(nrow(countsTable_regions2overlappingRegion1classes),ncol(countsTable_regions2overlappingRegion1classes),nsamples),
                                                      dimnames = list(rownames(countsTable_regions2overlappingRegion1classes),colnames(countsTable_regions2overlappingRegion1classes),1:nsamples))
  if(nsamples>0){
    # set RNGkind to avoid warning
    RNGkind("L'Ecuyer-CMRG")
    doParallel::registerDoParallel(nparallel)
    if(!is.null(randomSeed)){
      # set random seed now so no need to set it later
      doRNG::registerDoRNG(randomSeed)
    }
    
    res_list <- foreach::foreach(i=1:nsamples) %dorng% {
      message("[info correlateBedRegions] resampling ",i," of ",nsamples)
      
      if(resampleBedRegions1Flag){
        resampled_bed_table1 <- resampleBedRegions(bed_table = bed_table1,
                                                   genomev = genomev,
                                                   randomSeed = NULL,
                                                   samplingRegions = samplingRegions,
                                                   allowRegionsOverlap = resampleBedRegions1AllowOverlap)
        
      }else{
        resampled_bed_table1 <- bed_table1
      }
      
      if(resampleBedRegions2Flag){
        resampled_bed_table2 <- resampleBedRegions(bed_table = bed_table2,
                                                  genomev = genomev,
                                                  randomSeed = NULL,
                                                  samplingRegions = samplingRegions,
                                                  allowRegionsOverlap = resampleBedRegions2AllowOverlap)
        
      }else{
        resampled_bed_table2 <- bed_table2
      }
      
      # now get the stats
      res_sample_assign <- intersectBed(bed_table1 = resampled_bed_table1,
                                        bed_table2 = resampled_bed_table2)
      # combine and return
      returnObj <- list()
      returnObj$sampled_Regions1overlappingAnyRegion2 <- res_sample_assign$totalRegions1overlappingAnyRegion2
      returnObj$sampled_Regions2overlappingAnyRegion1 <- res_sample_assign$totalRegions2overlappingAnyRegion1
      returnObj$sampled_regions1overlappingRegion2classes <- as.matrix(res_sample_assign$countsTable_regions1overlappingRegion2classes[rownames(countsTable_regions1overlappingRegion2classes),colnames(countsTable_regions1overlappingRegion2classes),drop=F])
      returnObj$sampled_regions2overlappingRegion1classes <- as.matrix(res_sample_assign$countsTable_regions2overlappingRegion1classes[rownames(countsTable_regions2overlappingRegion1classes),colnames(countsTable_regions2overlappingRegion1classes),drop=F])
      return(returnObj)
    }
    
    # reorganise
    for(i in 1:length(res_list)){
      # i <- 1
      sampled_Regions1overlappingAnyRegion2[i] <- res_list[[i]]$sampled_Regions1overlappingAnyRegion2
      sampled_Regions2overlappingAnyRegion1[i] <- res_list[[i]]$sampled_Regions2overlappingAnyRegion1
      sampled_regions1overlappingRegion2classes[,,i] <- res_list[[i]]$sampled_regions1overlappingRegion2classes
      sampled_regions2overlappingRegion1classes[,,i] <- res_list[[i]]$sampled_regions2overlappingRegion1classes
      
    }
  }
  
  # return object
  returnObj <- list()
  returnObj$annotatedBedRegions1 <- annotatedBedRegions1
  returnObj$annotatedBedRegions2 <- annotatedBedRegions2
  returnObj$summaryOverlaps <- summaryTable
  returnObj$nsamples <- nsamples
  returnObj$countsTable_regions1overlappingRegion2classes <- countsTable_regions1overlappingRegion2classes
  returnObj$countsTable_regions2overlappingRegion1classes <- countsTable_regions2overlappingRegion1classes
  
  if(nsamples>0){
    message("[info correlateBedRegions] calculating p-values... ")
    
    returnObj$sampled_Regions1overlappingAnyRegion2 <- sampled_Regions1overlappingAnyRegion2
    returnObj$sampled_Regions2overlappingAnyRegion1 <- sampled_Regions2overlappingAnyRegion1
    
    if(altHypothesis=="greaterthan"){
      returnObj$pvalue_Regions1overlappingAnyRegion2 <- sum(totalRegions1overlappingAnyRegion2<=sampled_Regions1overlappingAnyRegion2)/nsamples
      returnObj$pvalue_Regions2overlappingAnyRegion1 <- sum(totalRegions2overlappingAnyRegion1<=sampled_Regions2overlappingAnyRegion1)/nsamples
    }else if(altHypothesis=="lowerthan"){
      returnObj$pvalue_Regions1overlappingAnyRegion2 <- sum(totalRegions1overlappingAnyRegion2>=sampled_Regions1overlappingAnyRegion2)/nsamples
      returnObj$pvalue_Regions2overlappingAnyRegion1 <- sum(totalRegions2overlappingAnyRegion1>=sampled_Regions2overlappingAnyRegion1)/nsamples
    }else{
      message("[info correlateBedRegions] unknown alternative hypothesis ",altHypothesis,", please use greaterthan or lowerthan.")
    }
    
    # calculate all p-values the easy way
    returnObj$sampled_regions1overlappingRegion2classes <- sampled_regions1overlappingRegion2classes
    returnObj$sampled_regions2overlappingRegion1classes <- sampled_regions2overlappingRegion1classes
    
    pvalue_regions1overlappingRegion2classes <- matrix(nrow = nrow(countsTable_regions1overlappingRegion2classes),ncol = ncol(countsTable_regions1overlappingRegion2classes),
                                                       dimnames = list(rownames(countsTable_regions1overlappingRegion2classes),colnames(countsTable_regions1overlappingRegion2classes)))
    pvalue_regions2overlappingRegion1classes <- matrix(nrow = nrow(countsTable_regions2overlappingRegion1classes),ncol = ncol(countsTable_regions2overlappingRegion1classes),
                                                       dimnames = list(rownames(countsTable_regions2overlappingRegion1classes),colnames(countsTable_regions2overlappingRegion1classes)))
    if(altHypothesis=="greaterthan"){
      for(i in 1:nrow(countsTable_regions1overlappingRegion2classes)){
        # i <- 1
        for(j in 1:ncol(countsTable_regions1overlappingRegion2classes)){
          # j <- 1
          pvalue_regions1overlappingRegion2classes[i,j] <- sum(countsTable_regions1overlappingRegion2classes[i,j]<=sampled_regions1overlappingRegion2classes[i,j,])/nsamples
        }
      }
      for(i in 1:nrow(countsTable_regions2overlappingRegion1classes)){
        # i <- 1
        for(j in 1:ncol(countsTable_regions2overlappingRegion1classes)){
          # j <- 1
          pvalue_regions2overlappingRegion1classes[i,j] <- sum(countsTable_regions2overlappingRegion1classes[i,j]<=sampled_regions2overlappingRegion1classes[i,j,])/nsamples
        }
      }
      returnObj$pvalue_regions1overlappingRegion2classes <- pvalue_regions1overlappingRegion2classes
      returnObj$pvalue_regions2overlappingRegion1classes <- pvalue_regions2overlappingRegion1classes
    }else if(altHypothesis=="lowerthan"){
      for(i in 1:nrow(countsTable_regions1overlappingRegion2classes)){
        # i <- 1
        for(j in 1:ncol(countsTable_regions1overlappingRegion2classes)){
          # j <- 1
          pvalue_regions1overlappingRegion2classes[i,j] <- sum(countsTable_regions1overlappingRegion2classes[i,j]>=sampled_regions1overlappingRegion2classes[i,j,])/nsamples
        }
      }
      for(i in 1:nrow(countsTable_regions2overlappingRegion1classes)){
        # i <- 1
        for(j in 1:ncol(countsTable_regions2overlappingRegion1classes)){
          # j <- 1
          pvalue_regions2overlappingRegion1classes[i,j] <- sum(countsTable_regions2overlappingRegion1classes[i,j]>=sampled_regions2overlappingRegion1classes[i,j,])/nsamples
        }
      }
      returnObj$pvalue_regions1overlappingRegion2classes <- pvalue_regions1overlappingRegion2classes
      returnObj$pvalue_regions2overlappingRegion1classes <- pvalue_regions2overlappingRegion1classes
    }else{
      message("[info correlateBedRegions] unknown alternative hypothesis ",altHypothesis,", please use greaterthan or lowerthan.")
    }
    
    # calculate also the mean/median/sd of counts
    returnObj$median_Regions1overlappingAnyRegion2 <- median(sampled_Regions1overlappingAnyRegion2)
    returnObj$mean_Regions1overlappingAnyRegion2 <- mean(sampled_Regions1overlappingAnyRegion2)
    returnObj$sd_Regions1overlappingAnyRegion2 <- sd(sampled_Regions1overlappingAnyRegion2)
    returnObj$median_Regions2overlappingAnyRegion1 <- median(sampled_Regions2overlappingAnyRegion1)
    returnObj$mean_Regions2overlappingAnyRegion1 <- mean(sampled_Regions2overlappingAnyRegion1)
    returnObj$sd_Regions2overlappingAnyRegion1 <- sd(sampled_Regions2overlappingAnyRegion1)
    
    median_regions1overlappingRegion2classes <- matrix(nrow = nrow(countsTable_regions1overlappingRegion2classes),ncol = ncol(countsTable_regions1overlappingRegion2classes),
                                                       dimnames = list(rownames(countsTable_regions1overlappingRegion2classes),colnames(countsTable_regions1overlappingRegion2classes)))
    mean_regions1overlappingRegion2classes <- matrix(nrow = nrow(countsTable_regions1overlappingRegion2classes),ncol = ncol(countsTable_regions1overlappingRegion2classes),
                                                     dimnames = list(rownames(countsTable_regions1overlappingRegion2classes),colnames(countsTable_regions1overlappingRegion2classes)))
    sd_regions1overlappingRegion2classes <- matrix(nrow = nrow(countsTable_regions1overlappingRegion2classes),ncol = ncol(countsTable_regions1overlappingRegion2classes),
                                                   dimnames = list(rownames(countsTable_regions1overlappingRegion2classes),colnames(countsTable_regions1overlappingRegion2classes)))
    median_regions2overlappingRegion1classes <- matrix(nrow = nrow(countsTable_regions2overlappingRegion1classes),ncol = ncol(countsTable_regions2overlappingRegion1classes),
                                                       dimnames = list(rownames(countsTable_regions2overlappingRegion1classes),colnames(countsTable_regions2overlappingRegion1classes)))
    mean_regions2overlappingRegion1classes <- matrix(nrow = nrow(countsTable_regions2overlappingRegion1classes),ncol = ncol(countsTable_regions2overlappingRegion1classes),
                                                     dimnames = list(rownames(countsTable_regions2overlappingRegion1classes),colnames(countsTable_regions2overlappingRegion1classes)))
    sd_regions2overlappingRegion1classes <- matrix(nrow = nrow(countsTable_regions2overlappingRegion1classes),ncol = ncol(countsTable_regions2overlappingRegion1classes),
                                                   dimnames = list(rownames(countsTable_regions2overlappingRegion1classes),colnames(countsTable_regions2overlappingRegion1classes)))
    for(i in 1:nrow(countsTable_regions1overlappingRegion2classes)){
      # i <- 1
      for(j in 1:ncol(countsTable_regions1overlappingRegion2classes)){
        # j <- 1
        median_regions1overlappingRegion2classes[i,j] <- median(sampled_regions1overlappingRegion2classes[i,j,])
        mean_regions1overlappingRegion2classes[i,j] <- mean(sampled_regions1overlappingRegion2classes[i,j,])
        sd_regions1overlappingRegion2classes[i,j] <- sd(sampled_regions1overlappingRegion2classes[i,j,])
      }
    }
    for(i in 1:nrow(countsTable_regions2overlappingRegion1classes)){
      # i <- 1
      for(j in 1:ncol(countsTable_regions2overlappingRegion1classes)){
        # j <- 1
        median_regions2overlappingRegion1classes[i,j] <- median(sampled_regions2overlappingRegion1classes[i,j,])
        mean_regions2overlappingRegion1classes[i,j] <- mean(sampled_regions2overlappingRegion1classes[i,j,])
        sd_regions2overlappingRegion1classes[i,j] <- sd(sampled_regions2overlappingRegion1classes[i,j,])
      }
    }
    returnObj$median_regions1overlappingRegion2classes <- median_regions1overlappingRegion2classes
    returnObj$mean_regions1overlappingRegion2classes <- mean_regions1overlappingRegion2classes
    returnObj$sd_regions1overlappingRegion2classes <- sd_regions1overlappingRegion2classes
    returnObj$median_regions2overlappingRegion1classes <- median_regions2overlappingRegion1classes
    returnObj$mean_regions2overlappingRegion1classes <- mean_regions2overlappingRegion1classes
    returnObj$sd_regions2overlappingRegion1classes <- sd_regions2overlappingRegion1classes
  }
  return(returnObj)
}





