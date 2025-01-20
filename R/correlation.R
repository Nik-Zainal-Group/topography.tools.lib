
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
#' @param resamplePositionsFlag if TRUE then positions will be resampled to determine the NULL distribution.
#' Note that at least one between resamplePositionsFlag and resampleBedRegionsFlag should be TRUE
#' @param resampleBedRegionsFlag if TRUE then bed_table will be resampled to determine the NULL distribution.
#' Note that at least one between resamplePositionsFlag and resampleBedRegionsFlag should be TRUE
#' @param resampleBedRegionsAllowOverlap allow overlapping bed regions when resampling bed_table
#' @param genomev hg19 or hg38
#' @param samplingRegions supply your own sampling regions for resampling positions and/or bed_table.
#' @param randomSeed set a random seed for the resampling 
#' @param returnResampledPositions if TRUE and if resampling of positions was performed, return a list with all the positions data frames resampled
#' @param returnResampledBedRegions if TRUE and if resampling of bed regions was performed, return a list with all the bed_tables data frames resampled
#' @param resampled_positions_list supply your own resampled position tables. Use it only if you know what you are doing. Typically useful for multiple
#' correlations testing so that the NULL distribution can be calculated only once
#' @param resampled_bed_regions_list supply your own resampled bed tables. Use it only if you know what you are doing. Typically useful for multiple
#' correlations testing so that the NULL distribution can be calculated only once 
#' @param nparallel how many parallel processes to use when running the resampling 
#' @param verbose set to FALSE to suppress the info messages. Warning and error messages will still be shown
#' @return object with correlation statistics and additional data
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
                                             returnResampledPositions=FALSE,
                                             returnResampledBedRegions=FALSE,
                                             resampled_positions_list=NULL,
                                             resampled_bed_regions_list=NULL,
                                             nparallel=1,
                                             verbose = TRUE){
  
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
  
  # checks on precomputed resampling
  precomputed_resampled_positions_list <- NULL
  if(!is.null(resampled_positions_list) & resamplePositionsFlag){
    precomputed_resampled_positions_list <- resampled_positions_list
    if(nsamples != length(precomputed_resampled_positions_list)){
      message("[error correlatePositionsWithBedRegions] attempting to use precomputed resampled positions,",
              " however nsamples and the length of resampled_positions_list differ.")
      return(NULL)
    }
  }
  
  precomputed_resampled_bed_regions_list <- NULL
  if(!is.null(resampled_bed_regions_list) & resampleBedRegionsFlag){
    precomputed_resampled_bed_regions_list <- resampled_bed_regions_list
    if(nsamples != length(precomputed_resampled_bed_regions_list)){
      message("[error correlatePositionsWithBedRegions] attempting to use precomputed resampled bed regions,",
              " however nsamples and the length of resampled_bed_regions_list differ.")
      return(NULL)
    }
  }
  
  # sort
  positions <- sortPositions(positions)
  bed_table <- sortBed(bed_table)
  
  # find overlaps
  res_assign <- intersectPositionsAndBedRegions(positions = positions,
                                                bed_table  = bed_table,
                                                verbose = verbose)
  
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
  resampled_positions_list <- list()
  resampled_bed_regions_list <- list()
  
  if(nsamples>0){
    if(verbose) message("[info correlatePositionsWithBedRegions] resampling... ")
    
    # set RNGkind to avoid warning
    RNGkind("L'Ecuyer-CMRG")
    doParallel::registerDoParallel(nparallel)
    if(!is.null(randomSeed)){
      # set random seed now so no need to set it later
      doRNG::registerDoRNG(randomSeed)
    }
    
    res_list <- foreach::foreach(i=1:nsamples) %dorng% {
      if(verbose) message("[info correlatePositionsWithBedRegions] resampling ",i," of ",nsamples)
      if(resamplePositionsFlag){
        if(is.null(precomputed_resampled_positions_list)){
          resampled_positions <- resamplePositions(positions = positions,
                                                   genomev = genomev,
                                                   randomSeed = NULL,
                                                   samplingRegions = samplingRegions)
        }else{
          resampled_positions <- precomputed_resampled_positions_list[[i]]
        }
      }else{
        resampled_positions <- positions
      }

      # I should resample the original regions and then extend
      if(resampleBedRegionsFlag){
        if(is.null(precomputed_resampled_bed_regions_list)){
          resampled_bed_table <- resampleBedRegions(bed_table = bed_table,
                                                    genomev = genomev,
                                                    randomSeed = NULL,
                                                    samplingRegions = samplingRegions,
                                                    allowRegionsOverlap = resampleBedRegionsAllowOverlap)
        }else{
          resampled_bed_table <- precomputed_resampled_bed_regions_list[[i]]
        }
      }else{
        resampled_bed_table <- bed_table
      }
      
      # now get the stats
      res_sample_assign <- intersectPositionsAndBedRegions(positions = resampled_positions,
                                                           bed_table = resampled_bed_table,
                                                           verbose = verbose)
      # combine and return
      returnObj <- list()
      returnObj$sampled_PostionsInAnyRegion <- res_sample_assign$totalPostionsInAnyRegion
      returnObj$sampled_RegionsAtAnyPosition <- res_sample_assign$totalRegionsAtAnyPosition
      returnObj$sampled_positionsInRegionClasses <- as.matrix(res_sample_assign$countsTable_positionsInRegionClasses[rownames(countsTable_positionsInRegionClasses),colnames(countsTable_positionsInRegionClasses),drop=F])
      returnObj$sampled_regionsAtPositionClasses <- as.matrix(res_sample_assign$countsTable_regionsAtPositionClasses[rownames(countsTable_regionsAtPositionClasses),colnames(countsTable_regionsAtPositionClasses),drop=F])
      if(resamplePositionsFlag & returnResampledPositions) returnObj$resampled_positions <- resampled_positions
      if(resampleBedRegionsFlag & returnResampledBedRegions) returnObj$resampled_bed_table <- resampled_bed_table
      return(returnObj)
    }
    
    # reorganise
    for(i in 1:length(res_list)){
      # i <- 1
      sampled_PostionsInAnyRegion[i] <- res_list[[i]]$sampled_PostionsInAnyRegion
      sampled_RegionsAtAnyPosition[i] <- res_list[[i]]$sampled_RegionsAtAnyPosition
      sampled_positionsInRegionClasses[,,i] <- res_list[[i]]$sampled_positionsInRegionClasses
      sampled_regionsAtPositionClasses[,,i] <- res_list[[i]]$sampled_regionsAtPositionClasses
      if(!is.null(res_list[[i]]$resampled_positions)) resampled_positions_list[[i]] <- res_list[[i]]$resampled_positions
      if(!is.null(res_list[[i]]$resampled_bed_table)) resampled_bed_regions_list[[i]] <- res_list[[i]]$resampled_bed_table
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
    if(verbose) message("[info correlatePositionsWithBedRegions] calculating p-values... ")
    
    returnObj$sampled_PostionsInAnyRegion <- sampled_PostionsInAnyRegion
    returnObj$sampled_RegionsAtAnyPosition <- sampled_RegionsAtAnyPosition
    
    if(altHypothesis=="greaterthan"){
      returnObj$pvalue_PostionsInAnyRegion <- sum(totalPostionsInAnyRegion<=sampled_PostionsInAnyRegion)/nsamples
      returnObj$pvalue_RegionsAtAnyPosition <- sum(totalRegionsAtAnyPosition<=sampled_RegionsAtAnyPosition)/nsamples
    }else if(altHypothesis=="lowerthan"){
      returnObj$pvalue_PostionsInAnyRegion <- sum(totalPostionsInAnyRegion>=sampled_PostionsInAnyRegion)/nsamples
      returnObj$pvalue_RegionsAtAnyPosition <- sum(totalRegionsAtAnyPosition>=sampled_RegionsAtAnyPosition)/nsamples
    }else{
      message("[warning correlatePositionsWithBedRegions] unknown alternative hypothesis ",altHypothesis,", please use greaterthan or lowerthan.")
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
      message("[warning correlatePositionsWithBedRegions] unknown alternative hypothesis ",altHypothesis,", please use greaterthan or lowerthan.")
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
    
    # check if I need to return the resampled positions and/or bed regions
    if(length(resampled_positions_list)>0) returnObj$resampled_positions_list <- resampled_positions_list
    if(length(resampled_bed_regions_list)>0) returnObj$resampled_bed_regions_list <- resampled_bed_regions_list
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
#' @param returnResampledBedRegions1 if TRUE and if resampling of bed_table1 was performed, return a list with all the bed_table1 data frames resampled
#' @param returnResampledBedRegions2 if TRUE and if resampling of bed_table2 was performed, return a list with all the bed_table2 data frames resampled
#' @param resampled_bed_regions1_list supply your own resampled bed_table1. Use it only if you know what you are doing. Typically useful for multiple
#' correlations testing so that the NULL distribution can be calculated only once
#' @param resampled_bed_regions2_list supply your own resampled bed_table2. Use it only if you know what you are doing. Typically useful for multiple
#' correlations testing so that the NULL distribution can be calculated only once 
#' @param nparallel how many parallel processes to use when running the resampling 
#' @param verbose set to FALSE to suppress the info messages. Warning and error messages will still be shown
#' @return object with correlation statistics and additional data
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
                                returnResampledBedRegions1=FALSE,
                                returnResampledBedRegions2=FALSE,
                                resampled_bed_regions1_list=NULL,
                                resampled_bed_regions2_list=NULL,
                                nparallel=1,
                                verbose = TRUE){
  
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
    message("[warning correlateBedRegions] both genomev and samplingRegions have been specified,",
            " genomev will be ignore and the samplingRegions table will be used to sample the positions.")
    genomev <- NULL
  }
  
  # checks on precomputed resampling
  precomputed_resampled_bed_regions1_list <- NULL
  if(!is.null(resampled_bed_regions1_list) & resampleBedRegions1Flag){
    precomputed_resampled_bed_regions1_list <- resampled_bed_regions1_list
    if(nsamples != length(precomputed_resampled_bed_regions1_list)){
      message("[error correlateBedRegions] attempting to use precomputed resampled bed regions of bed_table1,",
              " however nsamples and the length of resampled_bed_regions1_list differ.")
      return(NULL)
    }
  }
  
  precomputed_resampled_bed_regions2_list <- NULL
  if(!is.null(resampled_bed_regions2_list) & resampleBedRegions2Flag){
    precomputed_resampled_bed_regions2_list <- resampled_bed_regions2_list
    if(nsamples != length(precomputed_resampled_bed_regions2_list)){
      message("[error correlateBedRegions] attempting to use precomputed resampled bed regions of bed_table2,",
              " however nsamples and the length of resampled_bed_regions2_list differ.")
      return(NULL)
    }
  }
  
  # sort
  bed_table1 <- sortBed(bed_table1)
  bed_table2 <- sortBed(bed_table2)
  
  # find overlaps
  res_assign <- intersectBed(bed_table1 = bed_table1,
                             bed_table2 = bed_table2,
                             verbose = verbose)
  
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
  resampled_bed_regions1_list <- list()
  resampled_bed_regions2_list <- list()
  if(nsamples>0){
    if(verbose) message("[info correlateBedRegions] resampling... ")
    
    # set RNGkind to avoid warning
    RNGkind("L'Ecuyer-CMRG")
    doParallel::registerDoParallel(nparallel)
    if(!is.null(randomSeed)){
      # set random seed now so no need to set it later
      doRNG::registerDoRNG(randomSeed)
    }
    
    res_list <- foreach::foreach(i=1:nsamples) %dorng% {
      if(verbose) message("[info correlateBedRegions] resampling ",i," of ",nsamples)
      
      if(resampleBedRegions1Flag){
        if(is.null(precomputed_resampled_bed_regions1_list)){
          resampled_bed_table1 <- resampleBedRegions(bed_table = bed_table1,
                                                     genomev = genomev,
                                                     randomSeed = NULL,
                                                     samplingRegions = samplingRegions,
                                                     allowRegionsOverlap = resampleBedRegions1AllowOverlap)
        }else{
          resampled_bed_table1 <- precomputed_resampled_bed_regions1_list[[i]]
        }
      }else{
        resampled_bed_table1 <- bed_table1
      }
      
      if(resampleBedRegions2Flag){
        if(is.null(precomputed_resampled_bed_regions1_list)){
          resampled_bed_table2 <- resampleBedRegions(bed_table = bed_table2,
                                                    genomev = genomev,
                                                    randomSeed = NULL,
                                                    samplingRegions = samplingRegions,
                                                    allowRegionsOverlap = resampleBedRegions2AllowOverlap)
        }else{
          resampled_bed_table2 <- precomputed_resampled_bed_regions2_list[[i]]
        }
      }else{
        resampled_bed_table2 <- bed_table2
      }
      
      # now get the stats
      res_sample_assign <- intersectBed(bed_table1 = resampled_bed_table1,
                                        bed_table2 = resampled_bed_table2,
                                        verbose = verbose)
      # combine and return
      returnObj <- list()
      returnObj$sampled_Regions1overlappingAnyRegion2 <- res_sample_assign$totalRegions1overlappingAnyRegion2
      returnObj$sampled_Regions2overlappingAnyRegion1 <- res_sample_assign$totalRegions2overlappingAnyRegion1
      returnObj$sampled_regions1overlappingRegion2classes <- as.matrix(res_sample_assign$countsTable_regions1overlappingRegion2classes[rownames(countsTable_regions1overlappingRegion2classes),colnames(countsTable_regions1overlappingRegion2classes),drop=F])
      returnObj$sampled_regions2overlappingRegion1classes <- as.matrix(res_sample_assign$countsTable_regions2overlappingRegion1classes[rownames(countsTable_regions2overlappingRegion1classes),colnames(countsTable_regions2overlappingRegion1classes),drop=F])
      if(resampleBedRegions1Flag & returnResampledBedRegions1) returnObj$resampled_bed_table1 <- resampled_bed_table1
      if(resampleBedRegions2Flag & returnResampledBedRegions2) returnObj$resampled_bed_table2 <- resampled_bed_table2
      return(returnObj)
    }
    
    # reorganise
    for(i in 1:length(res_list)){
      # i <- 1
      sampled_Regions1overlappingAnyRegion2[i] <- res_list[[i]]$sampled_Regions1overlappingAnyRegion2
      sampled_Regions2overlappingAnyRegion1[i] <- res_list[[i]]$sampled_Regions2overlappingAnyRegion1
      sampled_regions1overlappingRegion2classes[,,i] <- res_list[[i]]$sampled_regions1overlappingRegion2classes
      sampled_regions2overlappingRegion1classes[,,i] <- res_list[[i]]$sampled_regions2overlappingRegion1classes
      if(!is.null(res_list[[i]]$resampled_bed_table1)) resampled_bed_regions1_list[[i]] <- res_list[[i]]$resampled_bed_table1
      if(!is.null(res_list[[i]]$resampled_bed_table2)) resampled_bed_regions2_list[[i]] <- res_list[[i]]$resampled_bed_table2
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
    if(verbose) message("[info correlateBedRegions] calculating p-values... ")
    
    returnObj$sampled_Regions1overlappingAnyRegion2 <- sampled_Regions1overlappingAnyRegion2
    returnObj$sampled_Regions2overlappingAnyRegion1 <- sampled_Regions2overlappingAnyRegion1
    
    if(altHypothesis=="greaterthan"){
      returnObj$pvalue_Regions1overlappingAnyRegion2 <- sum(totalRegions1overlappingAnyRegion2<=sampled_Regions1overlappingAnyRegion2)/nsamples
      returnObj$pvalue_Regions2overlappingAnyRegion1 <- sum(totalRegions2overlappingAnyRegion1<=sampled_Regions2overlappingAnyRegion1)/nsamples
    }else if(altHypothesis=="lowerthan"){
      returnObj$pvalue_Regions1overlappingAnyRegion2 <- sum(totalRegions1overlappingAnyRegion2>=sampled_Regions1overlappingAnyRegion2)/nsamples
      returnObj$pvalue_Regions2overlappingAnyRegion1 <- sum(totalRegions2overlappingAnyRegion1>=sampled_Regions2overlappingAnyRegion1)/nsamples
    }else{
      message("[warning correlateBedRegions] unknown alternative hypothesis ",altHypothesis,", please use greaterthan or lowerthan.")
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
      message("[warning correlateBedRegions] unknown alternative hypothesis ",altHypothesis,", please use greaterthan or lowerthan.")
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
    
    # check if I need to return the resampled positions and/or bed regions
    if(length(resampled_bed_regions1_list)>0) returnObj$resampled_bed_regions1_list <- resampled_bed_regions1_list
    if(length(resampled_bed_regions2_list)>0) returnObj$resampled_bed_regions2_list <- resampled_bed_regions2_list
  }
  return(returnObj)
}


#' Multiple correlations between positions and/or regions
#'
#' Given a table with a set of reference entities, which are either positions or
#' bed regions, calculate the correlation with a list of other entities, which
#' again can be positions and/or regions. The type of the entities, either positions
#' or bed regions, will be determined automatically by checking the column names,
#' where positions should have column names chr and position, while bed regions
#' should have column names chr, start and end. Statistical significance of the
#' intersection between the reference entities and other entities sets will be
#' determined by resampling the reference entities, and optionally also the entities
#' to compare to.
#' 
#' 
#' @param referenceEntities data frame of positions or bed regions. If positions are used,
#' required columns are id, chr and position. If bed regions are used, required columns are
#' id, chr, start and end.
#' @param compareEntitiesList list of data frames, with each data frame either a table of positions
#' or bed regions. The list names will be used as the name of each entities set
#' @param referenceEntitiesName name to be use when referring to the reference entities
#' @param nsamples number of resampling used to determine the NULL distribution
#' @param altHypothesis greatherthan (the default) or lowerthan can be used 
#' @param resampleCompareEntities if TRUE the compare entities will be resampled
#' @param resampleReferenceEntitiesAllowOverlap if TRUE and if the reference entities are
#' bed regions, the resampling will allow bed regions overlap
#' @param resampleCompareEntitiesAllowOverlap  if TRUE and if the compare entities are
#' bed regions, the resampling will allow bed regions overlap
#' @param genomev hg19 or hg38
#' @param samplingRegions supply your own sampling regions for resampling positions and/or bed_table.
#' @param randomSeed set a random seed for the resampling 
#' @param nparallel how many parallel processes to use when running the resampling
#' @param verbose set to FALSE to suppress the info messages. Warning and error messages will still be shown
#' @return correlation statistics and annotations
#' @export
multipleCorrelations <- function(referenceEntities,
                                 compareEntitiesList,
                                 referenceEntitiesName = "referenceEntities",
                                 nsamples = 0,
                                 altHypothesis="greaterthan",
                                 resampleCompareEntities=FALSE,
                                 resampleReferenceEntitiesAllowOverlap = TRUE,
                                 resampleCompareEntitiesAllowOverlap = TRUE,
                                 genomev = "hg19",
                                 samplingRegions = NULL,
                                 randomSeed = NULL,
                                 nparallel = 1,
                                 verbose = TRUE){
  
  # check type of referenceEntities
  etype <- getEntitiesType(entities = referenceEntities)
  if(etype=="ambiguous" | etype=="unknown"){
    message("[error multipleCorrelations] referenceEntities type is ",etype,". ",
            "Please make sure it is either a positions table, with columns chr and position, ",
            "or a bed table, with columns chr, start, end.")
    return("NULL")
  }
  
  # if we are resampling, let's do it once only for the reference
  resampled_referenceEntities <- NULL
  if(nsamples>0) {
    if(verbose) message("[info multipleCorrelations] resampling referenceEntities...")
    # set RNGkind to avoid warning
    RNGkind("L'Ecuyer-CMRG")
    doParallel::registerDoParallel(nparallel)
    if(!is.null(randomSeed)){
      # set random seed now so no need to set it later
      doRNG::registerDoRNG(randomSeed)
    }
    
    if(etype=="positions"){
      resampled_referenceEntities <- foreach::foreach(i=1:nsamples) %dorng% {
        return(resamplePositions(positions = referenceEntities,
                                 genomev = genomev,
                                 samplingRegions = samplingRegions,
                                 randomSeed = NULL,
                                 verbose = FALSE))
      }
    }else if(etype=="bedRegions"){
      resampled_referenceEntities <- foreach::foreach(i=1:nsamples) %dorng% {
        return(resampleBedRegions(bed_table = referenceEntities,
                                  genomev = genomev,
                                  samplingRegions = samplingRegions,
                                  randomSeed = NULL,
                                  allowRegionsOverlap = resampleReferenceEntitiesAllowOverlap,
                                  verbose = FALSE))
      }
    }else{
      message("[error multipleCorrelations] cannot resample referenceEntities, type ",etype,
              " not implemented.")
      return(NULL)
    }

  }
  
  # set things up
  annotatedEntities <- referenceEntities
  annotatedCompareEntitiesList <- list()
  counts_refEntitiesWithCompEntities <- as.data.frame(matrix(c(rep(NA,length(compareEntitiesList)),nrow(referenceEntities)),
                                                             nrow = 1,ncol = length(compareEntitiesList) + 1,
                                                             dimnames = list(referenceEntitiesName,c(names(compareEntitiesList),"total"))),
                                                      stringsAsFactors = F)
  expected_refEntitiesWithCompEntities <- as.data.frame(matrix(c(rep(NA,length(compareEntitiesList)),nrow(referenceEntities)),
                                                               nrow = 1,ncol = length(compareEntitiesList) + 1,
                                                               dimnames = list(referenceEntitiesName,c(names(compareEntitiesList),"total"))),
                                                        stringsAsFactors = F)
  pvalues_refEntitiesWithCompEntities <- as.data.frame(matrix(NA,nrow = 1,ncol = length(compareEntitiesList),
                                                              dimnames = list(referenceEntitiesName,names(compareEntitiesList))),
                                                       stringsAsFactors = F)
  counts_compEntitiesWithRefEntities <- as.data.frame(matrix(c(rep(NA,length(compareEntitiesList)*2)),
                                                             nrow = length(compareEntitiesList),ncol = 2,
                                                             dimnames = list(names(compareEntitiesList),c(referenceEntitiesName,"total"))),
                                                      stringsAsFactors = F)
  expected_compEntitiesWithRefEntities <- as.data.frame(matrix(c(rep(NA,length(compareEntitiesList)*2)),
                                                               nrow = length(compareEntitiesList),ncol = 2,
                                                               dimnames = list(names(compareEntitiesList),c(referenceEntitiesName,"total"))),
                                                        stringsAsFactors = F)
  pvalues_compEntitiesWithRefEntities <- as.data.frame(matrix(c(rep(NA,length(compareEntitiesList))),
                                                              nrow = length(compareEntitiesList),ncol = 1,
                                                              dimnames = list(names(compareEntitiesList),c(referenceEntitiesName))),
                                                       stringsAsFactors = F)
  cei <- 0
  for(CE in names(compareEntitiesList)){
    # CE <- names(compareEntitiesList)[1]
    cei <- cei+1
    message("[info multipleCorrelations] testing correlation with ",CE,", ",cei," of ",length(compareEntitiesList))
    currente <- compareEntitiesList[[CE]]
    cetype <- getEntitiesType(entities = currente)
    if(cetype=="ambiguous" | cetype=="unknown"){
      message("[warning multipleCorrelations] skipping compare entities ",CE,". ",
              "Entities type is ",cetype,". ",
              "Please make sure it is either a positions table, with columns chr and position, ",
              "or a bed table, with columns chr, start, end.")
    }else if(etype=="positions" & cetype=="positions"){
      message("[warning multipleCorrelations] skipping compare entities ",CE,". ",
              "Both referenceEntities and compare entities ",CE," entities type is positions. ")
    }else{
      if(etype=="positions"){
        # then the cetype must be bedRegions
        
        res_corr_pos_extend <- correlatePositionsWithBedRegions(positions = annotatedEntities,
                                                                bed_table = currente,
                                                                resamplePositionsFlag = TRUE,
                                                                resampleBedRegionsFlag = resampleCompareEntities,
                                                                resampleBedRegionsAllowOverlap = resampleCompareEntitiesAllowOverlap,
                                                                resampled_positions_list = resampled_referenceEntities,
                                                                genomev = genomev,
                                                                nsamples = nsamples,
                                                                altHypothesis = altHypothesis,
                                                                samplingRegions = samplingRegions,
                                                                nparallel = nparallel,
                                                                randomSeed = randomSeed,
                                                                verbose = verbose)
        # update annotated reference entities
        annotatedEntities <- res_corr_pos_extend$annotatedPositions
        annotatedEntities$classAnnotation <- NULL
        colnames(annotatedEntities)[which(colnames(annotatedEntities)=="nMatches")] <- paste0("n",CE)
        colnames(annotatedEntities)[which(colnames(annotatedEntities)=="idAnnotation")] <- CE
        # save some info
        # - annotated positions
        annotatedBedRegions <- res_corr_pos_extend$annotatedBedRegions
        annotatedBedRegions$classAnnotation <- NULL
        colnames(annotatedBedRegions)[which(colnames(annotatedBedRegions)=="nMatches")] <- paste0("n",referenceEntitiesName)
        colnames(annotatedBedRegions)[which(colnames(annotatedBedRegions)=="idAnnotation")] <- referenceEntitiesName
        # - summary overlaps
        summaryOverlaps <- res_corr_pos_extend$summaryOverlaps
        # some info about significance
        if(nsamples>0){
          # - p-values
          pvalue1 <- ifelse(res_corr_pos_extend$pvalue_RegionsAtAnyPosition==0,paste0("<",1/nsamples),res_corr_pos_extend$pvalue_RegionsAtAnyPosition)
          pvalue2 <- ifelse(res_corr_pos_extend$pvalue_PostionsInAnyRegion==0,paste0("<",1/nsamples),res_corr_pos_extend$pvalue_PostionsInAnyRegion)
          significanceTable <- data.frame(nsamples=res_corr_pos_extend$nsamples,
                                          pvalueRefEntities=pvalue2,
                                          pvalueCompEntities=pvalue1,
                                          stringsAsFactors = F)
          # - expected overlaps
          expectedOverlaps <- data.frame(nsamples=res_corr_pos_extend$nsamples,
                                         expectedRefEntities=res_corr_pos_extend$mean_PostionsInAnyRegion,
                                         expectedCompEntities=res_corr_pos_extend$mean_RegionsAtAnyPosition,
                                         stringsAsFactors = F)
        }
        
        # collect
        annotatedCompareEntitiesList[[CE]] <- annotatedBedRegions
        counts_refEntitiesWithCompEntities[1,CE] <- summaryOverlaps[1,1]
        counts_compEntitiesWithRefEntities[CE,1] <- summaryOverlaps[2,1]
        if(nsamples>0){
          expected_refEntitiesWithCompEntities[1,CE] <- expectedOverlaps[1,2]
          pvalues_refEntitiesWithCompEntities[1,CE] <- significanceTable[1,2]
          expected_compEntitiesWithRefEntities[CE,1] <- expectedOverlaps[1,3]
          pvalues_compEntitiesWithRefEntities[CE,1] <- significanceTable[1,3]
          counts_compEntitiesWithRefEntities[CE,2] <- nrow(currente)
          expected_compEntitiesWithRefEntities[CE,2] <- nrow(currente)
        }
      }else{
        # then etype is bedRegions and cetype can be positions or bedRegions
        if(cetype=="positions"){

          res_corr_pos_extend <- correlatePositionsWithBedRegions(positions = currente,
                                                                  bed_table = annotatedEntities,
                                                                  resamplePositionsFlag = resampleCompareEntities,
                                                                  resampleBedRegionsFlag = TRUE,
                                                                  resampleBedRegionsAllowOverlap = resampleReferenceEntitiesAllowOverlap,
                                                                  resampled_bed_regions_list = resampled_referenceEntities,
                                                                  genomev = genomev,
                                                                  nsamples = nsamples,
                                                                  altHypothesis = altHypothesis,
                                                                  samplingRegions = samplingRegions,
                                                                  nparallel = nparallel,
                                                                  randomSeed = randomSeed,
                                                                  verbose = verbose)
          # update annotated reference entities
          annotatedEntities <- res_corr_pos_extend$annotatedBedRegions
          annotatedEntities$classAnnotation <- NULL
          colnames(annotatedEntities)[which(colnames(annotatedEntities)=="nMatches")] <- paste0("n",CE)
          colnames(annotatedEntities)[which(colnames(annotatedEntities)=="idAnnotation")] <- CE
          # save some info
          # - annotated positions
          annotatedPositions <- res_corr_pos_extend$annotatedPositions
          annotatedPositions$classAnnotation <- NULL
          colnames(annotatedPositions)[which(colnames(annotatedPositions)=="nMatches")] <- paste0("n",referenceEntitiesName)
          colnames(annotatedPositions)[which(colnames(annotatedPositions)=="idAnnotation")] <- referenceEntitiesName
          # - summary overlaps
          summaryOverlaps <- res_corr_pos_extend$summaryOverlaps
          # some info about significance
          if(nsamples>0){
            # - p-values
            pvalue1 <- ifelse(res_corr_pos_extend$pvalue_RegionsAtAnyPosition==0,paste0("<",1/nsamples),res_corr_pos_extend$pvalue_RegionsAtAnyPosition)
            pvalue2 <- ifelse(res_corr_pos_extend$pvalue_PostionsInAnyRegion==0,paste0("<",1/nsamples),res_corr_pos_extend$pvalue_PostionsInAnyRegion)
            significanceTable <- data.frame(nsamples=res_corr_pos_extend$nsamples,
                                            pvalueRefEntitiesWithPositions=pvalue1,
                                            pvaluePositionsInAnyHotspot=pvalue2,
                                            stringsAsFactors = F)
            # - expected overlaps
            expectedOverlaps <- data.frame(nsamples=res_corr_pos_extend$nsamples,
                                           expectedRefEntitiesWithPositions=res_corr_pos_extend$mean_RegionsAtAnyPosition,
                                           expectedPositionsInAnyHotspot=res_corr_pos_extend$mean_PostionsInAnyRegion,
                                           stringsAsFactors = F)
          }
          
          # collect
          annotatedCompareEntitiesList[[CE]] <- annotatedPositions
          counts_refEntitiesWithCompEntities[1,CE] <- summaryOverlaps[1,1]
          counts_compEntitiesWithRefEntities[CE,1] <- summaryOverlaps[2,1]
          if(nsamples>0){
            expected_refEntitiesWithCompEntities[1,CE] <- expectedOverlaps[1,2]
            pvalues_refEntitiesWithCompEntities[1,CE] <- significanceTable[1,2]
            expected_compEntitiesWithRefEntities[CE,1] <- expectedOverlaps[1,3]
            pvalues_compEntitiesWithRefEntities[CE,1] <- significanceTable[1,3]
            counts_compEntitiesWithRefEntities[CE,2] <- nrow(currente)
            expected_compEntitiesWithRefEntities[CE,2] <- nrow(currente)
          }
        }else{
          # cetype=="bedRegions"
          
          res_corr_extend <- correlateBedRegions(bed_table1 = annotatedEntities,
                                                 bed_table2 = currente,
                                                 resampleBedRegions1Flag = TRUE,
                                                 resampleBedRegions2Flag = resampleCompareEntities,
                                                 resampleBedRegions1AllowOverlap = resampleReferenceEntitiesAllowOverlap,
                                                 resampleBedRegions2AllowOverlap = resampleCompareEntitiesAllowOverlap,
                                                 resampled_bed_regions1_list = resampled_referenceEntities,
                                                 genomev = genomev,
                                                 nsamples = nsamples,
                                                 altHypothesis = altHypothesis,
                                                 samplingRegions = samplingRegions,
                                                 nparallel = nparallel,
                                                 randomSeed = randomSeed,
                                                 verbose = verbose)
          # update annotated reference entities
          annotatedEntities <- res_corr_extend$annotatedBedRegions1
          annotatedEntities$classAnnotation <- NULL
          colnames(annotatedEntities)[which(colnames(annotatedEntities)=="nMatches")] <- paste0("n",CE)
          colnames(annotatedEntities)[which(colnames(annotatedEntities)=="idAnnotation")] <- CE
          # save some info
          # - annotated regions
          annotatedRegions <- res_corr_extend$annotatedBedRegions2
          annotatedRegions$classAnnotation <- NULL
          colnames(annotatedRegions)[which(colnames(annotatedRegions)=="nMatches")] <- paste0("n",referenceEntitiesName)
          colnames(annotatedRegions)[which(colnames(annotatedRegions)=="idAnnotation")] <- referenceEntitiesName
          # - summary overlaps
          summaryOverlaps <- res_corr_extend$summaryOverlaps
          # some info about significance
          if(nsamples>0){
            # - p-values
            pvalue1 <- ifelse(res_corr_extend$pvalue_Regions1overlappingAnyRegion2==0,paste0("<",1/nsamples),res_corr_extend$pvalue_Regions1overlappingAnyRegion2)
            pvalue2 <- ifelse(res_corr_extend$pvalue_Regions2overlappingAnyRegion1==0,paste0("<",1/nsamples),res_corr_extend$pvalue_Regions2overlappingAnyRegion1)
            significanceTable <- data.frame(nsamples=res_corr_extend$nsamples,
                                            pvalueOverlappingRefEntities=pvalue1,
                                            pvalueOverlappingRegions=pvalue2,
                                            stringsAsFactors = F)
            # - expected overlaps
            expectedOverlaps <- data.frame(nsamples=res_corr_extend$nsamples,
                                           expectedOverlappingRefEntities=res_corr_extend$mean_Regions1overlappingAnyRegion2,
                                           expectedOverlappingRegions=res_corr_extend$mean_Regions2overlappingAnyRegion1,
                                           stringsAsFactors = F)
          }
          
          # collect
          annotatedCompareEntitiesList[[CE]] <- annotatedRegions
          counts_refEntitiesWithCompEntities[1,CE] <- summaryOverlaps[1,1]
          counts_compEntitiesWithRefEntities[CE,1] <- summaryOverlaps[2,1]
          if(nsamples>0){
            expected_refEntitiesWithCompEntities[1,CE] <- expectedOverlaps[1,2]
            pvalues_refEntitiesWithCompEntities[1,CE] <- significanceTable[1,2]
            expected_compEntitiesWithRefEntities[CE,1] <- expectedOverlaps[1,3]
            pvalues_compEntitiesWithRefEntities[CE,1] <- significanceTable[1,3]
            counts_compEntitiesWithRefEntities[CE,2] <- nrow(currente)
            expected_compEntitiesWithRefEntities[CE,2] <- nrow(currente)
          }
        }
      }
    }
    
    
  }  
  
  # finalise return object
  returnObj <- list()
  returnObj$nsamples <- nsamples
  returnObj$annotatedEntities <- annotatedEntities
  returnObj$annotatedCompareEntitiesList <- annotatedCompareEntitiesList
  returnObj$counts_refEntitiesWithCompEntities <- counts_refEntitiesWithCompEntities
  returnObj$counts_compEntitiesWithRefEntities <- counts_compEntitiesWithRefEntities
  if(nsamples>0){
    returnObj$expected_refEntitiesWithCompEntities <- expected_refEntitiesWithCompEntities
    returnObj$pvalues_refEntitiesWithCompEntities <- pvalues_refEntitiesWithCompEntities
    returnObj$expected_compEntitiesWithRefEntities <- expected_compEntitiesWithRefEntities
    returnObj$pvalues_compEntitiesWithRefEntities <- pvalues_compEntitiesWithRefEntities
  }
  return(returnObj)
}


