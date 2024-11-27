

# positions columns needed: chr, position, id
# bed_table columns needed: chr, start, end, optionally id
# extended indicates how many bps the bed regions should be extended
# left and right to include very close positions that do not strictly in the regions
# bed_table <- data.frame(chr=c(1,1,1,1,1),id=c(1,2,3,4,5),start=c(1,2,3,5,8),end=c(1,2,3,5,8),stringsAsFactors = F)
# positions <- data.frame(chr=c(1,1),position=c(3,4),id=c(1,2),stringsAsFactors = F)
# bed_table <- data.frame(chr=c(1,1,1,1,1),id=c(1,2,3,4,5),class=c(1,2,3,4,5),start=c(1,2,3,5,8),end=c(1,2,3,5,8),stringsAsFactors = F)
# positions <- data.frame(chr=c(1,1),position=c(3,4),id=c(1,2),class=c(1,2),stringsAsFactors = F)
# extended <- 1
# nsamples <- 1

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
  
  # check which chromosomes have overlap if any
  isBedRegionOverlap <- FALSE
  overlapChroms <- checkBedRegionsOverlap(bed_table)
  if(!is.null(overlapChroms)) isBedRegionOverlap <- TRUE
  
  # sort
  positions <- sortPositions(positions)
  bed_table <- sortBed(bed_table)
  # give each region its own class id
  # bed_table$class=paste0(bed_table$chr,"_",sprintf("%d",bed_table$start),"_",sprintf("%d",bed_table$end))
  # if(!("id" %in% colnames(bed_table))){
  #   bed_table$id <- 1:nrow(bed_table)
  # }
  # might need a copy for simulations later
  # bed_table_copy <- bed_table
  
  # some checks
  # if(!extendedBedRegionAllowOverlap & isBedRegionOverlap){
  #   message("[warning correlatePositionsWithBedRegions] extendedBedRegionAllowOverlap ",
  #           "set to FALSE but will be ignored because bed_table contains overlaps already.")
  #   extendedBedRegionAllowOverlap <- TRUE
  # }
  # if(!extendedBedRegionAllowOverlap & resampleBedRegionsAllowOverlap){
  #   message("[warning correlatePositionsWithBedRegions] extendedBedRegionAllowOverlap ",
  #           "set to FALSE but will be ignored because resampleBedRegionsAllowOverlap is TRUE.")
  #   extendedBedRegionAllowOverlap <- TRUE
  # }
  
  # extend regions if required
  # if(extended>0){
  #   if(!extendedBedRegionAllowOverlap){
  #     tmp_bed_table <- extendBedRegionsWithNoOverlap(bed_table = bed_table,
  #                                                    extended = extended)
  #   }else{
  #     tmp_bed_table <- extendBedRegions(bed_table = bed_table,
  #                                       extended = extended)
  #   }
  # }else{
  #   tmp_bed_table <- bed_table
  # }

  
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
  
  # add some annotations
  # bed_table$npositionsInRegion <- regionsCounts[bed_table$id]
  # bed_table$positionIds <- NA
  # for(i in 1:nrow(bed_table)){
  #   # i <- 1
  #   posids <- res_assign$annotatedPositions[res_assign$annotatedPositions$regionClassAnnotation==bed_table[i,"id"],"id"]
  #   if(length(posids)>0) bed_table[i,"positionIds"] <- paste(posids,collapse = ";")
  # }
  # # change a column name to match bed table
  # res_assign$annotatedPositions$regionid <- res_assign$annotatedPositions$regionClassAnnotation
  # res_assign$annotatedPositions$regionClassAnnotation <- NULL
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
    # for (i in 1:nsamples){
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
        # # extend if necessary
        # if(!resampleBedRegionsAllowOverlap & !extendedBedRegionAllowOverlap){
        #   resampled_bed_table <- extendBedRegionsWithNoOverlap(bed_table = resampled_bed_table,
        #                                                        extended = extended)
        # }else{
        #   resampled_bed_table <- extendBedRegions(bed_table = resampled_bed_table,
        #                                           extended = extended)
        # }

      }else{
        resampled_bed_table <- bed_table
      }
      
      # now get the stats
      res_sample_assign <- intersectPositionsAndBedRegions(positions = resampled_positions,
                                                           bed_table = resampled_bed_table)
      # 
      # 
      # regionsCounts_sample <- apply(res_sample_assign$countsTable, 2, function(x) sum(x))
      # nregionsWithPositions_sample <- sum(regionsCounts_sample[1:(length(regionsCounts_sample)-1)]>0)
      # # sampledNregionsWithPositions <- c(sampledNregionsWithPositions,nregionsWithPositions_sample)
      # npositionsInRegions_sample <- sum(regionsCounts_sample[1:(length(regionsCounts_sample)-1)])
      # # sampledNpositionsInRegions <- c(sampledNpositionsInRegions,npositionsInRegions_sample)
      # return(c(nregionsWithPositions_sample,npositionsInRegions_sample))
      returnObj <- list()
      returnObj$sampled_PostionsInAnyRegion <- res_sample_assign$totalPostionsInAnyRegion
      returnObj$sampled_RegionsAtAnyPosition <- res_sample_assign$totalRegionsAtAnyPosition
      returnObj$sampled_positionsInRegionClasses <- as.matrix(res_sample_assign$countsTable_positionsInRegionClasses[rownames(countsTable_positionsInRegionClasses),colnames(countsTable_positionsInRegionClasses),drop=F])
      returnObj$sampled_regionsAtPositionClasses <- as.matrix(res_sample_assign$countsTable_regionsAtPositionClasses[rownames(countsTable_regionsAtPositionClasses),colnames(countsTable_regionsAtPositionClasses),drop=F])
      return(returnObj)
    }
    
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
    
  }
  return(returnObj)
}




