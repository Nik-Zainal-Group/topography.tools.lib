
#' Sample a random position from a set of regions
#'
#' Given a table of bed regions, sample a random position.
#' 
#' @param samplingRegions data frame containing bed regions, with required columns chr, start, end, and optionally regionprob. If regionprob is missing, it will be proportional to the region size.
#' @return sampled position
#' @export
randomPositionInRegions <- function(samplingRegions){
  # check before continue
  if(nrow(samplingRegions)==0) return(NULL)
  # add regionprob if missing
  if(!"regionprob" %in% colnames(samplingRegions)){
    size <- samplingRegions$end - samplingRegions$start + 1
    samplingRegions$regionprob <- size/sum(size)
  }
  # sample the region
  sampledRegion <- sample(1:nrow(samplingRegions),size = 1,prob = samplingRegions$regionprob)
  possiblePositions <- samplingRegions$start[sampledRegion]:samplingRegions$end[sampledRegion]
  if(length(possiblePositions)==1){
    sampledPosition <- possiblePositions
  }else{
    sampledPosition <- sample(possiblePositions,size = 1)
  }
  
  return(data.frame(chr=samplingRegions$chr[sampledRegion],
                    position=sampledPosition,
                    stringsAsFactors = F))
}


#' Resample a set of positions
#'
#' Given a table of positions, resample them to obtain a set of random position.
#' The chromosome of each resampled position will be preserved. Positions will be resampled
#' within regions that are mappable and not blacklisted according to ENCODE.
#' Custom sampling regions can be specified with the parameter samplingRegions, which can
#' be used for example to assign custom probability to each region. Bear in mind that
#' for each position, the corresponding subsets of regions in the same chromosome will
#' be selected, so it would be good practice to have the regionprob for each chromosome sum to 1.
#' 
#' @param positions data.frame with positions, required columns are: chr, position
#' @param samplingRegions data frame containing bed regions, with required columns chr, start, end,
#' and optionally regionprob. If regionprob is missing, it will be proportional to the region sizes 
#' for each chromosome. If NULL, samplingRegions will be the set of regions that are mappable and 
#' not blacklisted according to ENCODE for the reference genome specified by genomev
#' @param genomev reference genome version, hg19 or hg38. This parameter will be used only if
#' samplingRegions is NULL, to load the samplingRegions for the corresponding reference genome
#' @param randomSeed set a random seed
#' @param verbose print additional output
#' @return sampled positions
#' @export
resamplePositions <- function(positions,
                              genomev=NULL,
                              samplingRegions=NULL,
                              randomSeed=NULL,
                              verbose=FALSE){
  
  # check you have regions to sample from
  if(is.null(genomev) & is.null(samplingRegions)){
    message("[error resamplePositions] both genomev and samplingRegions parameters are NULL,",
            " please specify a genome version or a samplingRegions table.")
    return(NULL)
  }
  # check if both genomev and samplingRegions have been specified and give a warning
  if(!is.null(genomev) & !is.null(samplingRegions)){
    message("[warning resamplePositions] both genomev and samplingRegions have been specified,",
            " genomev will be ignore and the samplingRegions table will be used to sample the positions.")
  }
  
  if(!is.null(randomSeed)){
    set.seed(randomSeed)
  }
  
  if(is.null(samplingRegions)){
    if(genomev=="hg19"){
      samplingRegions <- samplingRegions_hg19
    }else if(genomev=="hg38"){
      samplingRegions <- samplingRegions_hg38
    }else{
      message("[error resamplePositions] invalid genomev. Use hg19 or hg38.")
      return(NULL)
    }
    colnames(samplingRegions)[colnames(samplingRegions)=="chrregionprob"] <- "regionprob"
  }
  
  # add regionprob if missing
  if(!"regionprob" %in% colnames(samplingRegions)){
    samplingRegions$regionprob <- 0
    for(tmpchr in unique(samplingRegions$chr)){
      # tmpchr <- unique(samplingRegions$chr)[1]
      chrpos <- which(samplingRegions$chr==tmpchr)
      size <- samplingRegions$end[chrpos] - samplingRegions$start[chrpos] + 1
      samplingRegions$regionprob[chrpos] <- size/sum(size)
    }
  }
  
  resampled_positions <- NULL
  for(i in 1:nrow(positions)){
    # i <- 1
    if((i %% 50 == 0) & verbose){
      message("[info resamplePositions] processing row ",i," of ",nrow(positions))
    }
    newposition <- randomPositionInRegions(samplingRegions = samplingRegions[samplingRegions$chr==positions$chr[i],,drop=F])
    resampled_positions <- rbind(resampled_positions,newposition)
  }
  return(resampled_positions)
}



#' Resample bed regions
#'
#' Given a table of bed regions, resample them to obtain a set of random bed regions.
#' The chromosome of each resampled region will be preserved. Regions will be resampled
#' within regions that are mappable and not blacklisted according to ENCODE.
#' Custom sampling regions can be specified with the parameter samplingRegions, which can
#' be used for example to assign custom probability to each region. Bear in mind that
#' for each region, the corresponding subsets of regions in the same chromosome will
#' be selected. By default (allowRegionsOverlap=TRUE) the resampled bed regions can
#' overlap. If allowRegionsOverlap is set the FALSE, then the resampled bed regions
#' are not allowed to overlap. This is implemented by updating the samplingRegions
#' table, subtracting each resampled region. The subtraction of a resampled region
#' from samplingRegions will split or reduce the size of a sampling region, so the
#' corresponding regionprob will be split or reduced proportionally to the size reduction.
#' 
#' @param bed_table data.frame with positions, required columns are: chr, start, end
#' @param samplingRegions data frame containing bed regions, with required columns chr, start, end,
#' and optionally regionprob. If regionprob is missing, it will be proportional to the region sizes 
#' for each chromosome. If NULL, samplingRegions will be the set of regions that are mappable and 
#' not blacklisted according to ENCODE for the reference genome specified by genomev
#' @param genomev reference genome version, hg19 or hg38. This parameter will be used only if
#' samplingRegions is NULL, to load the samplingRegions for the corresponding reference genome
#' @param randomSeed set a random seed
#' @param verbose print additional output
#' @return sampled regions
#' @export
resampleBedRegions <- function(bed_table,
                               genomev=NULL,
                               samplingRegions=NULL,
                               randomSeed=NULL,
                               allowRegionsOverlap=TRUE,
                               verbose=FALSE){
  
  # check you have regions to sample from
  if(is.null(genomev) & is.null(samplingRegions)){
    message("[error resampleBedRegions] both genomev and samplingRegions parameters are NULL,",
            " please specify a genome version or a samplingRegions table.")
    return(NULL)
  }
  # check if both genomev and samplingRegions have been specified and give a warning
  if(!is.null(genomev) & !is.null(samplingRegions)){
    message("[warning resampleBedRegions] both genomev and samplingRegions have been specified,",
            " genomev will be ignore and the samplingRegions table will be used to sample the positions.")
  }
  
  if(!is.null(randomSeed)){
    set.seed(randomSeed)
  }
  
  if(is.null(samplingRegions)){
    if(genomev=="hg19"){
      samplingRegions <- samplingRegions_hg19
    }else if(genomev=="hg38"){
      samplingRegions <- samplingRegions_hg38
    }else{
      message("[error resampleBedRegions] invalid genomev. Use hg19 or hg38.")
      return(NULL)
    }
    colnames(samplingRegions)[colnames(samplingRegions)=="chrregionprob"] <- "regionprob"
  }
  
  # add regionprob if missing
  if(!"regionprob" %in% colnames(samplingRegions)){
    samplingRegions$regionprob <- 0
    for(tmpchr in unique(samplingRegions$chr)){
      # tmpchr <- unique(samplingRegions$chr)[1]
      chrpos <- which(samplingRegions$chr==tmpchr)
      size <- samplingRegions$end[chrpos] - samplingRegions$start[chrpos] + 1
      samplingRegions$regionprob[chrpos] <- size/sum(size)
    }
  }
  
  # add size if missing
  if(!"size" %in% colnames(samplingRegions)){
    samplingRegions$size <- 0
    for(tmpchr in unique(samplingRegions$chr)){
      # tmpchr <- unique(samplingRegions$chr)[1]
      chrpos <- which(samplingRegions$chr==tmpchr)
      size <- samplingRegions$end[chrpos] - samplingRegions$start[chrpos] + 1
      samplingRegions$size[chrpos] <- size
    }
  }
  
  # at this point, samplingRegions should have at least columns: chr, start, end, size, regionprob
  # get a copy that can be updated with specific column orders
  samplingRegions_copy <- samplingRegions[,c("chr","start","end","size","regionprob"),drop=F]
  
  resampled_bed_table <- NULL
  for(i in 1:nrow(bed_table)){
    # i <- 1
    if((i %% 50 == 0) & verbose){
      message("[info resampleBedRegions] processing row ",i," of ",nrow(bed_table))
    }
    
    # check the size of the region we need to resample
    rsize <- abs(bed_table$end[i] - bed_table$start[i]) + 1
    # check regions that might contain it in the same chromosome
    tmpSamplingRegions <- samplingRegions_copy[samplingRegions_copy$chr==bed_table$chr[i] & samplingRegions_copy$size>=rsize,,drop=F]
    if(nrow(tmpSamplingRegions)==0){
      # we have no compatible regions to sample from
      message("[warning resampleBedRegions] could not resample region in bed_table line ",i,
              ", either chromosome not availalbe or no regions large enough in samplingRegions. ",
              "If allowRegionsOverlap=FALSE, this may also depend on other resampled regions on the same chromosome.")
      newrow <- NULL
    }else{
      # scale the probability according to the available regions
      nonzeroprob <- tmpSamplingRegions$regionprob > 0
      tmpSamplingRegions$regionprob[nonzeroprob] <- tmpSamplingRegions$regionprob[nonzeroprob]/sum(tmpSamplingRegions$regionprob[nonzeroprob])
      tmpSamplingRegions$end <- tmpSamplingRegions$end - rsize + 1
      # get a random position as the start
      rstart <- randomPositionInRegions(samplingRegions = tmpSamplingRegions)
      rend <- rstart$position + rsize - 1
      newrow <- data.frame(chr = bed_table$chr[i],
                           start = rstart$position,
                           end = rend,
                           stringsAsFactors = F)
    }
    
    resampled_bed_table <- rbind(resampled_bed_table,newrow)
    
    # we now need to remove the new row from the samplingRegions
    if(!allowRegionsOverlap & !is.null(newrow)){
      # add ids
      newrow$id <- 1
      samplingRegions_copy$id <- 1:nrow(samplingRegions_copy)
      # we need to replace the samplingRegions for the current chromosome only
      # by definition, the sampled segment has to be included in one of the
      # samplingRegions segments
      res_int <- intersectBed_nonOverlapping(bed_table1 = samplingRegions_copy[samplingRegions_copy$chr==bed_table$chr[i],,drop=F],bed_table2 = newrow)
      intersectTable <- res_int$intersectionTable
      samplingRegion_id <- intersectTable$segment1[!is.na(intersectTable$segment2)]
      intersectTable <- intersectTable[intersectTable$segment1==samplingRegion_id,]
      # at this point, intersectTable contains 1,2 or 3 regions, depending on
      # whether the sampled regions stretches the whole region, is at the start/end, or
      # it is in the middle. Remove the shared region
      intersectTable <- intersectTable[intersectTable$bedtable!="shared",]
      # we can now remove the sampling region from samplingRegions and replace if there
      # are 1 or 2 segments left, redistribution the probability proportionally to the 
      # size of the regions left
      oldProb <- samplingRegions_copy$regionprob[as.integer(samplingRegion_id)]
      oldSize <- samplingRegions_copy$size[as.integer(samplingRegion_id)]
      samplingRegions_copy <- samplingRegions_copy[-as.integer(samplingRegion_id),,drop=F]
      if(nrow(intersectTable)>0){
        samplingRegions_copy$id <- NULL
        newsize <- intersectTable$end-intersectTable$start+1
        samplingRegions_copy <- rbind(samplingRegions_copy,
                                      data.frame(chr=intersectTable$chr,
                                                 start=intersectTable$start,
                                                 end=intersectTable$end,
                                                 size=newsize,
                                                 regionprob=oldProb*newsize/oldSize,
                                                 stringsAsFactors = TRUE))
        # sorting samplingRegions_copy should not be necessary
      }
      
    }

  }
  
  return(sortBed(resampled_bed_table))
}
