
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
  sampledPosition <- sample(samplingRegions$start[sampledRegion]:samplingRegions$end[sampledRegion],size = 1)
  return(data.frame(chr=samplingRegions$chr[sampledRegion],
                    position=sampledPosition,
                    stringsAsFactors = F))
}


#' Resample a set of positions
#'
#' Given a table of positions, resample them to obtain a set of random position.
#' The chromosome of each resampled position will be preserved. Positions will be resampled
#' within regions that are mappable and not blacklisted according to ENCODE.
#' Custom sampling regions can be specified by the parameter samplingRegions, which can
#' be used for example to assign custom probability to each region. Bear in mind that
#' for each position, the correspoinding subsets of regions in the same chromosome will
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
                              genomev,
                              samplingRegions=NULL,
                              randomSeed=NULL,
                              verbose=FALSE){
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

