

replicationTiming <- function(mutations,
                              genomev,
                              nresampling=5){
  # some checks
  supportedrefs <- c("hg19","hg38")
  if(!(genomev %in% supportedrefs)){
    message("[error replicationTiming] reference genome ",genomev," not supported. Please use one of the following: ",paste(supportedrefs,collapse = ", "))
    return(NULL)
  }
  
  # infer mutation type
  mutationType <- NULL
  requiredSmallVariantColumns <- c("chr","position")
  requiredSVbedpeColumns <- c("chrom1","start1","chrom2","start2")
  if(all(requiredSmallVariantColumns %in% colnames(mutations))){
    mutationType <- "smallVariants"
  }else if(all(requiredSVbedpeColumns %in% colnames(mutations))){
    mutationType <- "structuralVariants"
  }else{
    message("[error replicationTiming] missing required columns in mutations input table. Columns needed are \"",
            paste(requiredSmallVariantColumns,collapse = ", "),"\" for small variants, or \"",
            paste(requiredSVbedpeColumns,collapse = ", "),"\" for structural variants.")
    return(NULL)
  }
  
  # load timing table
  # timing_table_file <- paste0(rootInputTopographyLib,"/data/timing_table_annotated_",genomev,"_blacklistFiltered.tsv")
  # timing_regions <- readTable(timing_table_file)
  if(genomev=="hg19"){
    timing_regions <- timing_regions_hg19
  }else if(genomev=="hg38"){
    timing_regions <- timing_regions_hg38
  }
  timing_regions$class <- timing_regions$timinggroup
  timing_regions$chr <- timing_regions$chrom
  # get regions size
  timing_regions_size <- aggregate(timing_regions$sizeMinusN,list(TimingRegion=timing_regions$timinggroup),sum)
  timing_regions_size_vector <- timing_regions_size$x
  names(timing_regions_size_vector) <- timing_regions_size$TimingRegion
  
  positions <- mutations
  if(mutationType=="structuralVariants"){
    bp1 <- mutations[,c("chrom1","start1"),drop=F]
    bp2 <- mutations[,c("chrom2","start2"),drop=F]
    colnames(bp1) <- c("chr","position")
    colnames(bp2) <- c("chr","position")
    positions <- rbind(bp1,bp2)
  }
  
  if(!"id" %in% colnames(positions)){
    positions$id <- 1:nrow(positions)
  }
  if(!"id" %in% colnames(timing_regions)){
    timing_regions$id <- 1:nrow(timing_regions)
  }
  
  message("[info replicationTiming] calculating replication time of data")
  regionsSummary <- intersectPositionsAndBedRegions_nonOverlapping(positions = positions,
                                                                   bed_table = timing_regions,
                                                                   verbose = F)
  
  values <- regionsSummary$countsTable_regionsAtPositionClasses
  values <- values[as.character(1:10),"anyPosition"]
  normvalues <- values/timing_regions_size_vector
  normvalues <- normvalues/sum(normvalues)
  
  # now I need to resample the mutations to determine whether the bias is significant
  samplingTable <- NULL
  samplingTableNorm <- NULL
  if(nresampling>0){
    for(i in 1:nresampling){
      # i <- 1
      if(mutationType=="structuralVariants"){
        message("[info replicationTiming] random resampling of mutations ",i, " of ",nresampling)
        resampled_sv_bedpe <- resampleSV(sv_bedpe = mutations,
                                         genomev = genomev)
        # prepare the table
        bp1 <- resampled_sv_bedpe[,c("chrom1","start1"),drop=F]
        bp2 <- resampled_sv_bedpe[,c("chrom2","start2"),drop=F]
        colnames(bp1) <- c("chr","position")
        colnames(bp2) <- c("chr","position")
        tmppositions <- rbind(bp1,bp2)
        
      }else{
        # TODO resampleSNV function
        tmppositions <- NULL
      }
      
      if(!is.null(tmppositions)){
        if(!"id" %in% colnames(positions)){
          positions$id <- 1:nrow(positions)
        }
        
        # get the timing
        message("[info replicationTiming] calculating replication time of resampling ",i, " of ",nresampling)
        tmpRegionsSummary <- intersectPositionsAndBedRegions_nonOverlapping(positions = tmppositions,
                                                                            bed_table = timing_regions,
                                                                            verbose = F)

        tmpvalues <- tmpRegionsSummary$countsTable_regionsAtPositionClasses
        tmpvalues <- tmpvalues[as.character(1:10),"anyPosition"]
        tmpnormvalues <- tmpvalues/timing_regions_size_vector
        tmpnormvalues <- tmpnormvalues/sum(tmpnormvalues)
        
        samplingTable <- rbind(samplingTable,tmpvalues)
        samplingTableNorm <- rbind(samplingTableNorm,tmpnormvalues)
      }

    }
  }
  
  message("[info replicationTiming] done!")
  
  # return the results
  returnObj <- list()
  returnObj$mutationType <- mutationType
  returnObj$genomev <- genomev
  returnObj$annotatedMutations <- regionsSummary$annotatedPositions
  returnObj$countsTable <- regionsSummary$countsTable
  returnObj$timingCounts <- values
  returnObj$timingRegionsSizes <- timing_regions_size_vector
  returnObj$timingCountsNormalised <- normvalues
  if(!is.null(samplingTable)){
    returnObj$timingOfResampling <- samplingTable
    returnObj$timingOfResamplingMean <- apply(samplingTable,2,mean)
    returnObj$timingOfResamplingSD <- apply(samplingTable,2,sd)
    returnObj$timingOfResamplingSE <- returnObj$timingOfResamplingSD/sqrt(nresampling)
    returnObj$timingOfResamplingNormalised <- samplingTableNorm
    returnObj$timingOfResamplingNormalisedMean <- apply(samplingTableNorm,2,mean)
    returnObj$timingOfResamplingNormalisedSD <- apply(samplingTableNorm,2,sd)
    returnObj$timingOfResamplingNormalisedSE <- returnObj$timingOfResamplingNormalisedSD/sqrt(nresampling)
  }
  
  return(returnObj)
}


plotReplicationTiming <- function(repTimeObj,
                                  filename = NULL,
                                  main = NULL,
                                  barcolour = "darkgrey"){
  
  if(!is.null(filename)){
    cairo_pdf(filename = filename,width = 9,height = 5)
    par(mai=c(1,1,1,0.2),mfrow=c(1,2),cex=1)
  }
  # plot data
  bx <- barplot(repTimeObj$timingCountsNormalised,
                beside = T,
                names.arg = NA,
                col = barcolour,
                main = main,
                las=2,
                ylab = "normalised proportion",
                border = barcolour)
  # need to add the Early <---> Late xlabel
  ypos <- 0 - max(repTimeObj$timingCountsNormalised)*0.1
  par(xpd=TRUE)
  text(bx[1,1],ypos,labels = "Early")
  text(bx[nrow(bx),1],ypos,labels = "Late")
  text(bx[1,1]+(bx[nrow(bx),1]-bx[1,1])/2,ypos,labels="<--------------------------------->")
  par(xpd=FALSE)
  
  if(!is.null(repTimeObj$timingOfResamplingNormalisedMean)){
    # plot resampling mean
    lines(bx[,1],
          repTimeObj$timingOfResamplingNormalisedMean,
          lty=3,
          lwd=2)
    # plot resampling standard error
    tickSize <- (bx[2,1] - bx[1,1])*0.5
    for(i in 1:nrow(bx)){
      # i <- 1
      bottom <- repTimeObj$timingOfResamplingNormalisedMean[i] - repTimeObj$timingOfResamplingNormalisedSE[i]/2
      top <- repTimeObj$timingOfResamplingNormalisedMean[i] + repTimeObj$timingOfResamplingNormalisedSE[i]/2
      lines(rep(bx[i,1],2),
            c(bottom,top),
            lwd=2)
      lines(c(bx[i,1]-tickSize/2,bx[i,1]+tickSize/2),
            c(bottom,bottom),
            lwd=2)
      lines(c(bx[i,1]-tickSize/2,bx[i,1]+tickSize/2),
            c(top,top),
            lwd=2)
    }
  }
  
  
  # now again with actual counts
  bx <- barplot(repTimeObj$timingCounts,
                beside = T,
                names.arg = NA,
                col = barcolour,
                main = main,
                las=2,
                ylab = "mutation count",
                border = barcolour)
  # need to add the Early <---> Late xlabel
  ypos <- 0 - max(repTimeObj$timingCounts)*0.1
  par(xpd=TRUE)
  text(bx[1,1],ypos,labels = "Early")
  text(bx[nrow(bx),1],ypos,labels = "Late")
  text(bx[1,1]+(bx[nrow(bx),1]-bx[1,1])/2,ypos,labels="<--------------------------------->")
  par(xpd=FALSE)
  
  if(!is.null(repTimeObj$timingOfResamplingMean)){
    # plot resampling mean
    lines(bx[,1],
          repTimeObj$timingOfResamplingMean,
          lty=3,
          lwd=2)
    # plot resampling standard error
    tickSize <- (bx[2,1] - bx[1,1])*0.5
    for(i in 1:nrow(bx)){
      # i <- 1
      bottom <- repTimeObj$timingOfResamplingMean[i] - repTimeObj$timingOfResamplingSE[i]/2
      top <- repTimeObj$timingOfResamplingMean[i] + repTimeObj$timingOfResamplingSE[i]/2
      lines(rep(bx[i,1],2),
            c(bottom,top),
            lwd=2)
      lines(c(bx[i,1]-tickSize/2,bx[i,1]+tickSize/2),
            c(bottom,bottom),
            lwd=2)
      lines(c(bx[i,1]-tickSize/2,bx[i,1]+tickSize/2),
            c(top,top),
            lwd=2)
    }
  }
  
  if(!is.null(filename)) dev.off()
  
}
