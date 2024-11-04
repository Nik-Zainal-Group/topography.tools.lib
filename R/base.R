
#' @importFrom foreach %dopar%
NULL

#' Sort chromosomes
#'
#' Order a list of chromosome names.
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

