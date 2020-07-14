#' @title PTC.findPP
#'
#' @description  Given a vector of miRNA names and a vector of mRNAs names,
#' creates a set of miRNAs (plausible parents) that can biologically target
#' that mRNA as predicted by TargetScan. \cr For an adequate functioning
#' miRNAs names in \code{TScan} and \code{miRs} must belong to the same version.
#' miRBase version 21 is recommended (see \link{PTC.miRv21})
#' @param TScan A \code{matrix} containing miRNA-mRNA interactions predicted by
#' TargetScan.\cr
#' If \code{NULL} (default): Loads the file \code{TScan.rda}, which
#' contains interactions from TargetScan 7.0. miRNA names in this file
#' belongs to miRBase version 21.
#' @param miRs A vector containing the names of miRNAs to be verified. These names
#' corresponds to the names of predictor candidates.
#' @param mRs A vector containing the names of mRNAs to be verified. These names
#' corresponds to the names of target genes (response variables).
#' @export
#' @seealso \link{PTC},\link{PTC.GeneSel}, \link{TScan} \link{PTC.miRv21}
#' @return A \code{list} where each element corresponds to a mRNA. Each element contains
#' the set of miRNAs that can bind the mRNA.
#'
#' @examples \dontrun{
#' data(TCGA_BRCAdata)
#' data(TScan)
#' seqData<-PTC.ptime(TCGA_BRCAdata,TCGA_BRCAdata$mRNAs[,"VIM"])
#' SelData<-PTC.GeneSel(seqData, nmiR = 30, nmR = 1500)
#' PParents<-PTC.findPP(TScan, miRs=SelData$miRs, mRs=SelData$mRs)
#' }


PTC.findPP<- function(TScan=NULL,miRs, mRs){
  if(is.null(TScan)){
    TScan<-data(TScan)
  }
  PParents<-vector("list",length = length(mRs))
  names(PParents)<-mRs
  for (i in mRs) {
    index.miRs<-which(TScan[,2]==i)
    if(length(index.miRs)>0){
      temp<-TScan[index.miRs,1]
      for (j in miRs) {
        if(any(temp==j)){
          ifelse(is.null(PParents[[i]])
                 ,PParents[[i]]<-j
                 ,PParents[[i]]<-c(PParents[[i]],j)
          )
        }
      }#end for miRs
    }#end if
  }#end for mRs
  return(PParents)
}
