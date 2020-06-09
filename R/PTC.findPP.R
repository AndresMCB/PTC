#' @title PTC.findPP
#'
#' @description  For each mRNA, creates a set of miRNAs (plausible parents) that can biologically target that mRNA.
#' A total of \code{nmiR} miRNAs and \code{nmR} mRNAs with the largest gene expression Median Absolute Deviation (MAD)
#' are selected as predictor candidates and target mRNAs respectively.
#' \cr The set of plausible parents of each gene is a subset of the predictor candidates and contains those miRNA that can bind
#' the gene as predicted by TargetScan 7.0.
#' miRNAs names are converted to \code{miRBase v.21} during selection process.
#'
#' @param TScan A \code{matrix} containing miRNA-mRNA interactions predicted by
#' TargetScan 7.0. The first and second columns contain miRNAs names from
#' miRBase versions 18 and 17 respectively. The third column contains the names of the target mRNAs.
#' @param miRs A vector containing the names of miRNAs to be verified. These names
#' corresponds to the names of predictor candidates.
#' @param mRs A vector containing the names of mRNAs to be verified. These names
#' corresponds to the names of target genes (response variables).
#' @export
#' @seealso \link{PTC},\link{PTC.GeneSel}, \link{TScan}
#' @return A list containing the set of miRNAs that can bind each target mRNA:
#'
#' @examples \dontrun{
#' data(TCGA_BRCAdata)
#' data(TScan)
#' seqData<-PTC.ptime(TCGA_BRCAdata,TCGA_BRCAdata$mRNAs[,"VIM"])
#' SelData<-PTC.GeneSel(seqData, topk_miR = 30, topk_mR = 1500)
#' PParents<-PTC.findPP(TScan, miRs=SelData$miRs, mRs=SelData$mRs)
#' }


PTC.findPP<- function(TScan,miRs, mRs){

  PParents<-vector("list",length = length(mRs))
  names(PParents)<-mRs
  for (i in mRs) {
    index.miRs<-which(TScan[,"mRNA"]==i)
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
