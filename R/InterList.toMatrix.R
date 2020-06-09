#' @title InterList.toMatrix
#'
#' @description Given a \code{list} of mRNAs containing the set of their
#' inferred parents (the output of \link{Extract.Parents}), it returns
#' a \code{matrix} with the miRNA-mRNA interactions. Each row represents an interaction.
#' The first column contains miRNAs. The second column contains mRNAs
#' @inheritParams Find.Confirmed
#' @export
#' @seealso
#' @return A matrix with 2 columns containing the miRNA-mRNA interactions.
#' @examples \dontrun{
#' data(TCGA_BRCAdata)
#' data(TScan)
#' seqData<-PTC.ptime(TCGA_BRCAdata,TCGA_BRCAdata$mRNAs[,"VIM"])
#' SelData<-PTC.GeneSel(seqData, topk_miR = 30, topk_mR = 1500)
#' PParents<-PTC.findPP(TScan, miRs=SelData$miRs, mRs="SORCS1")
#' temp <-SelData$PParents[["SORCS1"]]
#' SORCS1.X.TScan=SelData$d[,temp]
#' SORCS1.Parents<-PTC.TestInvariance(Y=SelData$d[,"SORCS1"], X=SORCS1.X.TScan)
#' SORCS1.Parents<-Extract.Parents(SORCS1.Parents,SelData$PParents[["SORCS1"]])
#' SORCS1.relMatrix<-na.omit(InterList.toMatrix(SORCS1.Parents$Plist$Names))
#' }

InterList.toMatrix<-function(miRNAsList){
  InterMatrix<-matrix(nrow = length(unlist(miRNAsList)),ncol = 2)
  k<-1
  if(length(miRNAsList)>0){
    for (i in 1:length(miRNAsList)) {
      for (j in 1:length(miRNAsList[[i]])) {
        InterMatrix[k,1]<-miRNAsList[[i]][j]
        InterMatrix[k,2]<-names(miRNAsList[i])
        k<-k+1
      }
    }
  }
  colnames(InterMatrix)<-c("miR","mR")
  return(InterMatrix)
}
