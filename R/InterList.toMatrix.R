#' @title InterList.toMatrix
#'
#' @description Given a \code{list} where each element is a mRNA containing the set of
#' the miRNAs (parents) such mRNA can interact with(e.g. the output of \link{Extract.Parents}),
#' it returns a \code{matrix} with the miRNA-mRNA interactions. Each row represents an interaction.
#' The first column contains miRNAs. The second column contains mRNAs
#' @inheritParams Confirmed.fromList
#' @export
#' @seealso
#' @return A matrix with 2 columns containing the miRNA-mRNA interactions.
#' @examples \dontrun{
#' data(TCGA_BRCAdata)
#' data(TScan)
#' seqData<-PTC.ptime(TCGA_BRCAdata,TCGA_BRCAdata$mRNAs[,"VIM"])
#' SelData<-PTC.GeneSel(seqData, nmiR = 30, nmR = 1500)
#' PParents<-PTC.findPP(TScan, miRs=SelData$miRs, mRs="ONECUT2")
#' temp <-SelData$PParents[["ONECUT2"]]
#' ONECUT2.X.TScan=SelData$d[,temp]
#' set.seed(1)
#' ONECUT2.Parents<-PTC.TestInvariance(Y=SelData$d[,"ONECUT2", drop = F], X=ONECUT2.X.TScan)
#' ONECUT2.Parents<-Extract.Parents(ONECUT2.Parents,SelData$PParents["ONECUT2"])
#' ONECUT2.relMatrix<-InterList.toMatrix(ONECUT2.Parents$Names)
#' }

InterList.toMatrix<-function(InterList){
  if(length(which(is.na(InterList)))>0){
    InterList<-InterList[-which(is.na(InterList))]
  }
  index<-which(sapply(InterList, is.null))
  if(length(index>0)){
    InterList<-InterList[-index]
  }
  InterMatrix<-matrix(nrow = length(unlist(InterList)),ncol = 2)
  k<-1
  if(length(InterList)>0){
    for (i in 1:length(InterList)) {
      for (j in 1:length(InterList[[i]])) {
        InterMatrix[k,1]<-InterList[[i]][j]
        InterMatrix[k,2]<-names(InterList[i])
        k<-k+1
      }
    }
  }
  colnames(InterMatrix)<-c("miR","mR")
  return(InterMatrix)
}
