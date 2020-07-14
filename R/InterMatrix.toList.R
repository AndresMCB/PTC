#' @title InterMatrix.toList
#'
#' @description Given a matrix containing miRNA-mRNA interactions,
#' creates a \code{list} where each element corresponds to a mRNA and contains
#' the miRNAs related to that mRNA.
#'
#' @param InterMatrix A two column \code{matrix} containing miRNA-mRNA relationships.
#' The first column corresponds to miRNAs.
#' the second column corresponds to mRNAs.
#' @export
#' @seealso
#' @return A \code{list} where each element has the following structure:
#'    \item{\code{name:}}{ mRNA name.}
#'    \item{\code{data:}}{ Names of the miRNAs linked to the current mRNA}
#'
#'
#' @examples \dontrun{
#' data(TCGA_BRCAdata)
#' data(TScan)
#' VIM=TCGA_BRCAdata$mRNAs[,"VIM"]
#' GEData<-list(TCGA_BRCAdata$miRs,TCGA_BRCAdata$mRNAs)
#' names(GEData)<-c("miRs","mRNAs")
#' seqData<-PTC.ptime(GEData,VIM)
#' SelData<-PTC.GeneSel(seqData, nmiR = 30, nmR = 1500)
#' InteractionsM<-InterList.toMatrix(SelData$PParents)
#' InteractionsL<-InterMatrix.toList(PParentsM)
#' }

InterMatrix.toList<-function(InterMatrix){
  InterMatrix<-na.omit(InterMatrix)
  aux<-unique(InterMatrix[,2])
  InterList<-vector("list",length = length(aux))
  names(InterList)<-aux
  for (i in aux) {
    index<-which(InterMatrix[,2]==i)
    InterList[[i]]<-InterMatrix[index,1,drop=F]
  }
  return(InterList)
}
