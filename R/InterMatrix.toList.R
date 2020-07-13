#' @title InterMatrix.toList
#'
#' @description Given a matrix containing miRNA-mRNA pairs,
#' creates a /code{list} where each element corresponds to a mRNA and contains
#' all miRNA linked to that mRNA.
#'
#' @param InterMatrix a two column \code{matrix}. The first column corresponds to miRNA.
#' the second column corresponds to mRNA.
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
#' PParentsM<-InterList.toMatrix(na.omit(SelData$PParents))
#' PParentsL<-InterMatrix.toList(PParentsM)
#' PParentsL2<-SelData$PParents[-which(sapply(SelData$PParents, is.null))]
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
