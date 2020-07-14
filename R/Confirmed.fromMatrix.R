#' @title Confirmed.fromMatrix
#'
#' @description Given a matrix with miRNA-mRNA interactions, returns a \code{matrix} with all
#' of those interactions that are present in the \code{GroundT}.
#' @inheritParams Confirmed.fromList
#' @param InterMatrix A \code{matrix} containing miRNA-mRNA interactions.
#' The first columns corresponds to miRNAs. The second column corresponds to mRNA.
#' @export
#' @seealso \link{PTC}, \link{GroundT} \link{Confirmed.fromList}
#'  \link{InterList.toMatrix}
#' @return A \code{matrix} containing all interactions in \code{InterMatrix} that
#'  are in the \code{GroundT}
#' @examples \dontrun{
#'data(TCGA_BRCAdata)
#'data(GroundT)
#'test1<-PTC(miRNAs=TCGA_BRCAdata$miRs,mRNAs=TCGA_BRCAdata$mRNAs, VIM=TCGA_BRCAdata$mRNAs[,"VIM"])
#'aux<-InterList.toMatrix(test1$Names)
#'t1.Confirmed<-Confirmed.fromMatrix(aux,GroundT)
#' }


Confirmed.fromMatrix<-function(InterMatrix,GroundT){
  Confirmed<-matrix(data = NA,ncol = 2)
  colnames(Confirmed)<-c("miR","mR")
  aux<-1
  miRs<-unique(InterMatrix[,1])
  for (i in miRs) {
    indexes<-which(InterMatrix[,1]==i)
    mRs<-InterMatrix[indexes,2]
    for(j in mRs){
      if(any(GroundT[[j]]==i)){
        if(aux==1){Confirmed[aux,]<-c(i,j)}
        else{
          Confirmed<-rbind(Confirmed,c(i,j))
        }
        aux<-aux+1
      }
    }
  }#end for miRs
  return(Confirmed)
}
