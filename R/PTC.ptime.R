#' @title PTC.ptime
#' @inheritParams PTC
#' @description Orders the gene expression by VIM_Time order.
#' @param matchedData A \code{List} with two elements. A miRNA gene expressions matrix (1st element)
#' and A mRNA gene expressions matrix (2nd element).
#' Columns represent miRNAs and rows represent samples.
#' @export
#' @seealso \link[PTC]{PTC}
#' @return Pseudotime ordered matched data.
#' @examples \dontrun{
#'   data(TCGA_BRCAdata)
#'   Time_series<-PTC.ptime(TCGA_BRCAdata, TCGA_BRCAdata$mRNAs[,"VIM"])
#' }

PTC.ptime<-function(matchedData, VIM){

  VIM_Time<-order(VIM)
  matchedData[[1]]<-matchedData[[1]][VIM_Time,,drop=FALSE]
  matchedData[[2]]<-matchedData[[2]][VIM_Time,,drop=FALSE]

  return(matchedData)

}
