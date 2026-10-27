#' @title PTC.RankByContext
#' @description This function uses TS70 Conserved_Site_Context_Scores (from
#' \url{http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_70})
#' to rank a matrix of inferred miRNA-mRNA relationships.
#' @param TS70context A \code{matrix} containing conserved site context++ scores. The headers
#' of the matrix from TargetScan must be conserved when the file is imported into R.
#' @param relationshipsM A matrix with two columns containing miRNAs and mRNAs names.
#' First column represents miRNAs and second column mRNAs (e.g the output of \link{InterList.toMatrix}).
#' @export
#' @return A relationship matrix with four columns [Rank, miRNA, mRNA, Score].
#' The more negative the score, the better the rank of the corresponding miRNA-mRNA pair.
#' Pairs without a context score (\code{NaN}) are ranked last.
#' @seealso \link[PTC]{PTC}
#' @examples \dontrun{
#' data(TCGA_BRCAdata)
#' data(TScan)
#' data(TS7.0_Conserved_Site_Context_Scores)
#'
#' seqData<-PTC.ptime(TCGA_BRCAdata,TCGA_BRCAdata$mRNAs[,"VIM"])
#' SelData<-PTC.GeneSel(seqData, nmiR = 30, nmR = 1500)
#' PParents<-PTC.findPP(TScan, miRs=SelData$miRs, mRs="ONECUT2")
#' temp <-SelData$PParents[["ONECUT2"]]
#'
#' ONECUT2.X.TScan=SelData$d[,temp]
#' set.seed(1)
#' ONECUT2.Parents<-PTC.TestInvariance(Y=SelData$d[,"ONECUT2", drop =F], X=ONECUT2.X.TScan)
#' ONECUT2.Parents<-Extract.Parents(ONECUT2.Parents,SelData$PParents["ONECUT2"])
#' ONECUT2.relMatrix<-InterList.toMatrix(ONECUT2.Parents$Names)
#' ONECUT2.Summary<-PTC.RankByContext(TS7.0_Conserved_Site_Context_Scores
#'                                   ,ONECUT2.relMatrix)
#' }

PTC.RankByContext <-function(TS70context,relationshipsM){

  if(ncol(relationshipsM)==2 && nrow(relationshipsM)==0){
    # No inferred interactions: return an empty table with the usual columns
    return(cbind(rank=character(0), relationshipsM, Score=character(0)))
  }
  if(ncol(relationshipsM)==2){
    names(relationshipsM)<-c("miR","mR")
  }else{
    stop("relationshipsM has to be a matrix with 2 columns [miR,mR]")
  }

  # Only the TargetScan rows of the genes being ranked are needed, so keep those
  # before converting miRNA names (much faster than converting all 1.4 million rows)
  TS70context<-TS70context[TS70context[,"Gene.Symbol"] %in% relationshipsM[,2],,drop=FALSE]
  TS70context[,"miRNA"]<-PTC.miRv21(TS70context[,"miRNA"])

  #--------------------Scores
  indexes<-apply(relationshipsM,1,FUN = function(x,y) {
    output<-which(y[,"miRNA"]==x['miR']
                  &y[,"Gene.Symbol"]==x['mR'])
    return(output)
  }
  ,y=TS70context)

  #--------calculate Score as the average of context scores of each pair
  Score<-sapply(indexes,FUN = function(x,y) {
                  output<-mean(y[x,"context...score"])
                  return(output)
                  }
                ,y=TS70context)
  relationshipsM<-cbind(relationshipsM,Score)
  # Score is stored as text in the matrix, so sort on its numeric value:
  # most negative first, pairs without a context score (NaN) last
  relationshipsM<-relationshipsM[order(as.numeric(relationshipsM[,'Score'])),,drop=FALSE]
  relationshipsM<-cbind(1:nrow(relationshipsM),relationshipsM)
  colnames(relationshipsM)[1]<-"rank"
  rownames(relationshipsM)<-1:nrow(relationshipsM)
  return(relationshipsM)
}
