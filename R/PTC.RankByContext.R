#' @title PTC.RankByTScontext
#' @description This fuction uses TS70 Conserved_Site_Context_Scores (from
#' \url{http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_70}
#' to rank a matrix of inferred miRNA-mRNA relationships.
#' @param TS70context: A \code{matrix} containing conserved site context++ scores. The headers
#' of the matrix from TargetScan must be conserved when the file is imported into R.
#' @param relationshipsM: A matrix with two columns containing miRNAs and mRNAs names.
#' First column represents miRNAs and second column mRNAs (e.g the output of \link{InterList.toMatrix}).
#' @return A relationship matrix with four columns [Rank, miRNA, mRNA, Score].
#' The more negative the score, the best the rank of the corresponding miRNA-mRNA pair.
#' @seealso \link[PTC]{PTC}
#' @examples \dontrun{
#' data(TCGA_BRCAdata)
#' data(TScan)
#' data(TS7.0_Conserved_Site_Context_Scores)
#' seqData<-PTC.ptime(TCGA_BRCAdata,TCGA_BRCAdata$mRNAs[,"VIM"])
#' SelData<-PTC.GeneSel(seqData, topk_miR = 30, topk_mR = 1500)
#' PParents<-PTC.findPP(TScan, miRs=SelData$miRs, mRs="SORCS1")
#' temp <-SelData$PParents[["SORCS1"]]
#' SORCS1.X.TScan=SelData$d[,temp]
#' SORCS1.Parents<-PTC.TestInvariance(Y=SelData$d[,"SORCS1"], X=SORCS1.X.TScan)
#' SORCS1.Parents<-Extract.Parents(SORCS1.Parents,SelData$PParents[["SORCS1"]])
#' SORCS1.relMatrix<-na.omit(InterList.toMatrix(SORCS1.Parents$Plist$Names))
#' SORCS1.Summary<-PTC.RankByContext(TS7.0_Conserved_Site_Context_Scores
#'                                   ,SORCS1.relMatrix)
#' }

PTC.RankByContext <-function(TS70context,relationshipsM){

  if(ncol(relationshipsM)==2){
    names(relationshipsM)<-c("miR","mR")
  }else{
    stop("relationshipsM has to be a matrix with 2 columns [miR,mR]")
  }

  TS70context[,"miRNA"]<-PTC.miRv21(TS70context[,"miRNA"])

  #--------------------Scores
  indexes<-apply(relationshipsM,1,FUN = function(x,y) {
    output<-which(y[,"miRNA"]==x['miR']
                  &y[,"Gene.Symbol"]==x['mR'])
    return(output)
  }
  ,y=TS70context)

  #--------calculate Score as the average of context scores of each pair
  Score<-matrix(nrow = nrow(relationshipsM),ncol=1)
  colnames(Score)<-"Score"
  Score<-sapply(indexes,FUN = function(x,y) {
                  output<-mean(y[x,"context...score"])
                  return(output)
                  }
                ,y=TS70context)
  relationshipsM<-cbind(relationshipsM,Score)
  relationshipsM<-relationshipsM[order(relationshipsM[,'Score']
                                       ,decreasing = T),]
  relationshipsM<-cbind(1:nrow(relationshipsM),relationshipsM)
  colnames(relationshipsM)[1]<-"rank"
  rownames(relationshipsM)<-1:nrow(relationshipsM)
  return(relationshipsM)
}
