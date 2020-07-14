#' @title PTC.TestInvariance
#' @importFrom seqICP seqICP.s
#' @description Finds a set of miRNAs (predictors) that
#' are causal parents of a target gene. This set is defined as the union of sets
#' that are invariant. Invariance is determined by using the \code{decoupled test}
#' from [seqICP]{seqICP.s}.
#' @inheritParams PTC
#' @param Y [nx1] A \code{named vector} containing the sequential gene expression of a
#' target gene Y.\cr Please use \code{drop = FALSE} when assigning a mRNA gene expression from a \code{matrix}
#' to preserve the name of the mRNA.
#' @param X [nxp] A \code{named matrix} containing the sequential gene expression of
#' \code{p} plausible parents (miRNAs that can bind the target mRNA). \cr
#' for a correct functioning, column names must be the mRNA names.
#' @export
#'
#' @seealso \link{PTC}, \link{seqICP::seqICP.s}
#' @return
#' \item{\code{Parents}}{A set containing the indexes of the parents inferred by PTC.
#' These indexes correspond to the indexes of the miRNAs in the set \code{PParents}
#' (Plausible Parents) obtained from \link{PTC.GeneSel}}
#' @examples \dontrun{
#' data(TCGA_BRCAdata)
#' data(TScan)
#' seqData<-PTC.ptime(TCGA_BRCAdata,TCGA_BRCAdata$mRNAs[,"VIM"])
#' SelData<-PTC.GeneSel(seqData, nmiR = 30, nmR = 1500)
#' PParents<-PTC.findPP(TScan, miRs=SelData$miRs, mRs="ONECUT2")
#' temp <-SelData$PParents[["ONECUT2"]]
#' ONECUT2.X.TScan=SelData$d[,temp, drop = FALSE]
#' set.seed(1)
#' ONECUT2.Parents<-PTC.TestInvariance(Y=SelData$d[,"ONECUT2", drop = F], X=ONECUT2.X.TScan)
#' }
#' @references
#' A Pseudo-Temporal Causality Approach to Identifying miRNA-mRNA
#' Interactions During Biological Processes\cr
#' Andres M. Cifuentes-Bernal, Vu VH Pham, Xiaomei Li,
#' Lin Liu, Jiuyong Li, Thuc Duy Le \cr
#' bioRxiv 2020.07.07.192724;
#'  \url{https://doi.org/10.1101/2020.07.07.192724}
#'

PTC.TestInvariance<-function(Y, X, ngrid=2, alpha=0.02
                             , explore.all=TRUE, silent=TRUE
                             , complements = TRUE){

  test = "decoupled"
  X<-as.matrix(X)
  par.test = list(grid = as.integer(seq(0,nrow(X),length.out=ngrid+1))
                  ,complements = complements
                  ,link = sum
                  ,alpha = alpha
                  ,B =100)
  model="iid"
  par.model=list(pknown=FALSE, p=0, max.p=10)

  finish<-FALSE
  p<-ncol(X)

  S.Union<-numeric(0)
  aux<-p
  model.rejected<-TRUE
  if(is.list(Y)){
    Y<-as.matrix(unlist(Y,use.names = TRUE))
  }

  while(!finish)
  {
    Sets <- combn(x=1:p,aux)

    #Remove all subsets in current UNION to don't evaluate them
    if(length(S.Union)>0){
      sets.to.remove<-combn(S.Union,aux)
      message(paste(" (PTC) nSets to REMOVE in current level = {"
                    ,toString(ncol(sets.to.remove)),"}",sep=" "))
      for (i in 1:ncol(sets.to.remove)) {
        #intersect each col in Set with set.to.remove
        temp<-apply(Sets, 2,"intersect",sets.to.remove[,i])
        #find where in Sets is current set.to.remove
        temp<-sapply(temp, "length")
        index.temp<-which(temp==aux)
        #remove that set
        Sets<-Sets[,-index.temp,drop=FALSE]
      }
      rm(temp,index.temp,sets.to.remove)
    }


    nSets <- ncol(Sets)
    message(paste(" (PTC) nSets to explore in current level = {",toString(nSets),"}",sep=" "))
    for (j in 1:nSets) {

      if(!silent){
        message(paste(" Current evaluated set  = {",toString(Sets[,j]),"}",sep=" "))
      }

      res <- seqICP.s(X, Y, Sets[,j], test, par.test, model, par.model)
      if(res$p.value>par.test$alpha){
        model.rejected<-FALSE
        S.Union<-union(S.Union,Sets[,j])
        if(!silent){
          message(paste(" Accepted set  = {",toString(Sets[,j]),"}",sep=" "))
        }
        # if 1 set was "accepted", and dont want to explore all, finish
        if(explore.all==FALSE){
          finish=TRUE
          break
        }
      }
      if(length(S.Union)==p){
        finish<-TRUE
        break
      }
    }#end for
    aux<-aux-1
    if(aux==0)
    {
      finish=TRUE
    }
  }#end while

  if(!model.rejected){
    message(paste(colnames(Y),
                  " PTC.findParents  = {",toString(S.Union),"}",sep=" "))
  } else{
    S.Union <- numeric(0)
    message("----------------------PTC.findParents MODEL REJECTED------------------------")
  }
  Parents<-list(S.Union)

  names(Parents)<-colnames(Y)

  return(Parents)
}# end PTC.TestInvariance
