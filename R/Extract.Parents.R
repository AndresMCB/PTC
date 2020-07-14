#' @title Extract.Parents
#'
#' @description Uses the outcomes of \link{PTC.GeneSel} and \link{PTC.TestInvariance}
#' to create a \code{list} of all mRNAs in \code{Predictors} with at least one causal parent.
#' @param PTC.outcome A \code{list} containing the results of the function \link{PTC.TestInvariance}.
#' Each element of the \code{list} must have the following form:
#' \itemize{
#'    \item{\code{name:}}{ mRNA name}
#'    \item{\code{data:}}{ indexes indicating the position of the inferred miRNAs in the vector of
#'    plausible parents.
#'    }
#'  }
#' @param Predictors  A \code{list} containing the set of miRNAs to be used as predictors
#'    (plausible parents) of each mRNA, obtained from \link{PTC.GeneSel}.
#' @export
#'
#' @seealso
#' @return A list of 3 elements.
#'    \item{Index}{ Indexes of miRNAs inferrred as causal parents of each mRNA.
#'    The indexes correspond to miRNAs position in the set of plausible parents of each mRNA.}
#'    \item{Names}{ A list where each element corresponds to a mRNA and contains
#'    the names of the miRNAs that are inferred as causal parents of that mRNA.}
#'    \item{genes}{ Names of all mRNAs with at least one parent inferred by PTC.}
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
#' }

Extract.Parents<- function(PTC.outcome,Predictors){
  PP<-lapply(PTC.outcome,
             FUN = function(x) {
               output<-length(x)>0
               return(output)
             })

  PP<-unlist(PP,recursive = FALSE)
  PP.index<-unlist(which(PP),use.names = FALSE)
  PP<-PTC.outcome[PP.index, drop=FALSE]

  Names<-vector("list",length = length(PP))
  names(Names)<-names(PP)
  for (i in names(PP)) {
    Names[[i]]<-Predictors[[i]][PP[[i]]]
  }
  return(list(Index=PP
              ,Names=na.omit(Names)
              ,genes=names(PP))
  )
}
