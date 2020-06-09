#' @title Find.Confirmed
#'
#' @description Returns a \code{list} with all experimentally confirmed interactions
#' inferred by PTC.
#' @param miRNAsList A \code{list} where each element has the following structure:
#'  \itemize{
#'    \item{\code{name:}}{ mRNA name.}
#'    \item{\code{data:}}{ Names of the miRNAs inferred as parents by PTC.}
#'  }
#' @param GroundT  A \code{list} where each element has the following structure:
#'  \itemize{
#'    \item{\code{name:}}{ mRNA name.}
#'    \item{\code{data:}}{ Names of the miRNAs whose interactions with the mRNA have
#'    been experimentally confirmed.}
#'  }
#' @export
#' @seealso \link{PTC}, \link{GroundT}
#' @return A \code{list} where each element has the following structure:
#'    \item{name}{ Name of a mRNA inferred by PTC with at least one confirmed interaction.}
#'    \item{data}{ Names of the miRNAs whose miRNA-mRNA interactions are in the \code{GroundT}.}
#' @examples \dontrun{
#'data(TCGA_BRCAdata)
#'data(GroundT)
#'test1<-PTC(miRNAs=TCGA_BRCAdata$miRs,mRNAs=TCGA_BRCAdata$mRNAs, VIM=TCGA_BRCAdata$mRNAs[,"VIM"])
#'t1.Confirmed<-Find.Confirmed(test1$Plist$Names,GroundT)
#' }

Find.Confirmed<-function(miRNAsList,GroundT){
  mRs<-names(miRNAsList)
  Confirmed<-vector("list")
  for (i in mRs) {
    for (j in miRNAsList[[i]]) {
      if(j%in%GroundT[[i]]){
        ifelse(is.null(Confirmed[[i]])
               ,Confirmed[[i]]<-j
               ,Confirmed[[i]]<-c(Confirmed[[i]],j)
        )
      }
    }#end for miR names
  }#end for genes
  return(Confirmed)
}
