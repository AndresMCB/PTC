#' @title PTC.miRv21
#' @description Takes a vector containing miRNAs names and returns
#' the \code{miRBase v.21} version of them.
#'
#' @importFrom miRBaseConverter miRNAVersionConvert
#' @param miRs miRNAs names to be changed to version \code{miRBase v.21}
#' @export
#'
#' @seealso \link{PTC}, \link[miRBaseConverter]{miRNAVersionConvert},
#' \link{PTC.GeneSel}
#' @return \code{nmiR} A vector with the miRNAs names
#' (version \code{miRBase v.21})
#'
#' @examples \dontrun{
#' data(TCGA_BRCAdata)
#' miRnamesv21<-PTC.miRv21(colnames(TCGA_BRCAdata$miRs))
#'  }

PTC.miRv21  <- function(miRs) {
  miRs<-gsub("\\..*","",miRs)
  aux<-miRNAVersionConvert(miRs
                           , targetVersion = "v21", exact = T
                           , verbose = T)
  #check if any conversion went wrong
  NA.index<-which(is.na(aux[,2]))
  #recover those miRs names
  if(length(NA.index)>0){
    message(paste(length(NA.index)," conversions failed, original name is kept"))
    aux[NA.index,2]<-miRs[NA.index]
  }
  #if there are multiple matched miRNAs, we choose the first one
  aux[,2]<-gsub("\\&.*","",aux[,2])
  miRs<-aux[,2]
  rm(aux,NA.index)
  return(miRs)
}
