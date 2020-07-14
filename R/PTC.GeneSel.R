#' @title PTC.GeneSel
#'
#' @description Creates a set of miRNAs (plausible parents) for each mRNA among \code{nmR} mRNAs that can biologically target that mRNA.
#' A total of \code{nmiR} miRNAs and \code{nmR} mRNAs with the largest gene expression Median Absolute Deviation (MAD)
#' are selected as predictor candidates and target mRNAs respectively.
#' \cr The set of plausible parents of each gene is a subset of the predictor candidates and contains those miRNA that can bind
#' the gene as predicted by TargetScan 7.0.
#' miRNAs names are converted to \code{miRBase v.21} during selection process
#'
#' @inheritParams PTC
#' @param seqData A list with two elements:
#'    \itemize{
#'       \item A matrix with miRNAs gene expression.
#'       \item A matrix with mRNAs  gene expression.
#'    }
#' Columns represent miRNAs/mRNAs and rows represent samples (time points). Columns names must
#' correspond to miRNAs/mRNAs names.
#' @export
#' @seealso \link{PTC}, \link{getDatabyMAD},
#' \link[miRBaseConverter]{miRNAVersionConvert}, \link{PTC.findPP}
#' @return A list containing four elements:
#'    \item{\code{miRs}}{ \code{nmiR} miRNAs names (version \code{miRBase v.21})
#'    selected by MAD}
#'    \item{\code{mRs}}{ \code{nmR} mRNAs names selected by MAD}
#'    \item{\code{d}}{ A matrix with \code{nmiR + nmR} columns,
#'    containing miRNAs:mRNAs gene expression. Columns [1 : \code{nmiR}] correspond to
#'    miRNAs data. Columns [\code{nmiR}+1 : \code{nmiR}+\code{nmR}]} correspond to
#'    mRNAs data. Rows correspond to samples.
#'    \item{\code{PParents}}{ A list containing the set of miRNAs to be used as predictors
#'    (plausible parents) of each mRNA. This set contains the miRNAs among the predictor
#'    candidates that can bind the target gene as predicted by TargetScan 7.0}
#' @examples \dontrun{
#' seqData<-PTC.ptime(TCGA_BRCAdata,TCGA_BRCAdata$mRNAs[,"VIM"])
#' SelData<-PTC.GeneSel(seqData, nmiR = 30, nmR = 1500)
#' }

PTC.GeneSel<-function(seqData, nmiR = 30, nmR = 1500){
  data(TScan, envir = environment())

  l <- getDatabyMAD(seqData, nmiR = 30, nmR = 1500)

##### ------Change selected miRNAs names to v21
  l$miRs<-PTC.miRv21(colnames(l$d[,1:nmiR]))

#######-----Find miRNAs that can bind each gene (TargetScan7.0)
  PParents<-PTC.findPP(TScan, miRs=l$miRs, mRs=l$mRs)
  l$PParents<-PParents
  return(l)
  }
