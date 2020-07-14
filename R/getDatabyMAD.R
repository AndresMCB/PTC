#' @title getDatabyMAD
#' @inheritParams PTC.GeneSel
#' @importFrom CancerSubtypes FSbyMAD
#' @author Vu VH Pham <vu_viet_hoang.pham@@mymail.unisa.edu.au>,
#'         Taosheng Xu <taosheng.x@@gmail.com>
#' @description Given a gene expression data including matched miRNAs-mRNAs,
#'    finds the \code{nmiR} miRNAs and the \code{nmR} mRNAs with the largest
#'    Median Absolute Deviation (MAD).
#' @return A \code{list} containing samples of miRs and mRs. The elements of the \code{list} are:
#'    \item{\code{d}}{ The data with rows being samples and columns being miRs and mRs.}
#'    \item{\code{miRs}}{ The names of the selected miRs.}
#'    \item{\code{mRs}}{ The names of the selected mRs.}
#' @export
#' @seealso \link[CancerSubtypes]{FSbyMAD}
#' @examples
#' data(TCGA_BRCAdata)
#' nmiR <- 30
#' nmR <- 1500
#' l <- getDatabyMAD(TCGA_BRCAdata, nmiR, nmR)
#' @references Taosheng Xu, Thuc Duy Le, Lin Liu, Ning Su,
#'  Rujing Wang, Bingyu Sun, Antonio Colaprico,
#'  Gianluca Bontempi, Jiuyong Li.
#'  CancerSubtypes: an R/Bioconductor package for molecular
#'  cancer subtype identification, validation, and visualization. Bioinformatics 33(19): 3131–3133 (2017).
#'  \url{https://doi.org/10.1093/bioinformatics/btx378}

getDatabyMAD <- function(seqData, nmiR, nmR) {
  if (!require(CancerSubtypes)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("CancerSubtypes")
  }
  library(CancerSubtypes)

  data(TScan, envir = environment())

  # Identify significant miRNAs by using function FSbyMAD in CancerSubtypes package
  # miRNAsData is a matriz with nmiR rows
  miRNAsData = FSbyMAD(t(seqData$miRs), value = nmiR)

  # miRNAsData is transposed
  miRNAsData <- t(miRNAsData)
  #Get column names (sample names)
  miRs <- colnames(miRNAsData)

  # Identify significant mRNAs by using function FSbyMAD in CancerSubtypes package
  mRNAsData = FSbyMAD(t(seqData$mRNAs), value=nmR)
  mRNAsData <- t(mRNAsData)
  # Remove duplicated data
  mRNAsData <- mRNAsData[,!duplicated(colnames(mRNAsData))]

  mRs <- colnames(mRNAsData)

  # Combine data
  d <- cbind(miRNAsData, mRNAsData)
  l = list(d = d, miRs = miRs, mRs = mRs)
  return(l)
}
