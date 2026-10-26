#' @title getDatabyMAD
#' @inheritParams PTC.GeneSel
#' @importFrom stats mad
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
#' @details The selection reproduces \code{FSbyMAD(cut.type = "topk")} from the
#'    CancerSubtypes package (now deprecated in Bioconductor), without drawing its histogram.
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
  # Identify significant miRNAs by MAD (as CancerSubtypes::FSbyMAD)
  # miRNAsData is a matrix with nmiR rows
  miRNAsData = selectByMAD(t(seqData$miRs), value = nmiR)

  # miRNAsData is transposed
  miRNAsData <- t(miRNAsData)
  # Get column names (miRNA names)
  miRs <- colnames(miRNAsData)

  # Identify significant mRNAs by MAD (as CancerSubtypes::FSbyMAD)
  mRNAsData = selectByMAD(t(seqData$mRNAs), value=nmR)
  mRNAsData <- t(mRNAsData)
  # Remove duplicated data
  mRNAsData <- mRNAsData[,!duplicated(colnames(mRNAsData))]

  mRs <- colnames(mRNAsData)

  # Combine data
  d <- cbind(miRNAsData, mRNAsData)
  l = list(d = d, miRs = miRs, mRs = mRs)
  return(l)
}

# Keeps the `value` rows of `Data` with the largest MAD. Same selection as
# CancerSubtypes::FSbyMAD(Data, cut.type = "topk", value) (Xu et al., 2017;
# GPL >= 2), which also drew a histogram of the MADs as a side effect.
selectByMAD <- function(Data, value) {
  mads = apply(Data, 1, mad)
  index = sort(mads, decreasing = TRUE, index.return = TRUE)
  if (value > nrow(Data)) {
    value = nrow(Data)
    cat("Warning: the feature selection number is beyond the original feature numnber")
  }
  index = index$ix[1:value]
  Data[index, ]
}
