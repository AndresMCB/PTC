#' @title PTC
#'
#' @description PTC estimates the causal parents of a set of \code{nmR} mRNAs,
#' given a set of \code{nmiR} predictors (miRNAs). This estimation assumes
#' a linear model \deqn{Y_t = \beta X_t + \epsilon}
#' With \eqn{Y_t, X_t} the sequential
#' data of the target and the predictors respectively, and \eqn{\epsilon} independent and identically
#' distributed errors for the n time points in the time series.
#' \eqn{Y_t, X_t} are the sequential data obtained after a pseudotime analysis.
#' @usage PTC(miRNAs, mRNAs, VIM, nmiR=30, nmR=1500 , ngrid=2, alpha=0.02, complements = TRUE
#' , explore.all=TRUE, silent=TRUE)\cr
#'
#' @param miRNAs A \code{matrix} containing miRNA gene expression.\cr
#' A total of \code{nmiR} miRNAs from this matrix are selected to be used
#' as predictors candidates (i.e. plausible parents).
#' Columns represent miRNAs, rows represent samples.
#' @param mRNAs A \code{matrix} containing mRNA gene expression.\cr
#' A total of \code{nmR} mRNAs from this matrix are
#' selected to be used as response variables.
#' Columns represent mRNAs, rows represent samples.
#' @param VIM VIM expression to be used for calculating VIM_Time.
#' @param nmiR Number of miRNAs to be selected as predictor candidates.
#' @param nmR Number of mRNAs to be selected as target variables.
#' @param ngrid Number of segments of the time series data used for creating
#' the different enviroments required for the statistical test.
#' \code{ngrid=2} by default.
#' @param alpha Significance level for the statistical test.
#' \code{alpha=0.02} by default.
#' @param complements If \code{TRUE} (default), each environment is compared against
#' its complement. If \code{FALSE} all environments are compared pairwise.
#' @param explore.all If \code{TRUE}(default), PTC explores all combinations of predictors
#' and returns the union set of all combinations that does not violate the invariance
#' property. If \code{FALSE} PTC returns the first set that does not violate the invariance
#' property. Exploration is made starting from the set of all predictors and reducing
#' the size of the set by one predictor at a time.
#' @param silent If \code{TRUE} (default), PTC displays the currently evaluated set.
#' If \code{FALSE}, PTC only displays the number of sets to be explored in
#' the current iteration.
#' @param TScan A 2 column matrix containing the miRNA-mRNA relationships.
#' for an appropriate functioning, miRNAs must be in lower case and mRNAs in Uppercase.
#' If \code{NULL} (default), PTC uses a pre-loaded matrix (file TScan) with TargetScan 7.0 Human.
#'
#' @author Andres Mauricio Cifuentes_Bernal, Vu VH Pham, Xiaomei Li, Lin Liu, JiuyongLi and Thuc Duy Le
#' @export
#' @seealso \link[PTC]{PTC.ptime}, \link[PTC]{PTC.GeneSel}, \link{TScan}
#' \link[PTC]{PTC.TestInvariance}.
#'
#' @return A \code{list} consisting of the following elements:
#'   \item{\code{Index}}{A \code{list} where each element is a target mRNA.
#'    For each target gene with at least one parent. The index of the parents.}
#'   \item{\code{names}}{A \code{list} where each element is a target mRNA.
#'    For each target gene with at least one parent. The name of the parents.}
#'   \item{\code{genes}}{The names of all target genes with at least one parent.}
#'   \item{\code{Summary}}{A matrix representing miRNA-mRNA regulatory interactions
#'   inferred by PTC. The columns of the matrix are:
#'        \enumerate{
#'           \item{\code{rank:}}{ Rank of the inferred interaction}
#'           \item{\code{miR:}}{ Names of miRNA (Parent)}
#'           \item{\code{mR:}}{ Names of mRNA (Child)}
#'           \item{\code{Score:}}{ Score as calculated by \link{PTC.RankByContext}}
#'           }
#'    }
#' @examples \dontrun{
#'    data(TCGA_BRCAdata)
#'    test1<-PTC(miRNAs=TCGA_BRCAdata$miRs, mRNAs=TCGA_BRCAdata$mRNAs
#'           , VIM=TCGA_BRCAdata$mRNAs[,"VIM"])
#' }
#'
#' @references
#' A Pseudo-Temporal Causality Approach to Identifying miRNA-mRNA
#' Interactions During Biological Processes\cr
#' Andres M. Cifuentes-Bernal, Vu VH Pham, Xiaomei Li,
#' Lin Liu, Jiuyong Li, Thuc Duy Le \cr
#' bioRxiv 2020.07.07.192724;
#'  \url{https://doi.org/10.1101/2020.07.07.192724}
#'

PTC<-function(miRNAs, mRNAs, VIM, nmiR=30, nmR=1500
              , ngrid=2, alpha=0.02, complements = TRUE
              , explore.all=TRUE, silent=TRUE, TScan=NULL)
  {
  #for compatibility, mRNAs are trasformed to uppercase
  colnames(mRNAs)<-toupper(colnames(mRNAs))

  data("TS7.0_Conserved_Site_Context_Scores",envir = environment())

  GEData<-list(miRNAs,mRNAs)
  names(GEData)<-c("miRs","mRNAs")
  seqData<-PTC.ptime(GEData,VIM)
  SelData<-PTC.GeneSel(seqData, nmiR = nmiR, nmR = nmR, TScan=TScan)

  PTC.outcome<-vector("list",length = length(SelData$mRs))
  names(PTC.outcome)<-SelData$mRs

  for (gene in SelData$mRs) {
    temp <-SelData$PParents[[gene]]
    X.TScan=SelData$d[,temp,drop = FALSE]
    if(length(X.TScan)>0){
      message(paste(" Current mRNA  = {",gene,"}",sep=" "))
      set.seed(1)
      PTC.outcome[gene]<-
        try(PTC.TestInvariance(Y = SelData$d[,gene,drop = FALSE]
                              ,X = X.TScan
                              ,ngrid = ngrid
                              ,alpha = alpha
                              ,explore.all = explore.all
                              ,silent=silent
                              ,complements=complements))
    }
  }
  Results<-list()

  Results<-Extract.Parents(PTC.outcome,SelData$PParents)
  Results$Summary<-InterList.toMatrix(Results$Names)
  Results$Summary<-PTC.RankByContext(TS7.0_Conserved_Site_Context_Scores
                                    ,Results$Summary)

  return(Results)
}
