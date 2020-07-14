#'@title wandTime.csv
#'
#'@description File containing the pseudotime obtained by using Wanderlust (Bendall et al) on
#' the dataset \link{SC_mRNAsdata}. Parameters were set to k=4, l=2, ng=2 and snn =0.
#'@name wandTime.csv
#'@docType data
#'@usage wandTime <- read.csv("./data/wandTime.csv", header=FALSE)
#'@format A \code{matrix} with 1 row and 19 columns with the pseudotime
#'obtained by Wanderlust using the dataset \link{SC_mRNAsdata}.
#'@references A Pseudo-Temporal Causality Approach to Identifying miRNA-mRNA
#' Interactions During Biological Processes\cr
#' Andres M. Cifuentes-Bernal, Vu VH Pham, Xiaomei Li,
#' Lin Liu, Jiuyong Li, Thuc Duy Le \cr
#' bioRxiv 2020.07.07.192724;
#'  \url{https://doi.org/10.1101/2020.07.07.192724}
#'
#'Wang, N., Zheng, J., Chen, Z. et al.
#'"Single-cell microRNA-mRNA co-sequencing reveals non-genetic heterogeneity
#'and mechanisms of microRNA regulation."
#'Nat Commun 10, 95 (2019).
#'https://doi.org/10.1038/s41467-018-07981-6
#'
#'Bendall SC, Davis KL, Amir el-AD, et al. \cr
#'Single-cell trajectory detection uncovers progression and
#'regulatory coordination in human B cell development. \cr
#'Cell. 2014;157(3):714-725. doi:10.1016/j.cell.2014.04.005
#'@keywords pseudotime Single_Cell
NULL
