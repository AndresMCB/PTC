% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PTC.ptime.R
\name{PTC.ptime}
\alias{PTC.ptime}
\title{PTC.ptime}
\usage{
PTC.ptime(matchedData, VIM)
}
\arguments{
\item{matchedData}{A \code{List} with two elements. A miRNA gene expressions matrix (1st element)
and A mRNA gene expressions matrix (2nd element).
Columns represent miRNAs and rows represent samples.}

\item{VIM}{VIM expression to be used for calculating VIM_Time.}
}
\value{
Pseudotime ordered matched data.
}
\description{
Orders the gene expression by VIM_Time order.
}
\examples{
\dontrun{
  data(TCGA_BRCAdata)
  Time_series<-PTC.ptime(TCGA_BRCAdata, TCGA_BRCAdata$mRNAs[,"VIM"])
}
}
\seealso{
\link[PTC]{PTC}
}
