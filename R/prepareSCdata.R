# Matched single-cell miRNA-mRNA data from GSE114071 (Wang et al., 2019),
# prepared as in inst/PTC.MainApplication.R. Used by the interface (runPTCApp)
# and by the tests. Returns a list with the miRs and mRNAs matrices
# (rows are cells, columns are genes).
prepareSCdata<-function(){
  data("SC_miRNAsdata", "SC_mRNAsdata", envir = environment())

  SC_miRNAsdata<-na.omit(SC_miRNAsdata[,-2])
  ## Reverse Log2 in miRNAs data
  SC_miRNAsdata[,-1]<-2^SC_miRNAsdata[,-1]
  ## Amplify miRNAs by the mean of total reads that
  ## map to miRNA in GSE114071
  SC_miRNAsdata[,-1]<-30.18*SC_miRNAsdata[,-1]

  SC_mRNAsdata<-na.omit(SC_mRNAsdata)
  mR.names<-SC_mRNAsdata[-1,1]
  mR.sample<-SC_mRNAsdata[1,-c(1,2)]
  SC_mRNAsdata<-mapply(SC_mRNAsdata[-1,-c(1,2),drop =FALSE], FUN=as.numeric, drop=FALSE)
  rownames(SC_mRNAsdata)<-mR.names
  colnames(SC_mRNAsdata)<-mR.sample
  ## Remove half-cells without matched miRNAs-mRNAs and transpose
  SC_mRNAsdata<-t(SC_mRNAsdata[,-c(6,21,22),drop=FALSE])
  ## Retain genes with mean expression greater than 0.1 and expressed in 20% of cells
  SC_mRNAsdata<-SC_mRNAsdata[,colMeans(SC_mRNAsdata) > 0.1 & colMeans(SC_mRNAsdata > 0) > 0.2]

  SC_miRNAsdata<-t(SC_miRNAsdata[,-c(21,22),drop=FALSE])
  ## Remove all duplicate miRs names by keeping only the 1st one
  SC.miR.Names<-unique(SC_miRNAsdata[1,])
  aux<-vector(mode="numeric",length = length(SC.miR.Names))
  names(aux)<-SC.miR.Names
  for (i in SC.miR.Names) {
    temp<-which(SC_miRNAsdata[1,]==i)
    aux[i]<-temp[1]
  }
  SC_miRNAsdata<-data.matrix(SC_miRNAsdata[,aux,drop=FALSE])
  colnames(SC_miRNAsdata)<-SC_miRNAsdata[1,]
  SC_miRNAsdata<-apply(SC_miRNAsdata[-1,],c(1,2),as.numeric)

  return(list(miRs = SC_miRNAsdata, mRNAs = SC_mRNAsdata))
}
