# Matched single-cell data, prepared as in inst/PTC.MainApplication.R
prepare_SCdata <- function() {
  SC_miRNAsdata <- PTC::SC_miRNAsdata
  SC_mRNAsdata <- PTC::SC_mRNAsdata

  SC_miRNAsdata<-na.omit(SC_miRNAsdata[,-2])
  SC_miRNAsdata[,-1]<-2^SC_miRNAsdata[,-1]
  SC_miRNAsdata[,-1]<-30.18*SC_miRNAsdata[,-1]

  SC_mRNAsdata<-na.omit(SC_mRNAsdata)
  mR.names<-SC_mRNAsdata[-1,1]
  mR.sample<-SC_mRNAsdata[1,-c(1,2)]
  SC_mRNAsdata<-mapply(SC_mRNAsdata[-1,-c(1,2),drop =F], FUN=as.numeric, drop=F)
  rownames(SC_mRNAsdata)<-mR.names
  colnames(SC_mRNAsdata)<-mR.sample
  SC_mRNAsdata<-t(SC_mRNAsdata[,-c(6,21,22),drop=FALSE])
  SC_mRNAsdata<- SC_mRNAsdata[,colMeans(SC_mRNAsdata) > 0.1 & colMeans(SC_mRNAsdata > 0) > 0.2]

  SC_miRNAsdata<-t(SC_miRNAsdata[,-c(21,22),drop=FALSE])
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

  list(miRs = SC_miRNAsdata, mRNAs = SC_mRNAsdata)
}
