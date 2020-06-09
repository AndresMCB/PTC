library(PTC)

######################Test with bulk data ##############
data(TCGA_BRCAdata)#optional, TCGA_BRCAdata is lazzy data (preloaded)
data(GroundT)

test1<-PTC(miRNAs=TCGA_BRCAdata$miRs,
          mRNAs=TCGA_BRCAdata$mRNAs, VIM=TCGA_BRCAdata$mRNAs[,"VIM"])
t1.Confirmed<-Find.Confirmed(test1$Plist$Names,GroundT)
t1.mConfirmed<-InterList.toMatrix(t1.Confirmed)


# ----- Formating SCdata to create matched data
SC_miRNAsdata<-na.omit(SC_miRNAsdata[,-2])
## Reverse Log2 in miRNAs data
SC_miRNAsdata[,-1]<-2^SC_miRNAsdata[,-1]
## Amplify miRNAs by the mean of total reads that
## map to miRNA in GSE114071
SC_miRNAsdata[,-1]<-30.18*SC_miRNAsdata[,-1]

SC_mRNAsdata<-na.omit(SC_mRNAsdata)
colnames(SC_mRNAsdata)<-SC_mRNAsdata[1,]
rownames(SC_mRNAsdata)<-SC_mRNAsdata[,1]
SC_mRNAsdata<-data.matrix(SC_mRNAsdata[-1,-c(1,2)])

## Remove half-cells without matched miRNAs-mRNAs
## and transpose
SC_mRNAsdata<-t(SC_mRNAsdata[,-c(6,21,22),drop=FALSE])
#suggested retain only genes with mean expression greater than 1 and expressed in 20% of cells
SC_mRNAsdata<- SC_mRNAsdata[,colMeans(SC_mRNAsdata) > 0.1 & colMeans(SC_mRNAsdata > 0) > 0.2]


SC_miRNAsdata<-t(SC_miRNAsdata[,-c(21,22),drop=FALSE])

# Remove all duplicate miRs names by letting only the 1st one
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

##         PTC Single Cell data
test2<-PTC(miRNAs=SC_miRNAsdata, mRNAs=SC_mRNAsdata
           , VIM=SC_mRNAsdata[,"VIM"])
