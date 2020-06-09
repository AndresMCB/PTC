#========================ConfirmedInTopK==============================
#' This fuction finds the number of confirmed interactions, given
#' a matrix of interacions miR-mR, and a Ground Truth
#' @param GT: a list where each element has
#'  name=gene, data= confirmed miRs able to bind that gene
#' @param relationshipsM: A matrix [n,2]
#' @param topK: a vector including all topK
#' @return A matrix [n,3], where the additional column has a score
#' calculated as the mean of all context score of each miR-mR pair
#========================================================================
PTC.ConfirmedInTopK <-function(GT,relationshipsM,
                           topK, GT_miR_EMT=NULL,GT_mR_EMT=NULL){

  ConfirmedInTop<-vector("list",length = length(topK)+1)
  names(ConfirmedInTop)<-c(paste("Top ",topK, sep=""),"summary")
  n<-3
  if(!is.null(GT_miR_EMT)){n=n+1}
  if(!is.null(GT_mR_EMT)){n=n+1}

  ConfirmedInTop[["summary"]]<-matrix(nrow=length(topK), ncol=n)

  for (i in 1:length(topK)) {
    if(topK[i]=="all"){topK[i]<-nrow(relationshipsM)}

    ConfirmedInTop[[i]][["IN_GT"]]<-
      Find.Confirmed2(relationshipsM[1:topK[i],1:2],GT)


    ConfirmedInTop[["summary"]][i,1]<-topK[i]
    ConfirmedInTop[["summary"]][i,2]<-nrow(na.omit(ConfirmedInTop[[i]][["IN_GT"]]))
    m<-nrow(na.omit(ConfirmedInTop[[i]][["IN_GT"]]))/as.numeric(topK[i])
    ConfirmedInTop[["summary"]][i,3]<-
      paste(round(100*m, 3), "%", sep="")
    aux<-0
    if(!is.null(GT_miR_EMT)){
      ConfirmedInTop[[i]][["miR.EMT"]]<-intersect(unique(relationshipsM[1:topK[i],1]),GT_miR_EMT[,1])
      ConfirmedInTop[[i]][["miR.EMT"]]<-data.matrix(ConfirmedInTop[[i]][["miR.EMT"]])
      colnames(ConfirmedInTop[[i]][["miR.EMT"]])<-"miR"

      aux<-aux+1
      ConfirmedInTop[["summary"]][i,3+aux]<-nrow(na.omit(ConfirmedInTop[[i]][["miR.EMT"]]))
    }


    if(!is.null(GT_mR_EMT)){
      ConfirmedInTop[[i]][["mR.EMT"]]<-intersect(unique(relationshipsM[1:topK[i],2]),GT_mR_EMT[,1])
      ConfirmedInTop[[i]][["mR.EMT"]]<-data.matrix(ConfirmedInTop[[i]][["mR.EMT"]])
      colnames(ConfirmedInTop[[i]][["mR.EMT"]])<-"mR"

      aux<-aux+1
      ConfirmedInTop[["summary"]][i,3+aux]<-nrow(na.omit(ConfirmedInTop[[i]][["mR.EMT"]]))
    }
  }
  colnames(ConfirmedInTop[["summary"]])<-rep("",ncol(ConfirmedInTop[["summary"]]))
  colnames(ConfirmedInTop[["summary"]])[1:3]<-c("Top","Confirmed","Percentage")
  aux<-0
  if(!is.null(GT_miR_EMT)){
    aux<-aux+1
    colnames(ConfirmedInTop[["summary"]])[3+aux]<-c("miR_EMT")
  }
  if(!is.null(GT_mR_EMT)){
    aux<-aux+1
    colnames(ConfirmedInTop[["summary"]])[3+aux]<-c("mR_EMT")
  }
  return(ConfirmedInTop)
}
