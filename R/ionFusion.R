#' @name ionFusion
#' @description Find isotopes, adduction and others
#' @author Fujian Zheng <zhengfj@dicp.ac.cn>
#' @param rawMatrix
#' @param tolMZppm
#' @param tolTr
#' @param ionMode
#' @return afterIonFusionRes
#' @example afterIonFusionRes <- ionFusion(rawMatrix=modifiedIdentificationRes, tolMZppm=10, tolTr=0.1, ionMode="positive")

ionFusion <- function(rawMatrix, tolMZppm, tolTr, ionMode){
  require(tcltk)
  rawMatrix$IsoAdductLable <- NA
  if (ionMode=="positive"){
    pb <- tkProgressBar("ionFusion", "Rate of progress %", 0, 100)
    for (i in c(1:(nrow(rawMatrix)-1))){
      info<- sprintf("Rate of progress %d%%", round(i*100/(nrow(rawMatrix)-1)))
      setTkProgressBar(pb, i*100/(nrow(rawMatrix)-1), sprintf("ionFusion", info),info)

      mzi <- rawMatrix$m.z[i]
      tri <- rawMatrix$RT[i]
      Indexi <- rawMatrix$Index[i]
      for (j in c(i:nrow(rawMatrix))){
        mzj <- rawMatrix$m.z[j]
        trj <- rawMatrix$RT[j]
        Indexj <- rawMatrix$Index[j]
        if (abs(mzj-mzi-1.007825)/mzi*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[j,"IsoAdductLable"] <- paste("[M+1] of",Indexi)
          }
          else{
            rawMatrix[j,"IsoAdductLable"] <- paste0(rawMatrix[j,"IsoAdductLable"],";",paste("[M+1] of",Indexi))
          }
        }
        if (abs(mzi-mzj-1.007825)/mzj*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[i,"IsoAdductLable"] <- paste("[M+1] of",Indexj)
          }
          else{
            rawMatrix[i,"IsoAdductLable"] <- paste0(rawMatrix[i,"IsoAdductLable"],";",paste("[M+1] of",Indexj))
          }
        }
        if (abs(mzj-mzi-21.981942)/mzi*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[j,"IsoAdductLable"] <- paste("[M+Na] of",Indexi)
          }
          else{
            rawMatrix[j,"IsoAdductLable"] <- paste0(rawMatrix[j,"IsoAdductLable"],";",paste("[M+Na] of",Indexi))
          }
        }
        if (abs(mzi-mzj-21.981942)/mzj*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[i,"IsoAdductLable"] <- paste("[M+Na] of",Indexj)
          }
          else{
            rawMatrix[i,"IsoAdductLable"] <- paste0(rawMatrix[i,"IsoAdductLable"],";",paste("[M+Na] of",Indexj))
          }
        }
        if (abs(mzj-mzi-37.955882)/mzi*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[j,"IsoAdductLable"] <- paste("[M+K] of",Indexi)
          }
          else{
            rawMatrix[j,"IsoAdductLable"] <- paste0(rawMatrix[j,"IsoAdductLable"],";",paste("[M+K] of",Indexi))
          }
        }
        if (abs(mzi-mzj-37.955882)/mzj*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[i,"IsoAdductLable"] <- paste("[M+K] of",Indexj)
          }
          else{
            rawMatrix[i,"IsoAdductLable"] <- paste0(rawMatrix[i,"IsoAdductLable"],";",paste("[M+K] of",Indexj))
          }
        }
        if (abs(mzj-mzi-17.026547)/mzi*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[j,"IsoAdductLable"] <- paste("[M+NH4] of",Indexi)
          }
          else{
            rawMatrix[j,"IsoAdductLable"] <- paste0(rawMatrix[j,"IsoAdductLable"],";",paste("[M+NH4] of",Indexi))
          }
        }
        if (abs(mzi-mzj-17.026547)/mzj*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[i,"IsoAdductLable"] <- paste("[M+NH4] of",Indexj)
          }
          else{
            rawMatrix[i,"IsoAdductLable"] <- paste0(rawMatrix[i,"IsoAdductLable"],";",paste("[M+NH4] of",Indexj))
          }
        }
        if (abs(mzi-mzj-18.01056)/mzi*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[j,"IsoAdductLable"] <- paste("[M+H-H2O] of",Indexi)
          }
          else{
            rawMatrix[j,"IsoAdductLable"] <- paste0(rawMatrix[j,"IsoAdductLable"],";",paste("[M+H-H2O] of",Indexi))
          }
        }
        if (abs(mzj-mzi-18.01056)/mzj*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[i,"IsoAdductLable"] <- paste("[M+H-H2O] of",Indexj)
          }
          else{
            rawMatrix[i,"IsoAdductLable"] <- paste0(rawMatrix[i,"IsoAdductLable"],";",paste("[M+H-H2O] of",Indexj))
          }
        }
        if (abs((mzj-1.007276)/2+1.007276-mzi)/mzi*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[j,"IsoAdductLable"] <- paste("[2M+H] of",Indexi)
          }
          else{
            rawMatrix[j,"IsoAdductLable"] <- paste0(rawMatrix[j,"IsoAdductLable"],";",paste("[2M+H] of",Indexi))
          }
        }
        if (abs((mzi-1.007276)/2+1.007276-mzj)/mzj*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[i,"IsoAdductLable"] <- paste("[2M+H] of",Indexj)
          }
          else{
            rawMatrix[i,"IsoAdductLable"] <- paste0(rawMatrix[i,"IsoAdductLable"],";",paste("[2M+H] of",Indexj))
          }
        }
      }
    }
    close(pb)
  }
  else if (ionMode=="negative"){
    pb <- tkProgressBar("ionFusion", "Rate of progress %", 0, 100)
    for (i in c(1:(nrow(rawMatrix)-1))){
      info<- sprintf("Rate of progress %d%%", round(i*100/(nrow(rawMatrix)-1)))
      setTkProgressBar(pb, i*100/(nrow(rawMatrix)-1), sprintf("ionFusion", info),info)

      mzi <- rawMatrix$m.z[i]
      tri <- rawMatrix$RT[i]
      Indexi <- rawMatrix$Index[i]
      for (j in c(i:nrow(rawMatrix))){
        mzj <- rawMatrix$m.z[j]
        trj <- rawMatrix$RT[j]
        Indexj <- rawMatrix$Index[j]
        if (abs(mzj-mzi-1.007825)/mzi*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[j,"IsoAdductLable"] <- paste("[M+1] of",Indexi)
          }
          else{
            rawMatrix[j,"IsoAdductLable"] <- paste0(rawMatrix[j,"IsoAdductLable"],";",paste("[M+1] of",Indexi))
          }
        }
        if (abs(mzi-mzj-1.007825)/mzj*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[i,"IsoAdductLable"] <- paste("[M+1] of",Indexj)
          }
          else{
            rawMatrix[i,"IsoAdductLable"] <- paste0(rawMatrix[i,"IsoAdductLable"],";",paste("[M+1] of",Indexj))
          }
        }
        if (abs(mzj-mzi-35.976678)/mzi*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[j,"IsoAdductLable"] <- paste("[M+Cl] of",Indexi)
          }
          else{
            rawMatrix[j,"IsoAdductLable"] <- paste0(rawMatrix[j,"IsoAdductLable"],";",paste("[M+Cl] of",Indexi))
          }
        }
        if (abs(mzi-mzj-35.976678)/mzj*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[i,"IsoAdductLable"] <- paste("[M+Cl] of",Indexj)
          }
          else{
            rawMatrix[i,"IsoAdductLable"] <- paste0(rawMatrix[i,"IsoAdductLable"],";",paste("[M+Cl] of",Indexj))
          }
        }
        if (abs(mzi-mzj-18.01056)/mzi*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[j,"IsoAdductLable"] <- paste("[M-H-H2O] of",Indexi)
          }
          else{
            rawMatrix[j,"IsoAdductLable"] <- paste0(rawMatrix[j,"IsoAdductLable"],";",paste("[M-H-H2O] of",Indexi))
          }
        }
        if (abs(mzj-mzi-18.01056)/mzj*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[i,"IsoAdductLable"] <- paste("[M-H-H2O] of",Indexj)
          }
          else{
            rawMatrix[i,"IsoAdductLable"] <- paste0(rawMatrix[i,"IsoAdductLable"],";",paste("[M-H-H2O] of",Indexj))
          }
        }
        if ((abs(mzj+1.007276)/2-1.007276-mzi)/mzi*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[j,"IsoAdductLable"] <- paste("[2M-H] of",Indexi)
          }
          else{
            rawMatrix[j,"IsoAdductLable"] <- paste0(rawMatrix[j,"IsoAdductLable"],";",paste("[2M-H] of",Indexi))
          }
        }
        if ((abs(mzi+1.007276)/2-1.007276-mzj)/mzj*1000000 < tolMZppm & abs(tri-trj) < tolTr){
          if (is.na(rawMatrix[j,"IsoAdductLable"])){
            rawMatrix[i,"IsoAdductLable"] <- paste("[2M-H] of",Indexj)
          }
          else{
            rawMatrix[i,"IsoAdductLable"] <- paste0(rawMatrix[i,"IsoAdductLable"],";",paste("[2M-H] of",Indexj))
          }
        }
      }
    }
    close(pb)
  }
  afterIonFusionRes <- rawMatrix[which(is.na(rawMatrix$IsoAdductLable)),]
  afterIonFusionRes <- afterIonFusionRes[,-ncol(afterIonFusionRes)]
  return(afterIonFusionRes)
}

