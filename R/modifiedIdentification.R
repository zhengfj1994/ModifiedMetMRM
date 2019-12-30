#' @name modifiedIdentification
#' @description Use NL to identify modified metabolites
#' @author Fujian Zheng <zhengfj@dicp.ac.cn>
#' @param allIonPairsMatrix
#' @param modificationList
#' @param tolMzppm
#' @return modifiedIdentificationRes
#' @example modifiedIdentificationRes <- modifiedIdentification(allIonPairsMatrix, modificationList="D:\\github\\ModifiedMetMRM\\示例文件\\Modification List.xlsx", tolMZppm=20)

modifiedIdentification <- function(allIonPairsMatrix, modificationList, tolMZppm){
  require(xlsx)
  require(tcltk)
  modificationData <- read.xlsx2(modificationList,sheetIndex = 1)
  pb <- tkProgressBar("modifiedIdentification", "Rate of progress %", 0, 100)
  for (i in c(1:nrow(allIonPairsMatrix))){
    info<- sprintf("Rate of progress %d%%", round(i*100/nrow(allIonPairsMatrix)))
    setTkProgressBar(pb, i*100/nrow(allIonPairsMatrix), sprintf("modifiedIdentification", info),info)

    NLi <- allIonPairsMatrix$neutralLoss[i]
    posi <- which(abs(as.numeric(as.character(modificationData$m.z)) - NLi)/as.numeric(as.character(modificationData$m.z))*1000000 < tolMZppm)
    if (length(posi)==1){
      allIonPairsMatrix[i,"modifiedType"] <- modificationData$Modified.Type[posi]
    }
    else if (length(posi)>1){
      modifiedTypeName <- as.character(modificationData$Modified.Type[posi[1]])
      for (j in c(2:length(posi))){
        modifiedTypeName <- paste0(modifiedTypeName,";",as.character(modificationData$Modified.Type[posi[j]]))
      }
      allIonPairsMatrix[i,"modifiedType"] <- modifiedTypeName
    }
  }
  close(pb)
  modifiedIdentificationRes <- na.omit(allIonPairsMatrix)
  return(modifiedIdentificationRes)
}
