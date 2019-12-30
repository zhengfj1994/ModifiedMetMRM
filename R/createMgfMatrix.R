#' @name createMgfMatrix
#' @description Use mgf file to create a matrix
#' @author Fujian Zheng <zhengfj@dicp.ac.cn>
#'
#' @param mgfFile
#'
#' @return mgfMatrix mgfData
#' @example mgfMatrix <- createMgfMatrix(mgfData)

createMgfMatrix <- function(mgfData){
  require(stringr)
  # mgfData <- scan(mgfFile, what = character(0), sep = "\n")
  beginNum <- grep("BEGIN IONS", mgfData)
  pepmass <- grep("PEPMASS=",mgfData)
  trInSecond <- grep("RTINSECONDS=",mgfData)
  endNum <- grep("END IONS", mgfData)
  mgfMatrix <- cbind(beginNum,trInSecond,pepmass,endNum)

  for (i in c(1:length(pepmass))){
    pepmassi <- gsub("[^0-9,.]", "", strsplit(mgfData[pepmass[i]],split = " ")[[1]][1])
    mgfMatrix[i,"pepmass"] <- pepmassi

    tri <- gsub("[^0-9,.]", "", mgfData[trInSecond[i]])
    mgfMatrix[i,"trInSecond"] <- tri
  }
  return(mgfMatrix)
}
##########
