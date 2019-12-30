#' @name createAllIonPairs
#' @description Use MS1 peak table and mgf(MS2) to create ion pairs.
#' @author Fujian Zheng <zhengfj@dicp.ac.cn>
#' @param MS1file: MS1 peak table
#' @param MS2path: the folder saved mgf files
#' @param tolMZppm: the tolerence of m/z between ms1 and ms2
#' @param tolTr: the tolerence of retention time between ms1 and ms2
#' @return allIonPairsMatrix
#' @example allIonPairsMatrix <- createAllIonPairs(MS1file="D:\\github\\ModifiedMetMRM\\示例文件\\POS_一级峰表\\POS-TIC.txt",
#                                                  MS2path="D:\\github\\ModifiedMetMRM\\示例文件\\POS_二级文件",
#                                                  tolMZppm=20,
#                                                  tolTr=0.2)

createAllIonPairs <- function(MS1file, MS2path, tolMZppm, tolTr){
  require(tcltk)
  require(stringr)
  require(readr)
  require(dplyr)
  MS1data <- read.table(MS1file,header=TRUE)
  mgfFiles <- list.files(MS2path)

  mgfMatrixInfo <- matrix(c(1:4),nrow=1,ncol=4,dimnames=list(c("r1"),c("beginNum","trInSecond","pepmass","endNum")))
  productIonInfo <- matrix(c(1:4),nrow=1,ncol=4,dimnames=list(c("r1"),c("productIon","intensityOfProductIon","neutralLoss","CE")))
  allIonPairsMatrix <- cbind(MS1data[0,], mgfMatrixInfo[0,], productIonInfo[0,])
  Index = RT = m.z = beginNum = trInSecond = pepmass = endNum = productIon = intensityOfproductIon = neutralLoss = CE = NULL

  for (mgfi in mgfFiles){
    mgfData <- scan(paste0(MS2path,'\\',mgfi), what = character(0), sep = "\n")
    mgfMatrix <- as.data.frame(createMgfMatrix(mgfData))
    CEi <- parse_number(mgfi)

    pb <- tkProgressBar(paste0("createAllIonPairs for ",mgfi) ,"Rate of progress %", 0, 100)
    for(j in c(1:nrow(MS1data))){
      info<- sprintf("Rate of progress %d%%", round(j*100/nrow(MS1data)))
      setTkProgressBar(pb, j*100/nrow(MS1data), sprintf(paste0("createAllIonPairs for ",mgfi), info),info)

      trj <- MS1data$RT[j]
      mzj <- MS1data$m.z[j]
      posj <- which(abs(as.numeric(as.character(mgfMatrix$pepmass)) - mzj)/mzj*1000000 < tolMZppm & abs(as.numeric(as.character(mgfMatrix$trInSecond)) - trj*60) < tolTr*60)
      if (length(posj)>=1){
        posj <- posj[which.min(abs(as.numeric(mgfMatrix$trInSecond[posj])-trj*60))]
        beginNumj <- as.numeric(as.character(mgfMatrix$beginNum[posj]))
        endNumj <- as.numeric(as.character(mgfMatrix$endNum[posj]))

        if (str_detect(mgfData[beginNumj+2],'PEPMASS=') & str_detect(mgfData[beginNumj+3],'RTINSECONDS=')){
          addNum <- 4
        }
        else if (str_detect(mgfData[beginNumj+3],'PEPMASS=') & str_detect(mgfData[beginNumj+4],'RTINSECONDS=')){
          addNum <- 5
        }
        else{
          packageStartupMessage("The mgf file is not create!!!")
          stop()
        }
        if (str_detect(mgfData[beginNumj+addNum],'END IONS')){
          next()
        }
        mgfDataMS2 <- ms2dataframe(mgfData[(beginNumj+addNum):(endNumj-1)])
        addMgfMatrix <- matrix(rep(as.matrix(mgfMatrix[posj,]),nrow(mgfDataMS2)),nrow = nrow(mgfDataMS2),byrow=T)
        colnames(addMgfMatrix) <- c("beginNum","trInSecond","pepmass","endNum")

        addMS1data <- matrix(rep(as.matrix(MS1data[j,]),nrow(mgfDataMS2)),nrow = nrow(mgfDataMS2),byrow=T)
        colnames(addMS1data) <- colnames(MS1data)

        allIonPairsMatrixi <- cbind(addMS1data,addMgfMatrix,mgfDataMS2)
        allIonPairsMatrixi$neutralLoss <- as.numeric(as.character(allIonPairsMatrixi$pepmass))-allIonPairsMatrixi$productIon
        allIonPairsMatrixi$CE <- CEi
        Index = c(Index,allIonPairsMatrixi$Index)
        RT = c(RT,allIonPairsMatrixi$RT)
        m.z = c(m.z,allIonPairsMatrixi$m.z)
        beginNum = c(beginNum,as.numeric(as.character(allIonPairsMatrixi$beginNum)))
        trInSecond = c(trInSecond,as.numeric(as.character(allIonPairsMatrixi$trInSecond)))
        pepmass = c(pepmass,as.numeric(as.character(allIonPairsMatrixi$pepmass)))
        endNum = c(endNum,as.numeric(as.character(allIonPairsMatrixi$endNum)))
        productIon = c(productIon,allIonPairsMatrixi$productIon)
        intensityOfproductIon = c(intensityOfproductIon,allIonPairsMatrixi$intensity)
        neutralLoss = c(neutralLoss,allIonPairsMatrixi$neutralLoss)
        CE = c(CE,allIonPairsMatrixi$CE)
      }
    }
    close(pb)
  }
  allIonPairsMatrix <- data.frame(Index,RT,m.z,beginNum,trInSecond,pepmass,endNum,productIon,intensityOfproductIon,neutralLoss,CE)
  return(allIonPairsMatrix)
}
