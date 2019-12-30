#' @name MRMIonPairFinder
#' @description MRM-Ion Pair Finder performed in R
#' @author Fujian Zheng <zhengfj@dicp.ac.cn>
#' @param MS1File MS1 peak detection result save in .csv filetype, the first column is m/z named 'mz', the second column is retention time(s) named 'tr', intensity of samples is located begin the third column.
#' @param MS2path The folder path which have mgf files.
#' @param tolMZppm: The tolerence of m/z between MS1 peak detection result and mgf files. 0.01 is suitable for Q-TOF.
#' @param tolTr(min): The tolerence of retention time between MS1 peak detection result and mgf files.
#' @param diffMS1MS2(Da): The smallest difference between product ion and precusor ion.
#' @param MS2Int: The smallest intensity of product ion.
#' @param resultPath: A csv file named "MRM transitions list.csv" will saved in the path.
#' @return
#' @example MRMIonPairFinder(MS1File="D:\\github\\ModifiedMetMRM\\Identification result.csv", MS2path="D:\\github\\ModifiedMetMRM\\示例文件\\筛选结果", tolMZppm=20, tolTr=0.2, diffMS1MS2=0, MS2Int=0, resultPath="D:\\github")
#' @references Analytical Chemistry 87.10(2015):5050-5055.

MRMIonPairFinder <- function(MS1File, MS2path, tolMZppm, tolTr, diffMS1MS2, MS2Int, resultPath){
  # Some packages used in the function
  ##########
  require(tcltk)
  require(stringr)
  require(readr)
  require(dplyr)
  ##########

  ##########
  before_pretreatment <- read.csv(file = MS1File)
  mz <- before_pretreatment$m.z
  tr <- before_pretreatment$RT
  ##########

  MS2_filename <- list.files(MS2path)
  data_ms1ms2 <- cbind(before_pretreatment[1,], mzinmgf=1, trinmgf=1, mz_ms2=1, int_ms2=1, CE=1)[-1,]  # Create data.frame to store information of ms1ms2 information

  # Reading and processing mgf files one by one.
  for (i_new in MS2_filename){
    mgfData <- scan(paste0(MS2path,'\\',i_new), what = character(0), sep = "\n")
    mgfMatrix <- createMgfMatrix(mgfData)  # create mgfMatrix
    CE <- parse_number(i_new) # get CE value in the filename of mgf

    # Delete the data with charge > 1
    ##########
    pb <- tkProgressBar(paste("Delete the data in", i_new, "with charge > 1"),"rate of progress %", 0, 100)
    for (i in c(1:length(mgfData))){
      info<- sprintf("rate of progress %d%%", round(i*100/length(mgfData)))
      setTkProgressBar(pb, i*100/length(mgfData), sprintf(paste("Delete the data in", i_new, "with charge > 1 (%s)"), info),info)
      # If the row of mgfData is contain the "CHARGE=",
      if (!is.na(mgfData[i]) & str_detect(mgfData[i],"CHARGE=")){
        if (!str_detect(mgfData[i],"CHARGE=1")){
          mgfData[mgfMatrix[tail(which(as.numeric(mgfMatrix[,"beginNum"]) < i),1),"beginNum"]:
                     mgfMatrix[which(as.numeric(mgfMatrix[,"endNum"]) > i)[1],"endNum"]] <- NA
        }
      }
    }
    close(pb)
    mgfData <- na.omit(mgfData)
    packageStartupMessage(paste("Deleting the data in", i_new, "with charge > 1 is finished."))
    ########

    # Delete the data by diffMS1MS2 and MS2Int
    ########
    mgfMatrix <- createMgfMatrix(mgfData)  # create mgfMatrix
    pb <- tkProgressBar(paste("Delete the data in", i_new, "by diffMS1MS2 and MS2Int"),"rate of progress %", 0, 100)
    for (i in c(1:length(mgfData))){
      info<- sprintf("rate of progress %d%%", round(i*100/length(mgfData)))
      setTkProgressBar(pb, i*100/length(mgfData), sprintf(paste("Delete the data in", i_new, "by diffMS1MS2 and MS2Int (%s)"), info),info)
      if (!grepl("[a-zA-Z]", mgfData[i])){
        mz_ms2 <- as.numeric(unlist(strsplit(mgfData[i], " "))[1])
        int_ms2 <- as.numeric(unlist(strsplit(mgfData[i], " "))[2])
        mz_ms1 <- as.numeric(mgfMatrix[tail(which(as.numeric(mgfMatrix[,"beginNum"]) < i),1),"pepmass"])
        if (mz_ms1-mz_ms2 <= diffMS1MS2 | int_ms2 <= MS2Int){
          mgfData[i] <- NA
        }
      }
    }
    close(pb)
    mgfData <- na.omit(mgfData)
    packageStartupMessage(paste("Deleting the data in", i_new, "by diffMS1MS2 and MS2Int is finished."))
    ########

    # Delete the data without useful MS2
    ########
    mgfMatrix <- as.data.frame(createMgfMatrix(mgfData))  # creat mgfMatrix
    pb <- tkProgressBar(paste("Delete the data in", i_new, "without useful MS2"),"rate of progress %", 0, 100)
    for (i in c(1:nrow(mgfMatrix))){
      info<- sprintf("rate of progress %d%%", round(i*100/nrow(mgfMatrix)))
      setTkProgressBar(pb, i*100/nrow(mgfMatrix), sprintf(paste("Delete the data in", i_new, "without useful MS2 (%s)"), info),info)
      if (as.numeric(as.character(mgfMatrix$endNum[i])) - as.numeric(as.character(mgfMatrix$beginNum[i])) < 5){
        mgfData[as.numeric(as.character(mgfMatrix$beginNum[i])):as.numeric(as.character(mgfMatrix$endNum[i]))] <- NA
      }
    }
    close(pb)
    mgfData <- na.omit(mgfData)
    packageStartupMessage(paste("Deleting the data in", i_new, "without useful MS2 is finished."))
    ########

    # Combine ms1 and ms2
    mgfMatrix <- as.data.frame(createMgfMatrix(mgfData))

    pb <- tkProgressBar("Combine ms1 and ms2","rate of progress %", 0, 100)
    for (i in c(1:nrow(mgfMatrix))){
      info<- sprintf("rate of progress %d%%", round(i*100/nrow(mgfMatrix)))
      setTkProgressBar(pb, i*100/nrow(mgfMatrix), sprintf("Combine ms1 and ms2", info),info)
      mzinmgf <- as.numeric(as.character(mgfMatrix$pepmass[i]))
      trinmgf <- as.numeric(as.character(mgfMatrix$trInSecond[i]))
      posi <- which(abs(before_pretreatment$m.z-mzinmgf)/before_pretreatment$m.z*1000000 < tolMZppm & abs(before_pretreatment$tr-trinmgf) < tolTr*60)
      if (length(posi)==1){
        ms1info <- before_pretreatment[posi,]
        for (j in mgfData[as.numeric(as.character(mgfMatrix$beginNum[i])):as.numeric(as.character(mgfMatrix$endNum[i]))]){
          if (grepl("[a-zA-Z]", j)){
            next()
          }else{
            mz_ms2 <- as.numeric(unlist(strsplit(j, " "))[1])
            int_ms2 <- as.numeric(unlist(strsplit(j, " "))[2])
            ms1ms2conb <- cbind(ms1info,mzinmgf,trinmgf,mz_ms2,int_ms2,CE)
            data_ms1ms2 <- rbind(data_ms1ms2,ms1ms2conb)
          }
        }
      }
      if (length(posi)>1){
        for (k in c(1:length(posi))){
          posk <- posi[k]
          ms1info <- before_pretreatment[posk,]
          for (j in mgfData[as.numeric(as.character(mgfMatrix$beginNum[i])):as.numeric(as.character(mgfMatrix$endNum[i]))]){
            if (grepl("[a-zA-Z]", j)){
              next()
            }else{
              mz_ms2 <- as.numeric(unlist(strsplit(j, " "))[1])
              int_ms2 <- as.numeric(unlist(strsplit(j, " "))[2])
              ms1ms2conb <- cbind(ms1info,mzinmgf,trinmgf,mz_ms2,int_ms2,CE)
              data_ms1ms2 <- rbind(data_ms1ms2,ms1ms2conb)
            }
          }
        }
      }
    }
  }

  data_ms1ms2_final <- data_ms1ms2[0,]
  # uniquedata_ms1ms2 <- distinct(data_ms1ms2[,1:ncol(before_pretreatment)])
  for (i in c(1:nrow(before_pretreatment))){
    posi <- which(data_ms1ms2$m.z==before_pretreatment$m.z[i] & data_ms1ms2$RT==before_pretreatment$RT[i])
    if (length(posi) >= 1){
      temp <- data_ms1ms2[posi,]
      posi <- which(temp$int_ms2 == max(temp$int_ms2))
      temp <- temp[posi[1],]
      data_ms1ms2_final <- rbind(data_ms1ms2_final,temp)
    }
    else{
      temp <- cbind(before_pretreatment[i,], mzinmgf=0, trinmgf=0, mz_ms2=0, int_ms2=0, CE=0)
      data_ms1ms2_final <- rbind(data_ms1ms2_final,temp)
    }
  }
  setwd(resultPath)
  write.csv(data_ms1ms2_final,file = "MRM transitions list.csv",row.names = FALSE)
  return(data_ms1ms2_final)
}
