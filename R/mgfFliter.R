#' @name mgfFilter
#' @description Filter mgf by the identification result
#' @author Fujian Zheng <zhengfj@dicp.ac.cn>
#' @param identificationResult
#' @param rawMS2path
#' @param resultMS2path
#' @return None
#' @example mgfFilter(identificationResult=afterIonFusionRes, rawMS2path="D:\\github\\ModifiedMetMRM\\示例文件\\POS_二级文件", resultMS2path="D:\\github\\ModifiedMetMRM\\示例文件\\筛选结果")

mgfFilter <- function(identificationResult, rawMS2path, resultMS2path){
  library(tcltk)
  library(stringr)
  mgfFiles <- list.files(rawMS2path)
  for (mgfi in mgfFiles){
    mgfRawDatai <- scan(paste0(rawMS2path,"\\",mgfi), what = character(0), sep = "\n")
    CE <- parse_number(mgfi)
    filterMgfi <- mgfRawDatai[0]
    identificationResulti <- identificationResult[which(identificationResult$CE==CE), ]
    identificationResulti <- identificationResulti[!duplicated(identificationResulti$beginNum), ]

    pb <- tkProgressBar(paste0("mgfFilter for ",mgfi) ,"Rate of progress %", 0, 100)
    for (j in c(1:nrow(identificationResulti))){
      info<- sprintf("Rate of progress %d%%", round(j*100/nrow(identificationResulti)))
      setTkProgressBar(pb, j*100/nrow(identificationResulti), sprintf(paste0("mgfFilter for ",mgfi), info),info)

      beginNumj <- as.numeric(as.character(identificationResulti$beginNum[j]))
      endNumj <- as.numeric(as.character(identificationResulti$endNum[j]))
      filterMgfi <- rbind(filterMgfi, as.matrix(mgfRawDatai[beginNumj:endNumj]))
    }
    close(pb)
    write.table(filterMgfi,file=paste0(resultMS2path,"\\Filter ",mgfi), row.names=F,col.names = F, quote=F ,sep="\t")
  }
}
