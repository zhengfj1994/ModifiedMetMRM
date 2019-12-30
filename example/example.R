# ModifiedMetMRM使用范例和说明
# 第一步：安装ModifiedMetMRM包

# 第二步：调用
library(ModifiedMetMRM)

# 第三步：依次运行主程序

# 创建所有的离子对，储存在矩阵中，MS1File是一级文件，txt格式；MS2path是存储二级文件的文件夹；tolMZppm是质荷比偏差，tolTr是保留时间偏差，单位是分钟。
allIonPairsMatrix <- createAllIonPairs(MS1file="D:\\github\\ModifiedMetMRM\\示例文件\\POS_一级峰表\\POS-TIC.txt",
                                       MS2path="D:\\github\\ModifiedMetMRM\\示例文件\\POS_二级文件",
                                       tolMZppm=20,
                                       tolTr=0.2)

# 根据中性丢失进行修饰代谢物的定性，需要上一步的结果，和修饰种类列表，和质荷比偏差。
modifiedIdentificationRes <- modifiedIdentification(allIonPairsMatrix,
                                                    modificationList="D:\\github\\ModifiedMetMRM\\示例文件\\Modification List.xlsx",
                                                    tolMZppm=20)

# 对结果进行离子融合
afterIonFusionRes <- ionFusion(rawMatrix=modifiedIdentificationRes,
                               tolMZppm=20,
                               tolTr=0.2,
                               ionMode="positive")

# 写出结果，结果中会包含重复值，因为存在多个不同的碰撞能。
write.csv(afterIonFusionRes, file="D:\\github\\ModifiedMetMRM\\Identification result.csv")

# 根据最终定性结果，筛选mgf文件
mgfFilter(identificationResult=afterIonFusionRes,
          rawMS2path="D:\\github\\ModifiedMetMRM\\示例文件\\POS_二级文件",
          resultMS2path="D:\\github\\ModifiedMetMRM\\示例文件\\筛选结果")

# MRM-Ion Pair Finder
MRMIonPairFinder(MS1File="D:\\github\\ModifiedMetMRM\\Identification result.csv",
                 MS2path="D:\\github\\ModifiedMetMRM\\示例文件\\筛选结果",
                 tolMZppm=20,
                 tolTr=0.2,
                 diffMS1MS2=-1500,
                 MS2Int=0,
                 resultPath="D:\\github")

