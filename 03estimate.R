# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma", version = "3.8")
rm(list=ls())

library(limma)
library(estimate)


inputFile="05Cox/trainSet_riskScore.xls"
rt=read.table(inputFile,header=T,check.names=F,sep="\t")

OV_mRNA <- read.table("01rawData/TCGA-LUAD_FPKM_clean.xls",row.names = 1,quote = "",check.names = F,sep = "\t",header = T)
OV_mRNA[1:3,1:5]
colnames(OV_mRNA) <- substring(colnames(OV_mRNA),1,16)
low <- rt$id[rt$risk =="Low Risk"]
high <- rt$id[rt$risk =="High Risk"]
# rt1 <- log(rt+1)
low_high_df <- OV_mRNA[,c(low,high)]
low_high_out <- cbind(id = rownames(low_high_df),low_high_df)
write.table(low_high_out,"08estimate/trainSet_riskScore_low_high.xls",row.names = F,quote = F,sep="\t")


library(estimate)
# estimate
filterCommonGenes(input.f="08estimate/trainSet_riskScore_low_high.xls", 
                  output.f="08estimate/commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "08estimate/commonGenes.gct",
              output.ds="08estimate/estimateScore.gct", 
              platform="illumina")


scores=read.table("08estimate/estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(id=colnames(scores),scores)
write.table(out,file="08estimate/estimateScores.xls",sep="\t",quote=F,col.names=F)

# 
# library(pheatmap)
# rt=read.csv("immunesurvival.csv",header=T,sep=",",check.names=F)
# risk <- read.csv("riskScore.csv",header = T,sep = ",",check.names = F)
# Data <- merge(risk,rt[-c(2:3)],by="id")
# write.table(Data,file="heatmap_input.csv",sep=",",quote=F,col.names=T,row.names = F)
# Data=read.csv("heatmap_input.csv",header=T,sep=",",check.names=F)
# Data <- Data[order(Data$risk,decreasing = T ),]
# heatmap1 <- Data[,c(1,4:9)]
# rownames(heatmap1) <- heatmap1$id
# heatmap1 <- t(heatmap1[,c(2:7)])
# heatmap1 <-  log2(heatmap1+1)
# 
# # 构建列注释信息
# annotation_col = data.frame(
#   #Risk = factor(ifelse(Data$risk =="high","High","Low"),levels = c("Low","High")),
#   Risk = ifelse(Data$risk =="high","High","Low"),
#   Gender =ifelse(Data$gender =="FEMALE","Female","Male"),
#   Status = ifelse(Data$fustat==0,"Alive","Dead"),
#   Stage = Data$stage,
#   ImmuneScore = ifelse(Data$ImmuneScore > median(Data$ImmuneScore),"High","Low")
#   
# )
# row.names(annotation_col) <- Data$id
# # 自定注释信息的颜色列表
# ann_colors = list(
#   Risk = c(Low = "#2874A6", High = "#CB4335"),
#   Gender = c(Female = "#1B9E77", Male = "#D95F02"),
#   Status = c(Alive = "#00CC33",Dead = "#DC7633"),
#   Stage = c(`Stage I`="#00897B", `Stage II`="#689F38", `Stage III`="#AFB42B",`Stage IV`="#FF8F00"),
#   ImmuneScore = c(Low = "#2874A6", High = "#CB4335")
# )
# pdf(file="ImmuneScore_heatmap.pdf",width = 9,height = 4)
# pheatmap(heatmap1,
#          cluster_cols = FALSE,
#          show_colnames=F,
#          annotation_col = annotation_col,annotation_colors =ann_colors,
#          color = colorRampPalette(c("#028846", "black", "red"))(50) )
# dev.off()

###########  画免疫分数在高低风险组的箱线图
library(tidyr)
library(ggpubr)
risk <- read.table("05Cox/trainSet_riskScore.xls",header = T,sep = "\t",check.names = F)
scores =read.table("08estimate/estimateScores.xls",header=T,sep="\t",check.names=F)
Data <- merge(risk,scores,by="id")
plotData <- Data[c("risk","StromalScore","ImmuneScore","ESTIMATEScore")]
plotData$risk <- factor(plotData$risk,levels = c("Low Risk","High Risk"))
# plotData <- gather(plotData,Immune,Score,StromalScore,ImmuneScore)
# plotData$risk <- ifelse(plotData$risk =="high","High_RS","Low_RS")
p1 <- ggviolin(plotData, x = "risk", y = "StromalScore",
               fill  = "risk", palette = c("#00AFBB","#E7B800"),
               add = c("boxplot"),add.params = list(fill="white"))+stat_compare_means(aes(group = risk)) +xlab("")

p2 <- ggviolin(plotData, x = "risk", y = "ImmuneScore",
               fill  = "risk", palette = c("#00AFBB","#E7B800"),
               add = c("boxplot"),add.params = list(fill="white"))+stat_compare_means(aes(group = risk)) +xlab("")
p3 <- ggviolin(plotData, x = "risk", y = "ESTIMATEScore",
               fill  = "risk", palette = c("#00AFBB","#E7B800"),
               add = c("boxplot"),add.params = list(fill="white"))+stat_compare_means(aes(group = risk)) +xlab("")
# 使用ggpubr包的函数ggarrange（）中在一页上进行组合展示
ggarrange(p1,p2,p3,labels = c("A", "B","C"),ncol = 3, nrow = 1,common.legend = T)

ggsave("08estimate/Score_violin.pdf",width = 12,height = 4)



