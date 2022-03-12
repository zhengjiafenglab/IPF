rm(list=ls())
options(stringsAsFactors = F)
## 加载R包
# library(GEOquery)
library("ggvenn")

methyltransferase <-  read.table("02GSVA/GSE70866_EMTScore_diff.xls",check.names = F,sep="\t",header = T)
DEGs <-  read.table("03estimate/GSE70866_ImmuneScore_diff.xls",check.names = F,sep="\t",header = T)
# DEGs <-   gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\", DEGs[,1])

ls <- list(`EMT-DEGs` = methyltransferase[,1] ,`Immune-DEGs`  = DEGs[,1])
ggvenn(ls ,     
       fill_color = c("blue", "yellow"),#填充颜色      
       fill_alpha = 0.3,#改变图形透明度       
       stroke_color="black",#边界线的颜色
       stroke_alpha = 0.2,#边界线的透明度  
       stroke_size = 0.2,
       # stroke_linetype="dashed",#将边界线转化为虚线       
       text_size = 4)#调节字体大小
ggsave(file="04venn/veen_intersectGenes.pdf",width = 5,height = 5)
intersectGenes <- intersect( methyltransferase[,1],DEGs[,1])
write.table(intersectGenes,file="04venn/veen_intersectGenes.xls",row.names=F,quote=F,sep = "\t",col.names = F)

dat <- read.table("01rawData/GSE70866_removeBatchEffect_expres.xls",head = TRUE,row.names=1,sep="\t")
dat[1:4,1:4]
exprSet_EMT <- dat[intersectGenes,]
# 保存在01rawData，分析数据
exprSet_EMT_out <- cbind(gene_id=rownames(exprSet_EMT),exprSet_EMT)
exprSet_EMT_out[1:3,1:5]
write.table(exprSet_EMT_out,file = "04venn/GSE70866_intersectGenes_exprset.xls",row.names = F,quote = F,sep = "\t")



#clinical
# clincal_raw <- read.delim("01rawData/TCGA-STAD_clinical.tsv",check.names = F,header = T,sep = "\t")
clincal_raw <- read.delim("01rawData/GSE70866_clinical.txt",check.names = F,header = T,sep = "\t")
# clincal <- data.frame(id = character(),`Overall Survival` = numeric(),`Survival status` = numeric(),`stage`=character(),stringsAsFactors = F )



exprSet <- read.table("04venn/GSE70866_intersectGenes_exprset.xls",row.names = 1,check.names = F,sep = "\t",header = T)
# NormalId <- colnames(exprSet)[-1][substr(colnames(exprSet)[-1],14,15)=="11"]
# TumorId <- colnames(exprSet)[substr(colnames(exprSet),14,15)=="01"]
# exprSet <- exprSet[,TumorId]

exprSet <- as.data.frame( t(exprSet))
exprSet$id <- rownames(exprSet)
# exprSet$id <- substr(rownames(exprSet),1,12)
exprSet_Time <- merge(clincal_raw,exprSet,by="id")
# exprSetTime <- exprSetTime[,c(ncol(exprSetTime),2:(ncol(exprSetTime)-1))]
write.table(exprSet_Time,file = "05Cox/GSE70866_intersectGenes_exprset_cox_input.xls",row.names = F,quote = F,sep = "\t")

