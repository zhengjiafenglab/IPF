
## 加载R包
# library(GEOquery)
library(limma)
library(dplyr)
library(tibble)
#高低风险组差异基因分析

sampleId <- read.table("05Cox/trainSet_riskScore.xls",sep="\t",header=T,check.names=F)
table(sampleId$risk)

data <- read.table("07CIBERSORT/GSE70866_removeBatchEffect_expres.txt",row.names = 1,check.names = F,sep = "\t",header = T)

lowId = colnames(fe_expr)[  colnames(data) %in% sampleId$id[sampleId$risk =="Low Risk"]] #肿瘤组
highId = colnames(fe_expr)[ colnames(data) %in% sampleId$id[sampleId$risk =="High Risk"]] #肿瘤组

GSEA_df <- data[,c(lowId,highId)]
# 292low   292high 
# rt1 <- cbind(id = sapply(rownames(rt1),function(x) strsplit(x,split = "\\|")[[1]][2]),rt1)
GSEA_df[1:3,1:2]
GSEA_df <- cbind(Name = rownames(GSEA_df) ,DESCRIPTION = "na",GSEA_df)

group <- c(rep("Low_Risk", 117), rep("High_Risk", 59))
group <- paste(group, collapse = " ")
group <- c(paste(c(176, 2, 1), collapse = " "), "# Low_Risk High_Risk", group)


write.table(file = "09GSEA/GSE70866_risk_group_exp.txt", GSEA_df, sep = "\t", col.names = T, row.names = F, quote = F)

write.table(file = "09GSEA/risk_group.cls", group, col.names = F, row.names = F, quote = F)
