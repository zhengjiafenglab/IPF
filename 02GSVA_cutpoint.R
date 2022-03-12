
## 加载R包
# library(GEOquery)
library(limma)
library(dplyr)
library(tibble)
library(pheatmap)
library(Rtsne)
library(GSVA)
exprSet <-  read.table("01rawData/GSE70866_removeBatchEffect_expres.xls",head = TRUE,row.names=1,sep="\t")
# set.seed(9)
# tsne<- Rtsne(t(exprSet), dims = 2, perplexity=9, verbose=TRUE, max_iter = 500)
# 
# #二次聚类，因为tsne仅仅是更好的展示聚类的结果
# set.seed(9)
# cl_i <- kmeans(tsne$Y, centers = 2, nstart =10,algorithm="Lloyd")
# data1 <- data.frame(colnames(exprSet), tsne$Y,cl_i$cluster)
# colnames(data1) <- c("id", 'Coordinate1','Coordinate2','Cluster')
# # data1$Cluster <- ifelse(data1$Cluster==1,'Cluster2','Cluster1')
# data1$Cluster <- ifelse(data1$Cluster==1,'Cluster1','Cluster2')
# ggplot(data1,aes(Coordinate1,Coordinate2,colour=Cluster))+
#   geom_point(size=3.2)+scale_colour_manual(values =  c("#FF0033","#1B9E77"))+
#   theme_bw()+theme(panel.grid =element_blank())+theme(axis.title = element_text(size=18),
#                                                       axis.text = element_text(size=18),
#                                                       legend.text=element_text(size=18),
#                                                       legend.title = element_text(size=18),
#                                                       plot.title = element_text(size = 20,hjust = 0.5))
# ggsave("02GSVA/sample_cluster.pdf",width = 6.5,height = 5)

gs <- read.table("01rawData/EMT_Genes.csv",sep = ",",check.names = F,header = TRUE)

exprSet_FerrDb <- exprSet[rownames(exprSet) %in%  gs$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,]
# exprSet_FerrDb <- as.matrix(exprSet_FerrDb)
# #一步完成聚类
# library(ConsensusClusterPlus)
# title="F:\\November\\YQ189-8"
# results = ConsensusClusterPlus(exprSet_FerrDb,maxK=6,reps=50,pItem=0.8,pFeature=1,
#                                title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
pheatmap(exprSet_FerrDb,scale = "row")


gs <- as.list(gs)

# 这一句就完成了GSVA分析
gsva_es <- gsva(as.matrix(exprSet), gs)

# 预览GSVA分析返回的矩阵
head(gsva_es)
gsva_es_df <-  as.data.frame(t(gsva_es))
gsva_es_df <- cbind(id=rownames(gsva_es_df),gsva_es_df)
colnames(gsva_es_df)[2] <- "score"
write.table(gsva_es_df,file = "02GSVA/GSE70866_EMT_score.xls",row.names = F,quote = F,sep = "\t")



clincal_raw <- read.delim("01rawData/GSE70866_clinical.txt",check.names = F,header = T,sep = "\t")


fe_expr <- read.table("02GSVA/GSE70866_EMT_score.xls",head = TRUE,row.names=1,sep="\t")
# group <- substr( rownames(fe_expr),14,16)
# table(group)
fe_expr$id <- rownames(fe_expr)
# fe_expr$id <- substr(rownames(fe_expr),1,12)
fe_exprTime <- merge(clincal_raw,fe_expr,by="id")
# fe_exprTime <- fe_exprTime[,c(ncol(fe_exprTime),2:(ncol(fe_exprTime)-1))]
write.table(fe_exprTime,file = "02GSVA/GSE70866_EMT_optimal_survial_input.xls",row.names = F,quote = F,sep = "\t")

# 寻找最佳生存率
library(survminer)
# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(fe_exprTime, time = "futime", event = "fustat",minprop = 0.1,
                         variables = "score")

# 2. Plot cutpoint for DEPDC1
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
pdf("02GSVA/optimal_cutpoint.pdf",width = 6.5,height = 5)
# plot(res.cut, "score", palette = "npg")
plot(res.cut, "score", palette =  c("#FF0033","#1B9E77"))
dev.off()
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
res.cat_out <-cbind( id =fe_exprTime$id,res.cat) 
write.table(res.cat_out,file = "02GSVA/GSE70866_EMT_optimal_survial_group.xls",row.names = F,quote = F,sep = "\t")
# 4. Fit survival curves and visualize
library("survival")
fit <- survfit(Surv(futime, fustat) ~score, data = res.cat)
pdf("02GSVA/optimal_cutpoint_survfit.pdf",width = 5,height = 6)
ggsurvplot(fit, data = res.cat, 
           palette = c("#FF0033","#1B9E77"),
           xlab = "Time (Days)",
           risk.table = TRUE, 
           conf.int = TRUE,
           pval=TRUE)
dev.off()


library(umap)
library(Rcpp)
library(factoextra)
dat <- read.table("01rawData/GSE70866_removeBatchEffect_expres.xls",head = TRUE,row.names=1,sep="\t")
dat[1:4,1:4]
gs <- read.table("01rawData/EMT_Genes.csv",sep = ",",check.names = F,header = TRUE)
dat_Autophagy <- dat[rownames(dat) %in% gs$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,]

groups <- read.table("02GSVA/GSE70866_EMT_optimal_survial_group.xls",header = TRUE,sep="\t",row.names = 1)
# groups <- groups[colnames(dat),3]
group<-groups[match( colnames(dat),rownames(groups) ),c(2,3)]
# group$score <- factor(group$score,levels =c("low","high") )
group$labels <- ifelse(group$score =="low","#1B9E77","#FF0033")
iris.umap = umap(t(dat_Autophagy))
pdf("02GSVA/EMTGenes_score_umap.pdf",width = 6,height = 6)
plot(iris.umap$layout,col=group$labels,pch=16,asp = 1,
     xlab = "UMAP_1",ylab = "UMAP_2")
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")
# 添加图例
legend("topright",title = "EMT score",#bty = "n",
       legend = unique(groups$score),pch=16,
       col =  c("#1B9E77","#FF0033"))

dev.off()

# g=factor( groups$score )
# g
# pdf(file="02GSVA/pca_score.pdf",height=5,width=6.5)
# pca_plot(dat,g)
# dev.off()


Normal =  rownames(groups)[groups$score =="low"] #正常组
Tumor = rownames(groups)[groups$score =="high"] #肿瘤组

# `High stage`=sample[grepl("Ovarian carcinoma tissue",as.character(sample$source_name_ch1)),]$geo_accession#肿瘤组


group_list=c(rep("Normal",length(Normal)),
             rep("Tumor",length(Tumor)))  #分组信息
# group_list=factor(group_list,levels = c("Low stage","High stage"))
group_list=factor(group_list)
# ## 强制限定顺序
group_list <- relevel(group_list, ref="Normal")
table(group_list)
# group_list
# Normal  Tumor 
# 84     92 

exprSet <- dat[,c(Normal,Tumor)]
# 保存在01rawData，分析数据
# exprSet_out <- cbind(gene_id=rownames(exprSet),exprSet)
# exprSet_out[1:3,1:5]
# write.table(exprSet_out,file = "03estimate/GSE70866_ImmuneScore_exprset.xls",row.names = F,quote = F,sep = "\t")

# 4.判断是否需要进行数据转换
######################################
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

# 过去掉低表达的基因
exprSet=exprSet[rowMeans(exprSet)>1,]
# exprSet1 <- log2(exprSet)
exprSet1 <- exprSet
library(limma)
# 7.差异分析
#############################################
# 我们选取格式比较简单的。如果没有配对信息，差异分析这样做：
design=model.matrix(~ group_list)
fit=lmFit(exprSet1,design)
fit=eBayes(fit) 
res=topTable(fit,adjust='fdr',coef="group_listTumor",number=Inf)
res <- na.omit(res)

logFCCutoff <- 1
pvalueCutoff <- 0.05

write.table(res,file="02GSVA/GSE70866_EMTScore_all.xls",row.names=T,quote=F,sep = "\t",col.names = NA)
res.lncRNA <- res
outDiff.lncRNA =res.lncRNA[(abs(res.lncRNA$logFC)>logFCCutoff & res.lncRNA$adj.P.Val<pvalueCutoff),]
write.table(outDiff.lncRNA,file="02GSVA/GSE70866_EMTScore_diff.xls",row.names=T,quote=F,sep = "\t",col.names = NA)


outDiff <- outDiff.lncRNA
#绘制差异基因热图
library(pheatmap)
geneNum=50
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=rownames( outDiff )
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp= exprSet[hmGene,]
# hmExp_out <- cbind(id=rownames(hmExp),hmExp) 
# hmExp_out[1:3,1:5]
# write.csv(hmExp_out,file = "01rawData/TCGA-THCA_DEmRNA_FPKM_exprSet.csv",row.names = F,quote = F)

# hmExp = log2(hmExp+1)
# hmExp[hmExp > 8] =8
# hmExp[hmExp < 8] =8
max(hmExp)
min(hmExp)

# rownames(hmExp) <- gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\",rownames(hmExp))
# Normal  Tumor 
# 103    367 
`EMT Score`= factor(c(rep("low",84),rep("high",92)),levels = c("low","high"))
names(`EMT Score`)=colnames(hmExp)
`EMT Score`=as.data.frame(`EMT Score`)

# loc <- order(Type,colSums(hmExp),decreasing = T)
pdf(file="02GSVA/EMT_mRNA_heatmap.pdf",height=8,width=9)
pheatmap(hmExp, 
         annotation=`EMT Score`, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8,
         border=FALSE)
dev.off()




#定义显著性
library(ggplot2)
library(ggrepel)
res <- res.lncRNA
# labelGene <- read.table("01rawData/m6A_genes.tsv",header = F,sep="\t")$V1

# res$label <- ifelse(gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\",rownames(res)) %in% labelGene,gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\",rownames(res)),"" )

Significant=ifelse((res$adj.P.Val< 0.05 & abs(res$logFC)> 1), ifelse(res$logFC > 1,"Up","Down"), "Not")
#绘制火山图
p = ggplot(res, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#004BFB", "#BFBFBF", "#F91F10"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-1,1), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+xlab("log2FC")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+theme_bw()

# p1 <- p +geom_text_repel( data= res ,aes(logFC, -log10(adj.P.Val),label = label),
#                           size = 3,box.padding = unit(0.5, "lines"),
#                           point.padding = unit(0.8, "lines"), 
#                           segment.color = "black", 
#                           show.legend = FALSE)

#保存为图片
pdf("02GSVA/EMT_mRNA_vol.pdf",width=5.5,height=5)
print(p)
dev.off()




library(estimate)
# estimate
filterCommonGenes(input.f="01rawData/GSE70866_removeBatchEffect_expres.xls", 
                  output.f="03estimate/commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "03estimate/commonGenes.gct",
              output.ds="03estimate/estimateScore.gct", 
              platform="agilent")


scores=read.table("03estimate/estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
# rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(id=colnames(scores),scores)
write.table(out,file="03estimate/estimateScores.xls",sep="\t",quote=F,col.names=F)


clincal_raw <- read.delim("01rawData/GSE70866_clinical.txt",check.names = F,header = T,sep = "\t")


fe_expr <- read.table("03estimate/estimateScores.xls",head = TRUE,row.names=1,sep="\t")
# group <- substr( rownames(fe_expr),14,16)
# table(group)
fe_expr$id <- rownames(fe_expr)
# fe_expr$id <- substr(rownames(fe_expr),1,12)
fe_exprTime <- merge(clincal_raw,fe_expr,by="id")
# fe_exprTime <- fe_exprTime[,c(ncol(fe_exprTime),2:(ncol(fe_exprTime)-1))]
write.table(fe_exprTime,file = "03estimate/GSE70866_estimate_optimal_survial_input.xls",row.names = F,quote = F,sep = "\t")

# 寻找最佳生存率
library(survminer)
# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(fe_exprTime, time = "futime", event = "fustat",minprop = 0.1,
                         variables = "ImmuneScore")

# 2. Plot cutpoint for DEPDC1
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
pdf("03estimate/ImmuneScore_optimal_cutpoint.pdf",width = 6,height = 6)
# plot(res.cut, "score", palette = "npg")
plot(res.cut, "ImmuneScore", palette =  c("#FF0033","#1B9E77"))
dev.off()
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
res.cat_out <-cbind( id =fe_exprTime$id,res.cat) 
write.table(res.cat_out,file = "03estimate/GSE70866_ImmuneScore_optimal_survial_group.xls",row.names = F,quote = F,sep = "\t")
# 4. Fit survival curves and visualize
library("survival")
fit <- survfit(Surv(futime, fustat) ~ImmuneScore, data = res.cat)
pdf("03estimate/ImmuneScore_optimal_cutpoint_survfit.pdf",width = 5,height = 6)
ggsurvplot(fit, data = res.cat, 
           palette = c("#FF0033","#1B9E77"),
           xlab = "Time (Days)",
           risk.table = TRUE, 
           conf.int = TRUE,
           pval=TRUE)
dev.off()




library(umap)
library(Rcpp)
library(factoextra)
dat <- read.table("01rawData/GSE70866_removeBatchEffect_expres.xls",head = TRUE,row.names=1,sep="\t")
dat[1:4,1:4]

groups <- read.table("03estimate/GSE70866_ImmuneScore_optimal_survial_group.xls",header = TRUE,sep="\t",row.names = 1)
# groups <- groups[colnames(dat),3]
group<-groups[match( colnames(dat),rownames(groups) ),c(2,3)]
# group$score <- factor(group$score,levels =c("low","high") ,labels = c("#1B9E77","#FF0033") )
group$labels <- ifelse(group$ImmuneScore =="low","#1B9E77","#FF0033")
iris.umap = umap(t(dat))
pdf("03estimate/ImmuneScore_umap.pdf",width = 6,height = 6)
plot(iris.umap$layout,col=group$labels,pch=16,asp = 1,
     xlab = "UMAP_1",ylab = "UMAP_2")
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")
# 添加图例
legend("topright",title = "ImmuneScore",
       legend = unique(groups$ImmuneScore),pch=16,
       col =  c("#1B9E77","#FF0033"))

dev.off()


Normal =  rownames(groups)[groups$ImmuneScore =="low"] #正常组
Tumor = rownames(groups)[groups$ImmuneScore =="high"] #肿瘤组

# `High stage`=sample[grepl("Ovarian carcinoma tissue",as.character(sample$source_name_ch1)),]$geo_accession#肿瘤组


group_list=c(rep("Normal",length(Normal)),
             rep("Tumor",length(Tumor)))  #分组信息
# group_list=factor(group_list,levels = c("Low stage","High stage"))
group_list=factor(group_list)
# ## 强制限定顺序
group_list <- relevel(group_list, ref="Normal")
table(group_list)
# group_list
# Normal  Tumor 
# 67    109 

exprSet <- dat[,c(Normal,Tumor)]
# 保存在01rawData，分析数据
# exprSet_out <- cbind(gene_id=rownames(exprSet),exprSet)
# exprSet_out[1:3,1:5]
# write.table(exprSet_out,file = "03estimate/GSE70866_ImmuneScore_exprset.xls",row.names = F,quote = F,sep = "\t")

# 4.判断是否需要进行数据转换
######################################
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

# 过去掉低表达的基因
exprSet=exprSet[rowMeans(exprSet)>1,]
# exprSet1 <- log2(exprSet)
exprSet1 <- exprSet
library(limma)
# 7.差异分析
#############################################
# 我们选取格式比较简单的。如果没有配对信息，差异分析这样做：
design=model.matrix(~ group_list)
fit=lmFit(exprSet1,design)
fit=eBayes(fit) 
res=topTable(fit,adjust='fdr',coef="group_listTumor",number=Inf)
res <- na.omit(res)

logFCCutoff <- 1
pvalueCutoff <- 0.05

write.table(res,file="03estimate/GSE70866_ImmuneScore_all.xls",row.names=T,quote=F,sep = "\t",col.names = NA)
res.lncRNA <- res
outDiff.lncRNA =res.lncRNA[(abs(res.lncRNA$logFC)>logFCCutoff & res.lncRNA$adj.P.Val<pvalueCutoff),]
write.table(outDiff.lncRNA,file="03estimate/GSE70866_ImmuneScore_diff.xls",row.names=T,quote=F,sep = "\t",col.names = NA)
# 
# #绘制差异lncRNA热图
# library(pheatmap)
# geneNum=50
# outDiff.lncRNA=outDiff.lncRNA[order(as.numeric(as.vector(outDiff.lncRNA$logFC))),]
# diffGeneName=rownames( outDiff.lncRNA )
# diffLength=length(diffGeneName)
# hmGene=c()
# if(diffLength>(geneNum*2)){
#   hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
# }else{
#   hmGene=diffGeneName
# }
# hmExp= exprSet[hmGene,]
# hmExp_out <- cbind(id=rownames(hmExp),hmExp) 
# hmExp_out[1:3,1:5]
# write.csv(hmExp_out,file = "01rawData/TCGA-THCA_DElncRNA_FPKM_exprSet.csv",row.names = F,quote = F)
# 
# hmExp = log2(hmExp+1)
# # hmExp[hmExp > 8] =8
# # hmExp[hmExp < 8] =8
# max(hmExp)
# min(hmExp)
# 
# rownames(hmExp) <- gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\",rownames(hmExp))
# # Normal  Tumor 
# # 103    367 
# Type= factor(c(rep("Normal",58),rep("PTC",381)),levels = c("Normal","PTC"))
# names(Type)=colnames(hmExp)
# Type=as.data.frame(Type)
# 
# # loc <- order(Type,colSums(hmExp),decreasing = T)
# pdf(file="02DEGs/lncRNA_heatmap.pdf",height=8,width=9)
# pheatmap(hmExp, 
#          annotation=Type, 
#          color = colorRampPalette(c("blue", "white", "red"))(50),
#          cluster_cols =F,
#          show_colnames = F,
#          scale="row",
#          fontsize = 8,
#          fontsize_row=6,
#          fontsize_col=8,
#          border=FALSE)
# dev.off()
# 
# 
# #定义显著性
# library(ggplot2)
# 
# Significant=ifelse((res.lncRNA$adj.P.Val< 0.05 & abs(res.lncRNA$logFC)> 0.5), ifelse(res.lncRNA$logFC > 0.5,"Up","Down"), "Not")
# #绘制火山图
# p = ggplot(res.lncRNA, aes(logFC, -log10(adj.P.Val)))+
#   geom_point(aes(col=Significant))+
#   scale_color_manual(values=c("#004BFB", "#BFBFBF", "#F91F10"))+
#   labs(title = " ")+
#   geom_vline(xintercept=c(-0.5,0.5), colour="black", linetype="dashed")+
#   geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
#   theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+xlab("log2FC")
# p=p+theme_bw()
# #保存为图片
# pdf("02DEGs/lncRNA_vol.pdf",width=5.5,height=5)
# print(p)
# dev.off()


outDiff <- outDiff.lncRNA
#绘制差异基因热图
library(pheatmap)
geneNum=50
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=rownames( outDiff )
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp= exprSet[hmGene,]
# hmExp_out <- cbind(id=rownames(hmExp),hmExp) 
# hmExp_out[1:3,1:5]
# write.csv(hmExp_out,file = "01rawData/TCGA-THCA_DEmRNA_FPKM_exprSet.csv",row.names = F,quote = F)

# hmExp = log2(hmExp+1)
# hmExp[hmExp > 8] =8
# hmExp[hmExp < 8] =8
max(hmExp)
min(hmExp)

rownames(hmExp) <- gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\",rownames(hmExp))
# Normal  Tumor 
# 103    367 
Immune= factor(c(rep("low",67),rep("high",109)),levels = c("low","high"))
names(Immune)=colnames(hmExp)
Immune=as.data.frame(Immune)

# loc <- order(Type,colSums(hmExp),decreasing = T)
pdf(file="03estimate/immune_mRNA_heatmap.pdf",height=8,width=9)
pheatmap(hmExp, 
         annotation=Immune, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8,
         border=FALSE)
dev.off()


#定义显著性
library(ggplot2)
library(ggrepel)
res <- res.lncRNA
# labelGene <- read.table("01rawData/m6A_genes.tsv",header = F,sep="\t")$V1

# res$label <- ifelse(gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\",rownames(res)) %in% labelGene,gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\",rownames(res)),"" )

Significant=ifelse((res$adj.P.Val< 0.05 & abs(res$logFC)> 1), ifelse(res$logFC > 1,"Up","Down"), "Not")
#绘制火山图
p = ggplot(res, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#004BFB", "#BFBFBF", "#F91F10"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-1,1), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+xlab("log2FC")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+theme_bw()

# p1 <- p +geom_text_repel( data= res ,aes(logFC, -log10(adj.P.Val),label = label),
#                           size = 3,box.padding = unit(0.5, "lines"),
#                           point.padding = unit(0.8, "lines"), 
#                           segment.color = "black", 
#                           show.legend = FALSE)

#保存为图片
pdf("03estimate/immune_mRNA_vol.pdf",width=5.5,height=5)
print(p)
dev.off()







