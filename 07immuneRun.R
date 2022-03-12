#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(ggalluvial)
library(vioplot)
library(limma)
library(ggplot2)
pFilter=0.05 
source("CIBERSORT.R")
dat <- read.table("01rawData/GSE70866_removeBatchEffect_expres.xls",head = TRUE,row.names=1,sep="\t")
dat[1:4,1:4]
data_group= read.table("05Cox/trainSet_riskScore.xls",sep="\t",header=T,check.names=F)
Low <- data_group$id[data_group$risk == "Low Risk"]
HIGH <- data_group$id[data_group$risk == "High Risk"]
dat_out <- cbind(gene_id=rownames(dat),dat[,c(Low,HIGH )])
dat_out[1:3,1:5]
write.table(dat_out,file = "07CIBERSORT/GSE70866_removeBatchEffect_expres.txt",row.names = F,quote = F,sep = "\t")


# nperm给的是置换的次数，QN如果是芯片设置为T，如果是测序就设置为F，测序数据最好是TPM
results=CIBERSORT("07CIBERSORT/LM22.txt", "07CIBERSORT/GSE70866_removeBatchEffect_expres.txt", perm=100, QN=TRUE)

pFilter <- 0.05
#读取免疫结果文件，并对数据进行整理
immune=read.table("07CIBERSORT/CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

data=t(immune)
# col=rainbow(nrow(data),s=0.7,v=0.7)
col <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999","#ADD1E5",
         "#F0027F","#6A3D9A","#00FFFF","#808000","#0000FF","#CC33FF","#00FF33","#FFD9EC","#181717")



#生成图形
pdf(file = "07CIBERSORT/immuneCell_barplot.pdf",height=8,width=14)
par(las=1,mar=c(2,5,6,12),mgp=c(3,0.1,0),cex.axis=1.5)
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,tick=F,labels=F)
# par(srt=60,xpd=T);text(a1,-0.02,colnames(data),adj=1.1,cex=0.6);par(srt=0)
par(srt=60,xpd=T);par(srt=0)



ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])


legend(par('usr')[2]*0.98,par('usr')[4]*1,legend=rownames(data),col=col,pch=15,bty="n",cex=1)
par(xpd=T)
type <- ifelse(colnames(data) %in% Low ,"Low","High")
boxcolor <- ifelse(type=="Low","#00BFC4","#F8766D")
plot()
points(a1, rep(par('usr')[4]*1.02,length(a1)), pch = 15, col = boxcolor,lwd=1, cex=1.0)
legend(par('usr')[2]*0.45,par('usr')[4]*1.1,legend=c("Low","High"),col=c("#00BFC4","#F8766D"),pch=15,bty="n",cex=1,ncol=2)

dev.off()







library(ggcorrplot)
data= as.data.frame((immune))
corr <- cor(data)
ggcorrplot(corr,lab = TRUE)
ggsave(filename="07CIBERSORT/immune_cell_cor.pdf",width = 13,height = 10)

data_group= data_group[,c("id","risk")]
library(reshape2)
library(ggpubr)
# 加载用于处理数据格式的reshape2包
cellType=colnames(data)
# 从data矩阵中提取物种分类信息
data_frame=data.frame(t(data), cellType)
# 新建数据框
colnames(data_frame)=gsub("\\.","-",colnames(data_frame))

data_frame=melt(data_frame, id='cellType')  
names(data_frame)[2]='id'
# 重命名variable为sample_id，保持与data_group的样品变量名一致
data_frame=merge(data_frame, data_group, by='id')
# 根据样品变量名，给data_frame添加分组信息，如下：
p1 <- ggboxplot(data_frame, x = "risk", y = "value",
               fill  = "risk", palette = c("#48b34e", "#f91f10"),
               add = c("boxplot"),add.params = list(fill="white"),
               facet.by = "cellType")+stat_compare_means(label = "p.format") +xlab("")


ggsave(p1, filename="07CIBERSORT/immune_cell.pdf",width = 12,height = 10)









immune=read.table("07CIBERSORT/CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

# rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(immune))
immune=avereps(immune)


 
# tmb=read.table("07ETV7/sampleType.csv",sep=",",header=T,check.names=F,row.names=1)
# tmb=as.matrix(tmb)
# row.names(tmb)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(tmb))
# tmb=avereps(tmb)
# lowTmb=tmb[tmb[,"TMB"]<= median(tmb[,"TMB"]),]
# highTmb=tmb[tmb[,"TMB"]>median(tmb[,"TMB"]),]
# lowTmbName=names(lowTmb)
# highTmbName=names(highTmb)
# 根据表达量分为高低表达组
# 高低风险样本
risk <- read.table("05Cox/trainSet_riskScore.xls",sep="\t",header=T,check.names=F)
# risk$id=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",risk$id)
Normal= risk$id[risk$risk=="Low Risk"]
Tumor=  risk$id[risk$risk == "High Risk"] #肿瘤组
# # 
# # # tmb=avereps(tmb)
# # 
# # immuneName <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\",rownames(immune) )
immuneName <- rownames(immune)
lowTmbImm <-  rownames(immune[immuneName %in% Normal,])
highTmbImm <- rownames(immune[immuneName %in% Tumor,])

# # lowTmbImm=intersect(row.names(immune),lowTmbName)
# # highTmbImm=intersect(row.names(immune),highTmbName)
rt=rbind(immune[lowTmbImm,],immune[highTmbImm,])
rt[1:3,1:5]
rt1 <- cbind(id=rownames(rt),rt)
write.csv(rt1,file = "07ETV7/clinical_heatmap.csv",quote = F,row.names = F)

lowTmbNum=length(lowTmbImm)
highTmbNum=length(highTmbImm)
# > highTmbNum
# [1] 59
# > lowTmbNum
# [1] 116

pdf("07CIBERSORT/immuneCell_vioplot.pdf",height=8,width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")


for(i in 1:ncol(rt)){
  if(sd(rt[1:lowTmbNum,i])==0){
    rt[1,i]=0.001
  }
  if(sd(rt[(lowTmbNum+1):(lowTmbNum+highTmbNum),i])==0){
    rt[(lowTmbNum+1),i]=0.001
  }
  lowTmbData=rt[1:lowTmbNum,i]
  highTmbData=rt[(lowTmbNum+1):(lowTmbNum+highTmbNum),i]
  vioplot(lowTmbData,at=3*(i-1),lty=1,add = T,col = "#004BFB")
  vioplot(highTmbData,at=3*(i-1)+1,lty=1,add = T,col = "#F91F10")
  wilcoxTest=wilcox.test(lowTmbData,highTmbData)
  p=wilcoxTest$p.value
  mx=max(c(lowTmbData,highTmbData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
par(srt=60,xpd=T);par(srt=0)
legend(par('usr')[2]*0.9,par('usr')[4]*1,legend=c("Low","High"),col=c("#004BFB","#F91F10"),pch=15,bty="n",cex=1)
dev.off()




### 画差异基因和免疫细胞的相关性热图
immune=read.table("07CIBERSORT/CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<0.05,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

### 读入表达量数据
tcga_expr <- read.table("07CIBERSORT/GSE70866_removeBatchEffect_expres.txt",sep="\t", row.names = 1,header = T )
tcga_expr <- as.data.frame(t(tcga_expr))

tcga_expr <- tcga_expr[,rownames(immune)]
tcga_expr[1:5,1:6]
# 批量计算相关性
genelist <- rownames(tcga_expr)
genelist <- c("IL1R2","S100A12","CCL8")
immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(immune)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(immune[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
immuscore("IL1R2")
# 批量计算genelist跟免疫浸润相关性的结果
data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
data <- na.omit(data)
#保存到文件
write.table(data, "07CIBERSORT/immune_cells_gene_correlation.xls", quote = F, row.names = F,sep = "\t")
# 增加一列，区分p值的大小
# 使用两个ifelse实现三分类
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]

ggplot(data, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 8))+#调整y轴文字
  #调整legen
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
ggsave("07CIBERSORT/immune_cells_gene_cor.pdf",width = 12,height = 4)
