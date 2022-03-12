# 1.加载R包获取数据。
rm(list=ls())
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")
## 加载R包
library(GEOquery)
library(Biobase)
library(limma)

GEO_id <- "GSE70866"
## 加载R包
library(GEOquery)
## 下载数据，如果文件夹中有会直接读入.以后只要更改GSE号，就可以直接下载别的GEO数据
gset = getGEO(GEO_id, destdir="01rawData/",getGPL = F)
## 获取ExpressionSet对象，包括的表达矩阵和分组信息
gset=gset[[1]]


# 2.通过pData函数获取分组信息
####################################
pdata=pData(gset)
# ## 只要后36个,本次选择的这36个是配对的。
# ## 所以，别人的芯片我们也不是全要，我们只要适合自己的数据
# group_list=c(rep('before',18),rep('after',18))
# group_list=factor(group_list)
# ## 强制限定顺序
# group_list <- relevel(group_list, ref="before")
pd1=pdata[,apply(pdata, 2, function(x){
  length(unique(x))>2})]  #缩小范围
dim(pd1)
# 数据框再用apply循环去查找文章作者是用哪一列来分组的
apply(pd1,2,table)

# write.table(pdata,file =paste0(GEO_id, "GSE70866-GPL14550_clinical.tsv"),sep = "\t",quote = F)
write.table(pdata,file =paste0("01rawData/GSE70866-GPL14550_clinical.tsv"),sep = "\t",quote = F)

clincal_raw <- read.delim("01rawData/GSE70866-GPL14550_clinical.tsv",check.names = F,header = T,sep = "\t")
##获取原始数据
##用getGEOSUppFiles()函数获取原始数据#####
##rawdata = getGEOSuppFiles(“GSE20706”)
##已经下载好了原始数据并存放在GSE20706文件夹
##数据解压到新建的文件夹中
# dir.create("sampFile")##创建新的文件夹
samPath = "./01rawData/GSE70866_RAW/"
files <- list.files(samPath,full.names = TRUE)##显示文件夹中的文件，及原始文件顺序
# unzip(zipfile = files[2],exdir = samPath)##解压原始文件到sampath文件夹中,文件在第几个，[]中就是几

##获取原始数据，以下代码可行，并且删除了不需要的样本的txt格式文件
dat = list.files(samPath, pattern = "txt.gz")
dat
dat = dat[substr(dat,1,10) %in% clincal_raw$id]
# 9.28日,还剩去掉没有生存时间的样本，还剩余74个样本
dat=read.maimages(files= dat,path = samPath,green.only = TRUE,source="agilent")
#需要查单通道还是两通道

ep=dat$E
# ep=log2(ep)##如果没有进行log2处理的，需要进行log2数据处理
# 判断是否需要进行数据转换
ex <- ep
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
ep <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

par(cex = 0.9)
# rainbow是R的一个函数，用于产生彩虹色
cols <- rainbow(ncol(ep) * 1.2)
boxplot(ep, col = cols, xlab = "Sample", ylab = "Log intensity")

##前期处理
####背景校正和标准化处理####
RG = backgroundCorrect(dat, "auto", normexp.method = "rma", offset=50) 
E = normalizeBetweenArrays(RG, method="quantile")
##查看标准化前后
plotMA(dat,array = 1)
plotMA(RG,array = 1)
plotMA(E,array = 1)

plotDensities(dat)
plotDensities(RG)
plotDensities(E)
##绘制标准化后数据box分布图，查看表达矩阵的整体分布情况
ep=as.data.frame( E )
ep=ep[,6:length(colnames(ep))]
par(cex = 0.9)
# rainbow是R的一个函数，用于产生彩虹色
cols <- rainbow(ncol(ep) * 1.2)
boxplot(ep, col = cols, xlab = "Sample", ylab = "Log intensity")

expr <- as.data.frame( E )##得到标准化后的表达数据
expr$ID <- E[["genes"]]$ProbeName

##下载并提取GPL平台数据
# gpl <- getGEO("GPL14550", destdir="01rawData/") ##根据GPL号下载的是芯片设计的信息
# GPL=gpl@dataTable@table
GPL <- read.table('01rawData/GPL14550-9757.txt',header = TRUE,sep = '\t',quote = '', comment.char = '#',check.names = F,fill = T)
colnames(GPL)

library(dplyr)
colnames(GPL)
##得到探针-基因ID-genesymbol的注释文件
annotation <- as.data.frame(GPL) %>%
  dplyr::select("ID","GENE","GENE_SYMBOL") %>% ##跟进GPL的情况选择
  dplyr::filter(GENE_SYMBOL != "")
#colnames(annotation)[1]="ID"##根据需要调整

expr2 <- merge(annotation,expr, by='ID')##根据探针信息对表达数据进行genesymbol注释
##按geneID合并，相同者取平均值
exprSet2 <- aggregate(x = expr2[,9:ncol(expr2)],
                      by = list(expr2$GENE_SYMBOL),
                      FUN = mean)

exprSet6=data.frame(annotation[match(exprSet2$Group.1,annotation$GENE_SYMBOL),],exprSet2)  
exprSet6=exprSet6[,-c(1,2,4)]
exprSet6[1:3,1:6]
##对样本名称进行处理，仅保留GSM号
m=colnames(exprSet6)[2:length(colnames(exprSet6))]
m= substr(m, 1, 10)
colnames(exprSet6)[2:length(colnames(exprSet6))]=m
head(exprSet6)
write.table(exprSet6,file = "01rawData/GSE70866-GPL14550.xls",quote = F,row.names = F,sep = "\t")

save(exprSet6,file = "01rawData/GSE70866-GPL14550.Rdata")


if(F){
## 下载数据，如果文件夹中有会直接读入.以后只要更改GSE号，就可以直接下载别的GEO数据
# GEO_id <- "GSE20706"
# gset = getGEO(GEO_id, destdir="01rawData/",getGPL = F)
# ## 获取ExpressionSet对象，包括的表达矩阵和分组信息
# names(gset)
# gset=gset[[1]]
# exprSet1 <- exprs(gset)
# 
# # 通过pData函数获取分组信息
# ####################################
# pdata=pData(gset)


# 将所有样本的raw data文件放置在R工作目录下（Microarray文件夹），然后用list.celfiles函数查看所有文件，
# 再用read.celfiles函数读入R中，并将上述的pdata信息也一并写入
library(oligo)

file_CELs <- list.celfiles("01rawData/GSE27957_RAW/", listGzipped = TRUE, full.name = TRUE,pattern = "CEL.gz")
# rawAffy <- read.celfiles(filenames = file_CELs, phenoData = phenoData(gset), sampleNames = rownames(pData(gset)))
file_CELs <-file_CELs[ substr(file_CELs,24,32) %in% clincal_raw$id]#  45个样本
rawAffy <- read.celfiles(filenames = file_CELs)

# 粗略看下，rawAffy数据中包含了芯片上所有探针的强度值信息
head(exprs(rawAffy), n = 2)

# 芯片常用的标准化方法也有好几种（如：MAS5），这里就不一一尝试了，就拿最常用的RMA（Robust Multichip Average algorithm），
# 引用芯片教程，RMA标准化过程主要分为3步：
eset <- rma(rawAffy)
exprSet <- exprs(eset)
dim(exprSet)
##[1] 45101     6
par(cex = 0.9)
# rainbow是R的一个函数，用于产生彩虹色
cols <- rainbow(ncol(exprSet) * 1.2)
boxplot(exprSet, col = cols, xlab = "Sample", ylab = "Log intensity")
plotDensities(exprSet)

ex = normalizeBetweenArrays(exprSet, method="quantile")
# # 判断是否需要进行数据转换
# # ex <- ep
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0) ||
#   (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# 
# if (LogC) { ex[which(ex <= 0)] <- NaN
# exprSet1 <- log2(ex)
# print("log2 transform finished")}else{print("log2 transform not needed")}


exprSet1 <- ex
par(cex = 0.9)
# rainbow是R的一个函数，用于产生彩虹色
cols <- rainbow(ncol(exprSet1) * 1.2)
boxplot(exprSet1, col = cols, xlab = "Sample", ylab = "Log intensity")
plotDensities(exprSet1)

##下载并提取GPL平台数据
# gpl <- getGEO("GPL5175", destdir="01rawData/") ##根据GPL号下载的是芯片设计的信息
# GPL=gpl@dataTable@table
library(devtools)
install_github("jmzeng1314/idmap3")
library(idmap2)
library(stringr)
ids=idmap3::get_pipe_IDs('GPL11532') 
ids=get_soft_IDs('GPL11532') 
# gpl <- read.table('01rawData/GPL1',header = TRUE,sep = '\t',quote = '', comment.char = '#',check.names = F,fill = T)
colnames(gpl)
GPL <- gpl %>%  dplyr::select("ID","GB_LIST","gene_assignment")%>%
  tidyr::separate_rows(GB_LIST,sep=",") %>% 
  dplyr::filter(GB_LIST != "") 
  

GPL$GENE_SYMBOL <- gsub("(.*?)\\/\\/ (.*?) \\/.*","\\2\\",stringr::str_extract(GPL$gene_assignment, paste0( GPL$GB_LIST," // .*? //" )))

annotation <- GPL %>%  dplyr::select("ID","GENE_SYMBOL")%>%
  dplyr::filter(GENE_SYMBOL != "") %>% unique()


library(dplyr)
# colnames(GPL)
# ##得到探针-基因ID-genesymbol的注释文件
# annotation <- as.data.frame(GPL) %>%
#   dplyr::select("ID","GENE","GENE_SYMBOL") %>% ##跟进GPL的情况选择
#   dplyr::filter(GENE_SYMBOL != "")
# #colnames(annotation)[1]="ID"##根据需要调整
exprSet1[1:3,1:5]
expr <- as.data.frame(exprSet1) 
expr$ID <- rownames(expr)
expr2 <- merge(annotation,expr, by='ID')##根据探针信息对表达数据进行genesymbol注释
##按geneID合并，相同者取平均值
exprSet2 <- aggregate(x = expr2[,2:ncol(expr2)],
                      by = list(expr2$GENE_SYMBOL),
                      FUN = mean)

exprSet6=data.frame(annotation[match(exprSet2$Group.1,annotation$GENE_SYMBOL),],exprSet2)  
exprSet6=exprSet6[,-c(1,3,4)]
exprSet6[1:3,1:6]
##对样本名称进行处理，仅保留GSM号
m=colnames(exprSet6)[2:length(colnames(exprSet6))]
m= substr(m, 1, 9)
colnames(exprSet6)[2:length(colnames(exprSet6))]=m
head(exprSet6)
write.table(exprSet6,file = "01rawData/GSE27957_GPL5175.xls",quote = F,row.names = F,sep = "\t")

save(exprSet6,file = "01rawData/GSE27957_GPL5175.Rdata")
}

clincal_raw <- read.delim("01rawData/GSE70866-GPL17077_clinical.tsv",check.names = F,header = T,sep = "\t")
##获取原始数据
##用getGEOSUppFiles()函数获取原始数据#####
##rawdata = getGEOSuppFiles(“GSE20706”)
##已经下载好了原始数据并存放在GSE20706文件夹
##数据解压到新建的文件夹中
# dir.create("sampFile")##创建新的文件夹
samPath = "./01rawData/GSE70866_RAW/"
files <- list.files(samPath,full.names = TRUE)##显示文件夹中的文件，及原始文件顺序
# unzip(zipfile = files[2],exdir = samPath)##解压原始文件到sampath文件夹中,文件在第几个，[]中就是几

##获取原始数据，以下代码可行，并且删除了不需要的样本的txt格式文件
dat = list.files(samPath, pattern = "txt.gz")
dat
dat = dat[substr(dat,1,10) %in% clincal_raw$id]
# 9.28日,还剩去掉没有生存时间的样本，还剩余74个样本
dat=read.maimages(files= dat,path = samPath,green.only = TRUE,source="agilent")
#需要查单通道还是两通道

ep=dat$E
# ep=log2(ep)##如果没有进行log2处理的，需要进行log2数据处理
# 判断是否需要进行数据转换
ex <- ep
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
ep <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

par(cex = 0.9)
# rainbow是R的一个函数，用于产生彩虹色
cols <- rainbow(ncol(ep) * 1.2)
boxplot(ep, col = cols, xlab = "Sample", ylab = "Log intensity")

##前期处理
####背景校正和标准化处理####
RG = backgroundCorrect(dat, "auto", normexp.method = "rma", offset=50) 
E = normalizeBetweenArrays(RG, method="quantile")
##查看标准化前后
plotMA(dat,array = 1)
plotMA(RG,array = 1)
plotMA(E,array = 1)

plotDensities(dat)
plotDensities(RG)
plotDensities(E)
##绘制标准化后数据box分布图，查看表达矩阵的整体分布情况
ep=as.data.frame( E )
ep=ep[,6:length(colnames(ep))]
par(cex = 0.9)
# rainbow是R的一个函数，用于产生彩虹色
cols <- rainbow(ncol(ep) * 1.2)
boxplot(ep, col = cols, xlab = "Sample", ylab = "Log intensity")

expr <- as.data.frame( E )##得到标准化后的表达数据
expr$ID <- E[["genes"]]$ProbeName

##下载并提取GPL平台数据
# gpl <- getGEO("GPL14550", destdir="01rawData/") ##根据GPL号下载的是芯片设计的信息
# GPL=gpl@dataTable@table
GPL <- read.table('01rawData/GPL17077-17467.txt',header = TRUE,sep = '\t',quote = '', comment.char = '#',check.names = F,fill = T)
colnames(GPL)

library(dplyr)
colnames(GPL)
##得到探针-基因ID-genesymbol的注释文件
annotation <- as.data.frame(GPL) %>%
  dplyr::select("ID","GB_ACC","GENE_SYMBOL") %>% ##跟进GPL的情况选择
  dplyr::filter(GENE_SYMBOL != "")
#colnames(annotation)[1]="ID"##根据需要调整

expr2 <- merge(annotation,expr, by='ID')##根据探针信息对表达数据进行genesymbol注释
##按geneID合并，相同者取平均值
exprSet2 <- aggregate(x = expr2[,9:ncol(expr2)],
                      by = list(expr2$GENE_SYMBOL),
                      FUN = mean)

exprSet6=data.frame(annotation[match(exprSet2$Group.1,annotation$GENE_SYMBOL),],exprSet2)  
exprSet6=exprSet6[,-c(1,2,4)]
exprSet6[1:3,1:6]
##对样本名称进行处理，仅保留GSM号
m=colnames(exprSet6)[2:length(colnames(exprSet6))]
m= substr(m, 1, 10)
colnames(exprSet6)[2:length(colnames(exprSet6))]=m
head(exprSet6)
write.table(exprSet6,file = "01rawData/GSE70866-GPL17077.xls",quote = F,row.names = F,sep = "\t")

save(exprSet6,file = "01rawData/GSE70866-GPL17077.Rdata")

###########完结，下一步做批次效应

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)

pca_plot = function(dddd,ggggg){
  library("FactoMineR")
  library("factoextra")
  df.pca <- PCA(t(dddd), graph = FALSE)
  fviz_pca_ind(df.pca,
               #axes = c(2,3),
               geom.ind = "point",
               col.ind = ggggg ,
               addEllipses = TRUE,
               legend.title = "Groups"
  )
}
# 下面的 step1-output.Rdata 文件，大家可以去学习我的GEO课程
# 就知道如何制作啦。
load(file = '01rawData/GSE70866-GPL14550.Rdata')
dat1 <- exprSet1

load(file = '01rawData/GSE70866-GPL17077.Rdata')
dat2 <- exprSet1

  
  
dat <- merge(dat1,dat2,by="symbol")
group_list <- data.frame( samples = c(colnames(dat1)[-1],colnames(dat2)[-1]),type= rep(c("GSE70866-GPL14550","GSE70866-GPL17077"),c(length(colnames(dat1)[-1]), length(colnames(dat2)[-1]))) )


rownames(group_list) <-group_list$samples
group_list <- group_list[,-1]
# 每次都要检测数据
rownames(dat) <- dat[,1]
dat <- dat[,-1]
dat[1:4,1:4]

g=factor( group_list )
g
pdf(file="01rawData/pca_raw.pdf",height=5,width=6.5)
pca_plot(dat,g)
dev.off()

# sva去除批次效应
## 使用 sva 的 ComBat 函数
group_list <- data.frame( samples = c(colnames(dat1)[-1],colnames(dat2)[-1]),type= rep(c(1,2,3),c(62,50,64) ))
library(sva)
batch <- group_list

ex_b_sva = ComBat(dat=as.matrix(dat), 
                  batch=batch$type 
)
ex_b_sva[1:4,1:4]



pdf(file="01rawData/pca_removeBatchEffect.pdf",height=5,width=6.5)
pca_plot(ex_b_sva,g)
dev.off()

GSE28221<- cbind(id = rownames(ex_b_sva),ex_b_sva)
GSE28221[1:4,1:4]
write.table(GSE28221,file = "01rawData/GSE70866_removeBatchEffect_expres.xls",quote = F,row.names = F,sep = "\t")

###############################    PCA  umap
library(umap)
library(Rcpp)
library(factoextra)
dat <- read.table("01rawData/GSE70866_removeBatchEffect_expres.xls",head = TRUE,row.names=1,sep="\t")
dat[1:4,1:4]
gs <- read.table("01rawData/geneset.txt",sep = "\t",check.names = F,header = TRUE)
exprSet_FerrDb <- dat[rownames(dat) %in%  gs$HALLMARK_HYPOXIA,]
pheatmap(exprSet_FerrDb,scale = "row")


exprSet_FerrDb <- as.matrix(exprSet_FerrDb)
#一步完成聚类
library(ConsensusClusterPlus)
title="./02ConsensusClusterPlus"
results = ConsensusClusterPlus(exprSet_FerrDb,maxK=6,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")




cl_i <- kmeans(t(exprSet_FerrDb), centers = 2)
data1 <- data.frame(cl_i$cluster)
data1 <- cbind(id = rownames(data1) ,data1)
write.table(data1,file = "01rawData/group.xls",quote = F,row.names = F,sep = "\t")
# groups <- read.table("group1.txt",header = TRUE,sep="\t")
k
cl_i <- kmeans(t(dat), centers = 2)

result <- dist(t(dat),method = "euclidean")
result_hc <- hclust(d=result,method = "complete")
fviz_dend(result_hc,k=2,cex = 0.6)


colnames(data1) <- c('Coordinate1','Coordinate2','Cluster')
# data1$Cluster <- ifelse(data1$Cluster==1,'Cluster2','Cluster1')
data1$Cluster <- ifelse(data1$Cluster==1,'Cluster1','Cluster2')



iris.umap = umap(t(exprSet_FerrDb))

plot(iris.umap$layout,col=cl_i$cluster,pch=16,asp = 1,
     xlab = "UMAP_1",ylab = "UMAP_2")
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")
# 添加图例
legend("topright",title = "Species",inset = 0.01,
       legend = unique(groups$group),pch=16,
       col = unique(groups$group))

tsne<- Rtsne(dat[,1:4], dims = 2, perplexity=i, verbose=TRUE, max_iter = 500)

#二次聚类，因为tsne仅仅是更好的展示聚类的结果
set.seed(i)
cl_i <- kmeans(tsne$Y, centers = 2, nstart =10,algorithm="Lloyd")
data1 <- data.frame(tsne$Y,cl_i$cluster)
colnames(data1) <- c('Coordinate1','Coordinate2','Cluster')
# data1$Cluster <- ifelse(data1$Cluster==1,'Cluster2','Cluster1')
data1$Cluster <- ifelse(data1$Cluster==1,'Cluster1','Cluster2')





fit=lmFit(ex_b_sva,design)  
fit=eBayes(fit) 
options(digits = 4) 
topTable(fit,coef=2,adjust='BH') 
# 首先是瘾君子与正常人的差异分析
deg3=topTable(fit,coef=2,adjust='BH',number = Inf)
pca_plot(ex_b_sva,g)



table(group_list)
library(limma)
g=factor( group_list )
g
g=relevel(g,'con')
design=model.matrix(~g) 
fit=lmFit(dat,design) 
fit=eBayes(fit) 
options(digits = 4) 
topTable(fit,coef=2,adjust='BH') 
# 首先是瘾君子与正常人的差异分析
deg1=topTable(fit,coef=2,adjust='BH',number = Inf)
pca_plot(dat,g)



pd1=pdata[,apply(pdata, 2, function(x){
  length(unique(x))>2})]  #缩小范围
dim(pd1)
# 数据框再用apply循环去查找文章作者是用哪一列来分组的
apply(pd1,2,table)

write.table(pdata,file =paste0(GEO_id, "_clinical.tsv"),sep = "\t",quote = F)

Normal=rownames(pd1[grepl("normal ovarian surface epithelium \\(OSE\\)",as.character(pdata$source_name_ch1)),]) #正常组
Tumor=rownames(pd1[grepl( "papillary serous ovarian adenocarcinoma",as.character(pdata$source_name_ch1)),])#肿瘤组

group_list=c(rep('Normal',length(Normal)),
             rep('Tumor',length(Tumor))) #分组信息
group_list=factor(group_list)
## 强制限定顺序
group_list <- relevel(group_list, ref="Normal")
table(group_list)
# group_list=c(rep('Tumor',length(Tumor))) #分组信息
# group_list=factor(group_list)
## 强制限定顺序
table(group_list)

# 3.通过exprs函数获取表达矩阵并校正
####################################
exprSet=exprs(gset)
exprSet=exprSet[,c(Normal,Tumor)] #对表达矩阵取子集
# exprSet=exprSet[,c(Tumor)] #对表达矩阵取子集
#先简单看一下整体样本的表达情况
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)

# 需要人工校正一下，用的方法类似于Quntile Normalization
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
# boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)

######################################
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

# 5.获取注释信息
####################################
# platformMap <- data.table::fread("D:/BaiduNetdiskDownload/myfirstGEO/platformMap.txt")
# index = gset@annotation
# platformDB = paste0(platformMap$bioc_package[grep(index,platformMap$gpl)],".db")
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# # 
# # BiocManager::install("illuminaHumanv4.db")
# ## 加载R包 
# library("illuminaHumanv4.db")
# ## 获取探针对应的symbol信息
# ## 获取表达矩阵的行名，就是探针名称
# probeset <- rownames(exprSet)
# ## 使用lookup函数，找到探针在illuminaHumanv2.db中的对应基因名称
# ## 如果分析别的芯片数据，把illuminaHumanv2.db更换即可
# SYMBOL <-  annotate::lookUp(probeset,"illuminaHumanv4.db", "SYMBOL")
# ## 转换为向量
# symbol = as.vector(unlist(SYMBOL))
# # 制作probe2symbol转换文件
# probe2symbol = data.frame(probeset,symbol,stringsAsFactors = F)
# ##############################
# # 6.探针转换与基因去重
##############################
# 一个基因会对应对个探针，有些基因名称会是重复的，这些都需要处理。
# 对于，多个探针，我们选取在样本中平均表达量最高的探针作为对应基因的表达量。一下代码完成所有事情，而且可以复用。


library(dplyr)
library(tibble)
# anno <-data.table::fread("GPL/GPL5474_family.soft",skip ="ID")
# probe2symbol <- anno %>%
# 
#   dplyr::select("ID","Symbol") %>% dplyr::rename(probe_id = "ID",symbol="Symbol") %>%
# 
#   filter(symbol != "") %>%
# 
#   separate(gene_assignment,c("drop","symbol"),sep="//") %>%
# 
#   select(-drop)
gpl <- read.table('GPL/GPL570-55999.txt',header = TRUE,sep = '\t',quote = '', comment.char = '#',check.names = F,fill = T)
colnames(gpl)

probe2id <- gpl %>%

  dplyr::select('ID', "Gene Symbol") %>%
  # dplyr::select('ID','Symbol') %>%
  tidyr::separate_rows("Gene Symbol",sep = ' /// ') %>%
  dplyr::rename(probe_id = 'ID',symbol="Gene Symbol") %>%

  dplyr::filter(symbol != '',symbol != '')

probe2id$probe_id <- as.character(probe2id$probe_id)


exprSet1 <- exprSet   %>% as.data.frame() %>%
  tibble::rownames_to_column(var="probe_id") %>% 
  as_tibble() %>% 
  #合并探针的信息
  dplyr::inner_join(probe2id,by="probe_id") %>% 
  #dplyr::filter(trans_biotype == "protein_coding") %>% 
  #去掉多余信息
  dplyr::select(-probe_id) %>% 
  #重新排列
  dplyr::select(symbol,everything()) %>% 
  #去除symbol中的NA
  dplyr::filter(symbol != "NA") %>% 
  #以symbol分组就平均数
  dplyr::group_by(symbol) %>%
  dplyr::summarise_all(mean) %>%
  # 列名变成行名
  tibble::column_to_rownames(var = "symbol")

exprSet_out <- cbind(id=rownames(exprSet1),exprSet1)
exprSet_out[1:4,1:5]
write.table(exprSet_out, file= paste0(GEO_id, "_symbol_exprSet.txt"),row.names = F,quote = F,sep = "\t")

# exprSet_df <- read.table(paste0(GEO_id, "_symbol_exprSet.txt"),header = T,row.names = 1,sep = "\t")
# exprSet_df <- as.data.frame(t(exprSet_df))
# exprSet_df[1:4,1:5]
# exprSet_df$id <- rownames(exprSet_df)
# 
# clinical <- read.table(paste0(GEO_id,"_clinical_clean.txt"),header = T,sep = "\t")
# clincial_exprSet <- merge(clinical,exprSet_df,by="id")
# clincial_exprSet[1:3,1:4]
# write.table(clincial_exprSet,paste0(GEO_id, "clinical_exprSet.txt"),sep = "\t",row.names = F)
#############################################
# 7.差异分析
#############################################
# 我们选取格式比较简单的。如果没有配对信息，差异分析这样做：
design=model.matrix(~ group_list)
fit=lmFit(exprSet1,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef="group_listTumor",number=Inf,p.value=0.05)
write.csv(allDiff,file=paste0(GEO_id, "_diff.txt"),row.names = T)

# 我们这里实际上是有配对信息的，差异分析应该这样做
# pairinfo = factor(rep(1:18,2))
# design=model.matrix(~ pairinfo+group_list)
# fit=lmFit(exprSet,design)
# fit=eBayes(fit) 
# allDiff_pair=topTable(fit,adjust='BH',coef="group_listafter",number=Inf,p.value=0.05)

#########################################################
# 8.作图验证(非必要)
#######################################################
# 得到ggplot2喜欢的数据格式，行是观测，列是变量
# 首先基因名称需要在列，所以需要用t函数，行列转置。
data_plot = as.data.frame(t(exprSet1))
# data_plot = data.frame(pairinfo=rep(1:18,2),
#                        group=group_list,
#                        data_plot,stringsAsFactors = F)

data_plot = data.frame(
                       group=group_list,
                       data_plot,stringsAsFactors = F)

data_plot[1:3,1:5]
write.table(data_plot,file = paste0(GEO_id, "_exprSet_type.txt"),row.names = T,quote = F,sep = "\t")


allDiff["ETV7",]
allDiff["S1PR2",]
allDiff["S1PR3",]
allDiff["S1PR4",]
allDiff["S1PR5",]
allDiff["SPHK1",]
allDiff["SPHK2",]
allDiff["ALOX15",]
allDiff["FGF21",]
allDiff["MBTPS1",]
c('ETV7') %in% rownames(allDiff)

genes <- c('ALOX15','FGF21','S1PR3','S1PR4')
# 挑选了一个基因CAMKK2

data <- read.table("F:\\November\\YQB049\\YQB049数据补充\\GSE133624_reads-count-all-sample.txt\\GSE133624_reads-count-all-sample.txt",header = T,sep = "\t",row.names = 1)
data <- log2( data["ENSG00000010030",])
data_plot <- t(data)
data_plot = data.frame(
  group= rep(c("Normal","Tumor"),c(29,36)) ,
  data_plot,stringsAsFactors = F)
colnames(data_plot)[2] <- "ETV7"

cliTest<-wilcox.test(ETV7 ~ group, data=data_plot,exact=FALSE)
pValue=cliTest$p.value
pval=ifelse(pValue<0.001,"<0.001",paste0("=",sprintf("%.03f",pValue)))
b = boxplot(ETV7 ~ group, data = data_plot,outline = T, plot=F)
data_plot <- data_plot[(!data_plot$ETV7 %in% b$out),]
yMin=min(b$stats);yMax = max(b$stats/15+b$stats)
ySeg = max(b$stats/20+b$stats);ySeg2 = max(b$stats/22+b$stats)
n = ncol(b$stats)
tab1=table(data_plot$group)

pdf(file="GSE7476_ETV7_expression.pdf",width = 5,height = 4.5)
par(mar = c(4.5,6,3,3))
# 	boxplot(expression ~ clinical, data = data,names=xlabel,col=dotCol,
# 		ylab = "TMB",main=clinical,xlab="",
# 		cex.main=1.3, cex.lab=1.2, cex.axis=1.1,ylim=c(yMin,yMax),outline = FALSE)
#     segments(1,ySeg, n,ySeg);segments(1,ySeg, 1,ySeg2);segments(n,ySeg, n,ySeg2)
#     text((1+n)/2,ySeg,labels=paste0("p",pval),cex=1.2,pos=3)
res <- boxplot(ETV7 ~ group, data = data_plot,outline = F,
               ylab = "ETV7 Expression",xlab="" ,boxwex=0.5,col="white",border=c("#004BFB", "#F91F10"),
               cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax))


# points(jitter(rep(1:2, c(as.vector(tab1)[1],as.vector(tab1)[2])), 1)[1:as.vector(tab1)[1]],unlist(split(data_plot$ETV7, data_plot$group))[1:as.vector(tab1)[1]],col=c("#004BFB"),cex=0.5, pch=16)
# points(jitter(rep(1:2, c(as.vector(tab1)[1],as.vector(tab1)[2])), 1)[as.vector(tab1)[1]+1:as.vector(tab1)[2]],unlist(split(data_plot$ETV7, data_plot$group))[as.vector(tab1)[1]+1:as.vector(tab1)[2]],col=c("#F91F10"),cex=0.5, pch=16)
segments(1,ySeg, n,ySeg);
segments(1,ySeg, 1,ySeg2)
segments(n,ySeg, n,ySeg2)
text((1+n)/2,ySeg,labels=paste("p",pval,sep=""),cex=1.5,pos=3)
title("GSE7476")
dev.off()


getwd()



x = c(7.82, 8.24, 8.21,8.45,8.37,7.69, 8.35,8.47, 9.07)
y = c(7.39,7.05,7.23,5.55,5.41, 6.56, 6.46,6.18)
wilcox.test(x,y)




library(ggplot2)
p1 <-ggplot(data_plot, aes(group,ETV7,fill=group)) +
  geom_boxplot(aes(colour=group), size=2, alpha=0.5)+#geom_point(size=2, alpha=0.5)+
  xlab("") +
  ylab(paste("Expression of ","S1PR1"))+
  annotate(geom = "text",x=2,y=10,label="P<1.233e-05 log2FC=0.259") +
  theme(plot.title = element_text(hjust = 0.5),axis.title.x  = element_blank(),axis.text.x = element_text(face = "bold",size = 9),axis.text.y = element_text(size = 10))+
  theme_classic()


p2 <-ggplot(data_plot, aes(group,S1PR4,fill=group)) +
  geom_boxplot(aes(colour=group), size=2, alpha=0.5)+#geom_point(size=2, alpha=0.5)+
  xlab("") +
  ylab(paste("Expression of ","S1PR4"))+
  annotate(geom = "text",x=2,y=10.5,label="P<3.918e-05 log2FC=0.525") +
  theme(plot.title = element_text(hjust = 0.5),axis.title.x  = element_blank(),axis.text.x = element_text(face = "bold",size = 9),axis.text.y = element_text(size = 10))+
  theme_classic()

p3 <-ggplot(data_plot, aes(group,SPHK2,fill=group)) +
  geom_boxplot(aes(colour=group), size=2, alpha=0.5)+#geom_point(size=2, alpha=0.5)+
  xlab("") +
  ylab(paste("Expression of ","SPHK2"))+
  annotate(geom = "text",x=2,y=10.1,label="P<0.007 log2FC=-0.161") +
  theme(plot.title = element_text(hjust = 0.5),axis.title.x  = element_blank(),axis.text.x = element_text(face = "bold",size = 9),axis.text.y = element_text(size = 10))+
  theme_classic()

p4 <-ggplot(data_plot, aes(group,MBTPS1,fill=group)) +
  geom_boxplot(aes(colour=group), size=2, alpha=0.5)+#geom_point(size=2, alpha=0.5)+
  xlab("") +
  ylab(paste("Expression of ","MBTPS1"))+
  annotate(geom = "text",x=2,y=10.5,label="0.002 log2FC=0.390") +
  theme(axis.title.x  = element_blank(),axis.text.x = element_text(face = "bold",size = 9),axis.text.y = element_text(size = 10))+
  theme_classic()

library(cowplot)
plot_grid(plotlist=p1,p2,p3,p4, ncol=2,labels = LETTERS[1:4], align = "h")

# 批量地作图
library(dplyr)
library(tibble)
allDiff_arrange <- allDiff %>% 
  rownames_to_column(var="genesymbol") %>% 
  arrange(desc(abs(logFC)))
genes <- allDiff_arrange$genesymbol[1:8]

genes <- c('ALOX15','FGF21','S1PR3','S1PR4')


plotlist <- lapply(genes, function(x){
  data =data.frame(group=data_plot[,"group"],gene=data_plot[,x])
  logFC = allDiff["ALOX15",1]
  P.Value = allDiff["ALOX15",4]
  lablels = paste0()
  ggplot(data, aes(group,gene,fill=group)) +
    geom_boxplot() +
    geom_point(size=2, alpha=0.5) +
    #geom_line(aes(group=pairinfo), colour="black", linetype="11") +
    xlab("") +
    ylab(paste("Expression of ",x))+
    annotate(geom = "text",x=2,y=5,label="P<0.897 log2FC=-0.134") +
    theme(plot.title = element_text(hjust = 0.5),axis.title.x  = element_blank(),axis.text.x = element_text(face = "bold",size = 9),axis.text.y = element_text(size = 10))+
    theme_classic()
    
})

library(cowplot)
plot_grid(plotlist=plotlist, ncol=2,labels = LETTERS[1:4])



