# 1.加载R包获取数据。
rm(list=ls())
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")

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
# apply(pd1,2,table)

# write.table(pdata,file =paste0(GEO_id, "_clinical.tsv"),sep = "\t",quote = F)
clincal_raw <- read.delim("01rawData/GSE70866-GPL14550_clinical.tsv",check.names = F,header = T,sep = "\t")

# 3.通过exprs函数获取表达矩阵并校正
####################################
exprSet=exprs(gset)
exprSet=exprSet[,clincal_raw$id] #对表达矩阵取子集
dim(exprSet)
# exprSet=exprSet[,c(Tumor)] #对表达矩阵取子集
#先简单看一下整体样本的表达情况
boxplot(exprSet,outline=FALSE, notch=T, las=2)


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
gpl <- read.table('01rawData/GPL14550-9757.txt',header = TRUE,sep = '\t',quote = '', comment.char = '#',check.names = F,fill = T)
colnames(gpl)

probe2id <- gpl %>%

  dplyr::select('ID','GENE_SYMBOL') %>%
  # dplyr::select('ID','Symbol') %>%
  tidyr::separate_rows('GENE_SYMBOL',sep = ' /// ') %>%
  dplyr::rename(probe_id = 'ID',symbol='GENE_SYMBOL') %>%

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
  dplyr::summarise_all(mean) #%>%
  # 列名变成行名
  #tibble::column_to_rownames(var = "symbol")

exprSet_out <- cbind(id=rownames(exprSet1),exprSet1)
exprSet_out[1:4,1:5]
write.table(exprSet_out,file = "01rawData/GSE70866-GPL14550.xls",quote = F,row.names = F,sep = "\t")

save(exprSet1,file = "01rawData/GSE70866-GPL14550.Rdata")


# 1.加载R包获取数据。
rm(list=ls())
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")

GEO_id <- "GSE70866"
## 加载R包
library(GEOquery)
## 下载数据，如果文件夹中有会直接读入.以后只要更改GSE号，就可以直接下载别的GEO数据
gset = getGEO(GEO_id, destdir="01rawData/",getGPL = F)
## 获取ExpressionSet对象，包括的表达矩阵和分组信息
gset=gset[[2]]


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
# apply(pd1,2,table)

# write.table(pdata,file =paste0(GEO_id, "_clinical.tsv"),sep = "\t",quote = F)
clincal_raw <- read.delim("01rawData/GSE70866-GPL17077_clinical.tsv",check.names = F,header = T,sep = "\t")

# 3.通过exprs函数获取表达矩阵并校正
####################################
exprSet=exprs(gset)
exprSet=exprSet[,clincal_raw$id] #对表达矩阵取子集
dim(exprSet)
# exprSet=exprSet[,c(Tumor)] #对表达矩阵取子集
#先简单看一下整体样本的表达情况
boxplot(exprSet,outline=FALSE, notch=T, las=2)


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
gpl <- read.table('01rawData/GPL17077-17467.txt',header = TRUE,sep = '\t',quote = '', comment.char = '#',check.names = F,fill = T)
colnames(gpl)

probe2id <- gpl %>%
  
  dplyr::select('ID','GENE_SYMBOL') %>%
  # dplyr::select('ID','Symbol') %>%
  tidyr::separate_rows('GENE_SYMBOL',sep = ' /// ') %>%
  dplyr::rename(probe_id = 'ID',symbol='GENE_SYMBOL') %>%
  
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
  dplyr::summarise_all(mean) #%>%
  # 列名变成行名
  #tibble::column_to_rownames(var = "symbol")

exprSet_out <- cbind(id=rownames(exprSet1),exprSet1)
exprSet_out[1:4,1:5]
write.table(exprSet_out,file = "01rawData/GSE70866-GPL17077.xls",quote = F,row.names = F,sep = "\t")

save(exprSet1,file = "01rawData/GSE70866-GPL17077.Rdata")










exprSet_df <- read.table(paste0(GEO_id, "_symbol_exprSet.txt"),header = T,row.names = 1,sep = "\t")
exprSet_df <- as.data.frame(t(exprSet_df))
exprSet_df[1:4,1:5]
exprSet_df$id <- rownames(exprSet_df)

clinical <- read.table(paste0(GEO_id,"_clinical_clean.txt"),header = T,sep = "\t")
clincial_exprSet <- merge(clinical,exprSet_df,by="id")
clincial_exprSet[1:3,1:4]
write.table(clincial_exprSet,paste0(GEO_id, "clinical_exprSet.txt"),sep = "\t",row.names = F)
#############################################
# 7.差异分析
#############################################
# 我们选取格式比较简单的。如果没有配对信息，差异分析这样做：
design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef="group_listNormal",number=Inf,p.value=0.05)

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
data_plot = as.data.frame(t(exprSet))
# data_plot = data.frame(pairinfo=rep(1:18,2),
#                        group=group_list,
#                        data_plot,stringsAsFactors = F)

data_plot = data.frame(
                       group=group_list,
                       data_plot,stringsAsFactors = F)


write.table(data_plot,file = "PTCL_exprSet.xls",row.names = T,quote = F,sep = "\t")

allDiff["S1PR1",]
allDiff["S1PR2",]
allDiff["S1PR3",]
allDiff["S1PR4",]
allDiff["S1PR5",]
allDiff["SPHK1",]
allDiff["SPHK2",]
allDiff["ALOX15",]
allDiff["FGF21",]
allDiff["MBTPS1",]
c('SPHK1','SPHK2','ALOX15','FGF21','MBTPS1','S1PR1','S1PR2','S1PR3','S1PR4','S1PR5') %in% rownames(allDiff)

genes <- c('ALOX15','FGF21','S1PR3','S1PR4')
# 挑选了一个基因CAMKK2
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



