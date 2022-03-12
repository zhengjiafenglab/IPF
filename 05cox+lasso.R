rm(list=ls())
## 加载R包
# library(GEOquery)
library(limma)
library(dplyr)
library(tibble)

# 单因素Cox分析
library(survival)
pFilter=0.05
rt=read.table("05Cox/GSE70866_intersectGenes_exprset_cox_input.xls",header=T,check.names=F,sep="\t",row.names = 1)

# # rt$`OS-censor`[rt$`OS (months)` >365*5] =0
# # rt$`OS (months)`[rt$`OS (months)` >365*5] =365*5
# # rt$`OS (months)` <- rt$`OS (months)`/365
# set.seed(1751)
# trainId <- sample(rownames(rt),nrow(rt)*0.7,replace=F)
# trainSet <- rt[trainId,]
# validId <- rownames(rt)[ !(rownames(rt) %in% trainId)]
# validSet <- rt[validId,]
# #
# rt <- trainSet
# # 
# a= apply(rt[,-(1:8)],2,function(x) log2(x+1))
# rt =as.data.frame(cbind(rt[,(1:8)],a))
rt1 <- rt

outTab=data.frame()
# sigGenes <- c("futime","fustat","age","gender","pathologic_stage","pathologic_t","pathologic_m","pathologic_n")
# sigGenes=c("EFS(months)","EFS-censor","OS (months)", "OS-censor", "age","gender","efs milestone outcome", "os milestone outcome" )
# sigGenes=c("`OS (months)`", "`OS-censor`")
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  #if(sd(rt[,gene])<1){next}
  #if(grepl("-", gene)){next}
  cox=coxph(Surv(futime,fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    group=ifelse(rt[,gene]>median(rt[,gene]),"high","low")
    diff=survdiff(Surv(futime,fustat) ~group,data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    if(pValue<1){
      sigGenes=c(sigGenes,gene)
      outTab=rbind(outTab,
                   cbind(gene=gene,
                         KM.Pvalue=pValue,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         pvalue=coxP) )
    }
  }
}


# outTab$`Adjusted P` <- p.adjust(as.numeric(outTab$KM_Pvalue), method = "BH")
write.table(outTab,file="05Cox/trainSet_uniCox.xls",sep="\t",row.names=F,quote=F)   
surSigExp=rt1[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="05Cox/trainSet_uniSigExp.xls",sep="\t",row.names=F,quote=F)



#绘制森林图函数
bioForest=function(coxFile=null,forestFile=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <-  gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\", rownames(rt)) 
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width = 6,height = 5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(4,2))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,3,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.8-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.8-0.5*0.2,n+1,'Pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard Ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)+0.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard Ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=1.5)
  abline(v=1,col="black",lty=2,lwd=1.3)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.8)
  axis(1)
  dev.off()
}
############绘制森林图函数############


bioForest(coxFile="05Cox/trainSet_uniCox.xls",forestFile="05Cox/trainSet_uniForest.pdf")




# 多变量cox+LASSO惩罚
library(glmnet)                                         #引用包
rt=read.table("05Cox/trainSet_uniSigExp.xls",header=T,sep="\t",check.names=F,row.names=1)
rt1 <- rt
rt <- rt[rt$futime >0,]

# a= apply(rt[,-(1:4)],2,function(x) log2(x+1))
# rt =as.data.frame(cbind(rt[,(1:4)],a))


x = as.matrix(rt[,-(1:2)])
y = Surv(rt$futime,rt$fustat)
set.seed(8)
cvfit = cv.glmnet(x,
                  y,
                  family = "cox",
                  nfold = 10)
pdf("05Cox/cvfit.pdf",width = 6,height = 6)
plot(cvfit)
# abline(v = c(log(cvfit$lambda.min), log(cvfit$lambda.1se)),lty=2)+
text(x = log(cvfit$lambda.min),y = 10.2,
     paste('Lambda.min\n',round(cvfit$lambda.min,4)),cex=1,adj=0.9)
text(x = log(cvfit$lambda.1se),y = 10.2,
     paste('Lambda.lse\n',round(cvfit$lambda.1se,4)),cex=1,adj=0.9)
dev.off()


fit <- glmnet(x, y, family = "cox",nfold = 10)
pdf(file="05Cox/fit.pdf",width=6,height=5.5)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()


myCoefs <- coef(cvfit, s=cvfit$lambda.min)

# lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0)]
# sigGenes <- c("futime","fustat","age","gender","pathologic_stage","pathologic_t","pathologic_m","pathologic_n")
# sigGenes=c("EFS(months)","EFS-censor","OS (months)", "OS-censor", "age","gender","efs milestone outcome", "os milestone outcome" )
sigGenes=c("futime","fustat")
rt <- rt1[,c(sigGenes,lasso_fea)]
# rt$riskScore <-  apply(rt[,lasso_fea], 1, function(x) {x %*% myCoefs@x})
# trainSet$riskScore <- apply(trainSet[,coxGene], 1, function(x) {x %*% multiCoxSum$coefficients[,"coef"]})
# rt$risk=as.vector(ifelse(rt$riskScore>median(rt$riskScore),"High Risk","Low Risk"))

write.table(cbind(id=rownames(rt),rt),
            file="05Cox/trainSet_lasso_exprset.xls",sep="\t",quote=F,row.names=F)

lassoGene <-  cbind(Gene = lasso_fea,Coef = myCoefs[which(myCoefs != 0 )])
# lassoGene <-  cbind(Gene = lasso_fea,Coef = myCoefs[1:3])
write.table(lassoGene,file="05Cox/lasso_coef.xls",row.names = F,quote = F,sep = "\t")






# # 多因素分析
library(survival)                                         #引用包
rt=read.table("05Cox/trainSet_lasso_exprset.xls",header=T,sep="\t",check.names=F,row.names=1)    #读取输入文件
# # rt=read.table("03veen+Cox/trainSet_uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)    #读取输入文件
# 
# 

# # a= apply(rt[,-(1:4)],2,function(x) log2(x+1))
# # rt =as.data.frame(cbind(rt[,(1:4)],a))
# rt <- rt[,-(9:10)]
rt1 <- rt
# # rt <- rt[,-c(1,2,5:8)] 
#COX模型构建
multiCox=coxph(Surv(futime,fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

Gene=rownames(multiCoxSum$coefficients)[multiCoxSum$coefficients[,"Pr(>|z|)"] < 0.05]
Gene
Gene=gsub("`","",Gene)
multiCox=coxph(Surv(`futime`, `fustat`) ~. , data = rt[,c("futime","fustat",Gene)])
multiCoxSum=summary(multiCox)
multiCoxSum

#输出模型参数
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="05Cox/trainSet_multiCox.xls",sep="\t",row.names=F,quote=F)



#绘制多因素森林图函数
bioForest=function(coxFile=null,forestFile=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <-  gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\", rownames(rt))
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

  #输出图形
  pdf(file=forestFile, width = 6,height = 2.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(4,2))

  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,3,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.8-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.8-0.5*0.2,n+1,'Pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard Ratio',cex=text.cex,font=2,adj=1,)

  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)+0.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard Ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=1.5)
  abline(v=1,col="black",lty=2,lwd=1.3)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.8)
  axis(1)
  dev.off()
}
############绘制森林图函数############

bioForest(coxFile="05Cox/trainSet_multiCox.xls",forestFile="05Cox/trainSet_multiForest.pdf")





#输出病人风险值
# riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
# sigGenes=c("`OS (months)`","`OS-censor`","age","pathologic_stage","pathologic_T","pathologic_M","pathologic_N")
# sigGenes=c("OS (months)", "OS-censor" )
rt1=rt1[,c(sigGenes,coxGene)]

# risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
rt1$riskScore <- apply(rt1[,coxGene], 1, function(x) {x %*% multiCoxSum$coefficients[,"coef"]})
library(survminer)
#对数据集的基因进行bestSeparation统计
res.cut <- surv_cutpoint(rt1, time = "futime",
                         event = "fustat",
                         variables = c("riskScore"),
                         minprop = 0.3) #
# cutpoint statistic
# riskScore 3.342592  6.331819
res.cat <- surv_categorize(res.cut)
rt1$risk=as.vector(ifelse(res.cat$riskScore == "high","High Risk","Low Risk"))

# rt1$risk=as.vector(ifelse(rt1$riskScore>median(rt1$riskScore),"High Risk","Low Risk"))

write.table( cbind(id=rownames(rt1),rt1),
             file="05Cox/trainSet_riskScore.xls",
             sep="\t",row.names = F,
             quote=F)

# sigGenes=c("`OS (months)`","`OS-censor`","age","pathologic_stage","pathologic_T","pathologic_M","pathologic_N")
multiSigExp=rt1[,c(sigGenes,coxGene)]
multiSigExp=cbind(id=row.names(multiSigExp),multiSigExp)
write.table(multiSigExp,file="05Cox/trainSet_multiSigExp.xls",sep="\t",row.names=F,quote=F)

# # 测试集风险分数
# rt <- validSet
# 
# a= apply(rt[,-(1:8)],2,function(x) log2(x+1))
# validSet =as.data.frame(cbind(rt[,(1:8)],a))
# 
# validSet <- validSet[,c(sigGenes,lasso_fea)]
# # validSet$riskScore <- apply(validSet[,lasso_fea], 1, function(x) {x %*% multiCoxSum$coefficients[,"coef"]})
# validSet$riskScore <-  apply(validSet[,lasso_fea], 1, function(x) {x %*% myCoefs@x})
# 
# validSet$risk=as.vector(ifelse(validSet$riskScore>median(validSet$riskScore),"High Risk","Low Risk"))
# 
# write.table( cbind(id=rownames(validSet),validSet),
#              file="04LASSO/testSet_riskScore.xls",
#              sep="\t",row.names = F,
#              quote=F)



#GEO外部验证集
# load("D:/GEO/GSE28221-GPL5175_series_matrix.txt.gz")
clincial_exprSet <- read.table("01rawData/GSE28221_removeBatchEffect_expres.xls",header = T,row.names = 1,sep = "\t",check.names = F)
clincal_raw <- read.delim("01rawData/GSE28221_readme.txt",check.names = F,header = T,sep = "\t")
GEO_validSet1 <- clincial_exprSet[coxGene,]
GEO_validSet1 <- as.data.frame( t(GEO_validSet1))
GEO_validSet1$id <- rownames(GEO_validSet1)
# exprSet$id <- substr(rownames(exprSet),1,12)
exprSet_Time <- merge(clincal_raw,GEO_validSet1,by="id")
# exprSetTime <- exprSetTime[,c(ncol(exprSetTime),2:(ncol(exprSetTime)-1))]
write.table(exprSet_Time,file = "05Cox/GSE28221_modelGenes_exprsetTime.xls",row.names = F,quote = F,sep = "\t")
clincial_exprSet
exprSet_Time <- clincial_exprSet
exprSet_Time <-  read.table("05Cox/GSE28221_modelGenes_exprsetTime.xls",header = T,row.names = 1,sep = "\t",check.names = F)
# sigGenes=c("`OS (months)`","`OS-censor`","age","gender","pathologic_stage","pathologic_t","pathologic_m","pathologic_n")
sigGenes=c("futime","fustat")
# GEO_coxGene  <-  gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\", coxGene)
# coxGene %in% colnames(GEO_validSet1)
# "LDLCQ4"  %in% colnames(GEO_validSet1)
GEO_validSet <- exprSet_Time[,c(sigGenes,coxGene)]
# GEO_validSet <- exprSet_Time[,c(sigGenes,lasso_fea)]
GEO_validSet$riskScore <- apply(GEO_validSet[,coxGene], 1, function(x) {x %*% multiCoxSum$coefficients[,"coef"]})
# GEO_validSet$riskScore <-  apply(GEO_validSet[,lasso_fea], 1, function(x) {x %*% myCoefs@x})
# GEO_validSet$risk=as.vector(ifelse(GEO_validSet$riskScore>median(GEO_validSet$riskScore),"High Risk","Low Risk"))
library(survminer)
#对数据集的基因进行bestSeparation统计
res.cut <- surv_cutpoint(GEO_validSet, time = "futime",
                         event = "fustat",
                         variables = c("riskScore"),
                         minprop = 0.3) #
# cutpoint statistic
# riskScore 3.656151  1.752779
res.cat <- surv_categorize(res.cut)
GEO_validSet$risk=as.vector(ifelse(res.cat$riskScore == "high","High Risk","Low Risk"))

# rt1$risk=as.vector(ifelse(rt1$riskScore>median(rt1$riskScore),"High Risk","Low Risk"))

write.table( cbind(id=rownames(GEO_validSet),GEO_validSet),
             file="05Cox/GSE28221_validSet_riskScore.xls",
             sep="\t",row.names = F,
             quote=F)





