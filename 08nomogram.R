
library(survival)
library(regplot)
library(rms)
library(nomogramEx)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
pbc=read.table("05Cox/trainSet_riskScore.xls",sep="\t",header=T,row.names=1,check.names=F)       #??ȡ?????ļ?
# pbc$figo_Stage <- as.numeric(as.factor(pbc$figo_Stage))
# pbc$Type <- as.numeric(as.factor(pbc$Type ))
# pbc$group <- factor(pbc$group ,levels = c("PTEN wild","PTEN mutation"))
dd <- datadist(pbc)
options(datadist="dd")
options(na.action="na.delete")
summary(pbc$futime)
coxpbc <- cph(formula = Surv(futime,fustat) ~  IL1R2 +S100A12 +CCL8 ,data=pbc,x=T,y=T,surv = T,na.action=na.delete)  #,time.inc =2920

print(coxpbc)

surv <- Survival(coxpbc) 
surv1 <- function(x) surv(365,x)
surv2 <- function(x) surv(730,x)
surv3 <- function(x) surv(1095,x)
# surv10 <- function(x) surv(356*10,x)

x <- nomogram(coxpbc,fun = list(surv1,surv2,surv3),lp=T,
              funlabel = c('1-year survival Probability','2-year survival Probability','3-year survival Probability'),
              maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.5,0.3,0.1))

pdf("08nomogram/nomogram_classical.pdf",width = 10,height = 7)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

# 绘制calibration curve进行验证
f1 <- cph(formula =  Surv(futime,fustat) ~IL1R2 +S100A12 +CCL8,data=pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 
#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=50,B=1000) 

f3 <- cph(formula =  Surv(futime,fustat) ~ IL1R2 +S100A12 +CCL8,data=pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 730) 
#参数m=50表示每组50个样本进行重复计算  365*3 = 1095
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=730,m=50,B=1000) 

f5 <- cph(formula =  Surv(futime,fustat) ~ IL1R2 +S100A12 +CCL8,data=pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 
#参数m=50表示每组50个样本进行重复计算   365*5 = 1825
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1095,m=50,B=1000) 

pdf("08nomogram/calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#F0027F"),
     xlim = c(0,1),ylim= c(0,1),col = c("#F0027F"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#F0027F"), pch = 16)

mtext("")

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)


abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year"), #图例文字
       col =c("#2166AC","#F0027F","#B2182B"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()



pdf("08nomogram/calibration_1y.pdf",width = 8,height = 8)
plot(cal1,
     lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("#2166AC"),#error bar的颜色
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小
lines(cal1[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("#2166AC")) #连线的颜色
mtext("")
box(lwd = 1) #边框粗细
abline(0,1,lty = 3, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("#224444")#对角线的颜色
) 
dev.off()

pdf("08nomogram/calibration_2y.pdf",width = 8,height = 8)
plot(cal3,
     lwd = 2,
     lty = 1,
     errbar.col = c("#F0027F"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#F0027F"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal3[,c('mean.predicted',"KM")],
      type= 'b',
      lwd = 2,
      col = c("#F0027F"),
      pch = 16)
mtext("")
box(lwd = 1)
abline(0,1,lty= 3,
       lwd = 2,
       col =c("#224444"))
dev.off()


pdf("08nomogram/calibration_3y.pdf",width = 8,height = 8)
plot(cal5,
     lwd = 2,
     lty = 1,
     errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#B2182B"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal5[,c('mean.predicted',"KM")],
      type= 'b',
      lwd = 2,
      col = c("#B2182B"),
      pch = 16)
mtext("")
box(lwd = 1)
abline(0,1,lty= 3,
       lwd = 2,
       col =c("#224444"))
dev.off()



#C index
library(compareC)
f <- coxph(Surv(futime,fustat) ~ IL1R2+S100A12 +CCL8,data=pbc)###计算出的C-index
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index
plot(sum.surv)
# C      se(C) 
# 0.70866888 0.02848457 

