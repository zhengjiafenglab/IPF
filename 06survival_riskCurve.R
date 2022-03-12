
#生存分析
library(survival)
library(survminer)
rt <- read.table("05Cox/trainSet_riskScore.xls",header=T,sep="\t")
# rt <- read.table("04LASSO/TCGA-STAD_lncRNA_cox_input.xls",header=T,sep="\t",check.names = F)
# group <- ifelse(rt$ENSG00000237094.lincRNA.RP4.669L17.10 <= median(rt$mRNAsi ) , "Low mRNAsi","High mRNAsi")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,length(diff$n) - 1)
HR = (diff$obs[2]/diff$exp[2])/(diff$obs[1]/diff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
# 
# HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
# CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")

paste0(round(HR,3),"(" ,paste(round(low95,3), round(up95,3), sep = " - "),")" )

median(rt$futime)
pval <- ifelse(pValue < 0.001, "p < 0.001", 
               paste("p = ",round(pValue,3), sep = ""))
pval


fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

ggsurvplot(fit, data = rt ,
           ggtheme = theme_classic(), #想要网格就运行这行
           
           conf.int = T, #不画置信区间，想画置信区间就把F改成T
           #conf.int.style = "step",#置信区间的类型，还可改为ribbon
           risk.table = TRUE, # 绘制累计风险曲线
           #tables.theme = theme_void(),
           tables.theme = theme_classic(),
           tables.height = 0.2,
           #surv.median.line = "hv",
           censor = T, #不显示观察值所在的位置
           #ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           #ncensor.plot.height = 0.2,
           palette = c("#FF0033","#1B9E77"), #线的颜色对应高、低(自定义调色板) FF0033红色
           xlab = "Time (Days)",
           ylab = "Survival probability",
           legend.title = "",#基因名写在图例题目的位置
           font.legend = 11,#图例的字体大小
           #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
           legend=c(0.75,0.85),
           legend.labs = c("High Risk","Low Risk"),
           #在图例上标出高低分界点的表达量，和组内sample数量
           #legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
           #paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
           
           #标出pvalue、HR、95% CI
           #太小的p value标为p < 0.001
           # pval = paste(pval =pval ,
           #              HR, CI, sep = "\n")
           pval = pval)
dev.copy2pdf(file = "06ROC/trainSet_survival.pdf", width = 5.5,height = 5.5)
dev.off()


library(survivalROC)
pdf(file="06ROC/trainSet_ROC.pdf",width=5,height=5)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc1=survivalROC(Stime=rt$futime,
                 status=rt$fustat,
                 marker = rt$riskScore ,
                 predict.time =365 ,
                 #span = 0.25*nrow(rt)^(-0.20)##span,NNE法的namda
                 method="KM")
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#1B9E77",
     xlab="1-specificity", ylab="sensitivity",
     main="ROC curve, Method = KM",
     lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# abline(0,1,lty= 3)
# legend("bottomright", legend=paste0("1 year AUC = ",sprintf("%.3f",roc1$AUC),")"),lwd=2,bty="n",col='#FF0033')
roc3=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore,
                 predict.time =365*2, method="KM")

lines(roc3$FP, roc3$TP, type="l",col="#E8822E",xlim=c(0,1), ylim=c(0,1),lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)

roc5=survivalROC(Stime=rt$futime, status=rt$fustat, marker =rt$riskScore,
                 predict.time =365*3, method="KM")
lines(roc5$FP, roc5$TP, type="l",col="#FF0033",xlim=c(0,1), ylim=c(0,1),lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1,lty= 3)
legend(0.5,0.18,c(paste("1 year AUC = ",sprintf("%.3f",roc1$AUC)),paste("2 year AUC = ",sprintf("%.3f",roc3$AUC)),
                  paste("3 year AUC = ",sprintf("%.3f",roc5$AUC))),
       x.intersp=0.5, y.intersp=0.85,
       lty= 1 ,lwd= 2,col=c("#1B9E77","#E8822E","#FF0033"),
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)#
# legend(0.5,0.18,c(paste("1 year AUC = ",sprintf("%.3f",roc1$AUC)),paste("3 year AUC = ",sprintf("%.3f",roc3$AUC))),
#        x.intersp=0.5, y.intersp=0.85,
#        lty= 1 ,lwd= 2,col=c("#1B9E77","#E8822E"),
#        bty = "n",# bty框的类型
#        seg.len=1,cex=0.8)#
dev.off()

## 训练集风险曲线


library(pheatmap)

rt=read.table("05Cox/trainSet_riskScore.xls",sep="\t",header=T,row.names=1,check.names=F)       #??ȡ?????ļ?
rt=rt[order(rt$riskScore),]                                     #????riskScore????Ʒ????

#???Ʒ???????
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="Low Risk"])
highLength=length(riskClass[riskClass=="High Risk"])
line=rt[,"riskScore"]
# line[line>10]=10
pdf(file="06ROC/trainSet_riskScore.pdf",width = 10,height = 4)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#028846",lowLength),
           rep("red",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High Risk", "Low Risk"),bty="n",pch=19,col=c("red","#028846"),cex=1.2)
dev.off()


# ggData1 <- data.frame(riskScore=line,patient=1:length(line))
# p1 <- ggplot(data=ggData1,mapping = aes(x=patient,y=riskScore))+
#         geom_point(col=c(rep("#028846",lowLength), rep("red",highLength)))+
#         theme_classic()+scale_color_manual(name="HH",)+
#         theme(axis.line.x = element_blank(),
#               axis.text.x = element_blank(),
#               axis.ticks.x = element_blank(),
#               axis.title.x = element_blank(),
#               panel.border = element_blank(),
#               panel.background = element_blank())+
#         ylab("Risk Score")+
#         geom_vline(xintercept =lowLength,linetype =2)







#????????״̬ͼ
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="#028846"
pdf(file="06ROC/trainSet_survStat.pdf",width = 10,height = 4)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (Days)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","#028846"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

#???Ʒ?????ͼ
rt1=rt[c(3:(ncol(rt)-2))]
max(rt1) 
min(rt1)

# rt1 <- log2(rt1+1)
# rt1[rt1 > 3] =3
# hmExp[hmExp < 8] =8
rt1=t(rt1)
# rownames(rt1) <- gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\",rownames(rt1))


annotation=data.frame(Type=(rt[,ncol(rt)]))
rownames(annotation)=rownames(rt)
ann_colors = list(Type =c(`Low Risk`="#00BFC4", `High Risk`="#F8766D"))


pdf(file="06ROC/trainSet_heatmap.pdf",width = 10,height = 3)
pheatmap(rt1, 
         scale = "row",
         annotation_col=annotation, 
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         border=FALSE,
         color = colorRampPalette(c("#028846", "white", "red"))(50) )
dev.off()



########### GEO 验证集
#生存分析
library(survival)
library(survminer)
rt <- read.table("05Cox/GSE28221_validSet_riskScore.xls",header=T,sep="\t")
# group <- ifelse(rt$"PYGM" <= median(rt$"PYGM" ) , "Low","High")
# diff=survdiff(Surv(futime, fustat) ~group,data = rt)
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,length(diff$n) - 1)
HR = (diff$obs[2]/diff$exp[2])/(diff$obs[1]/diff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
# 
# HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
# CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")

paste0(round(HR,3),"(" ,paste(round(low95,3), round(up95,3), sep = " - "),")" )

median(rt$futime)
pval <- ifelse(pValue < 0.001, "p < 0.001", 
               paste("p = ",round(pValue,3), sep = ""))
pval


fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

ggsurvplot(fit, data = rt ,
           ggtheme = theme_classic(), #想要网格就运行这行
           
           conf.int = T, #不画置信区间，想画置信区间就把F改成T
           #conf.int.style = "step",#置信区间的类型，还可改为ribbon
           risk.table = TRUE, # 绘制累计风险曲线
           #tables.theme = theme_void(),
           tables.theme = theme_classic(),
           tables.height = 0.2,
           #surv.median.line = "hv",
           censor = T, #不显示观察值所在的位置
           #ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           #ncensor.plot.height = 0.2,
           palette = c("#FF0033","#1B9E77"), #线的颜色对应高、低(自定义调色板) FF0033红色
           xlab = "Time (Days)",
           ylab = "Survival probability",
           legend.title = "",#基因名写在图例题目的位置
           font.legend = 11,#图例的字体大小
           #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
           legend=c(0.75,0.85),
           legend.labs = c("High Risk","Low Risk"),
           #在图例上标出高低分界点的表达量，和组内sample数量
           #legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
           #paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
           
           #标出pvalue、HR、95% CI
           #太小的p value标为p < 0.001
           # pval = paste(pval =pval ,
           #              HR, CI, sep = "\n")
           pval = pval)
dev.copy2pdf(file = "06ROC/GSE28221_survival.pdf", width = 5.5,height = 5.5)
dev.off()



library(survivalROC)

pdf(file="06ROC/GSE28221_ROC.pdf",width=5,height=5)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc1=survivalROC(Stime=rt$futime,
                 status=rt$fustat,
                 marker = rt$riskScore ,
                 predict.time =1 ,
                 #span = 0.25*nrow(rt)^(-0.20)##span,NNE法的namda
                 method="KM")
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#1B9E77",
     xlab="1-specificity", ylab="sensitivity",
     main="ROC curve, Method = KM",
     lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# abline(0,1,lty= 3)
# legend("bottomright", legend=paste0("1 year AUC = ",sprintf("%.3f",roc1$AUC),")"),lwd=2,bty="n",col='#FF0033')
roc3=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore,
                 predict.time =1*2, method="KM")

lines(roc3$FP, roc3$TP, type="l",col="#E8822E",xlim=c(0,1), ylim=c(0,1),lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)

roc5=survivalROC(Stime=rt$futime, status=rt$fustat, marker =rt$riskScore,
                 predict.time =1*3, method="KM")
lines(roc5$FP, roc5$TP, type="l",col="#FF0033",xlim=c(0,1), ylim=c(0,1),lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1,lty= 3)
legend(0.5,0.18,c(paste("1 year AUC = ",sprintf("%.3f",roc1$AUC)),paste("2 year AUC = ",sprintf("%.3f",roc3$AUC)),
                  paste("3 year AUC = ",sprintf("%.3f",roc5$AUC))),
       x.intersp=0.5, y.intersp=0.85,
       lty= 1 ,lwd= 2,col=c("#1B9E77","#E8822E","#FF0033"),
       bty = "n",# bty框的类型
       seg.len=1,cex=0.8)#
# legend(0.5,0.18,c(paste("1 year AUC = ",sprintf("%.3f",roc1$AUC)),paste("3 year AUC = ",sprintf("%.3f",roc3$AUC))),
#        x.intersp=0.5, y.intersp=0.85,
#        lty= 1 ,lwd= 2,col=c("#1B9E77","#E8822E"),
#        bty = "n",# bty框的类型
#        seg.len=1,cex=0.8)#
dev.off()


## 验证集风险曲线
library(pheatmap)

# rt=read.table("04LASSO/GSE30219_validSet_riskScore.xls",sep="\t",header=T,row.names=1,check.names=F)       #??ȡ?????ļ?
rt=rt[order(rt$riskScore),]                                     #????riskScore????Ʒ????

#???Ʒ???????
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="Low Risk"])
highLength=length(riskClass[riskClass=="High Risk"])
line=rt[,"riskScore"]
# line[line>10]=10
pdf(file="06ROC/GSE28221_riskScore.pdf",width = 10,height = 4)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#028846",lowLength),
           rep("red",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High Risk", "Low Risk"),bty="n",pch=19,col=c("red","#028846"),cex=1.2)
dev.off()

#????????״̬ͼ
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="#028846"
pdf(file="06ROC/GSE28221_survStat.pdf",width = 10,height = 4)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (Days)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","#028846"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

#???Ʒ?????ͼ
rt1=rt[c(4:(ncol(rt)-2))]
max(rt1) 
min(rt1)
# rt1[rt1 > 80] =80
# rt1 <- log2(rt1+1)

# hmExp[hmExp < 8] =8
rt1=t(rt1)
rownames(rt1) <- gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\",rownames(rt1))


annotation=data.frame(Type=(rt[,ncol(rt)]))
rownames(annotation)=rownames(rt)
ann_colors = list(Type =c(`Low Risk`="#00BFC4", `High Risk`="#F8766D"))


pdf(file="06ROC/GSE28221_heatmap.pdf",width = 10,height = 3)

pheatmap(rt1, 
         scale = "row",
         annotation_col=annotation, 
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         border=FALSE,
         color = colorRampPalette(c("#028846", "white", "red"))(50) )
dev.off()

