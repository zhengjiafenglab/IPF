rm(list = ls())
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

                
rt=read.table("04venn/veen_intersectGenes.xls",sep="\t",header=F,check.names=F)
rt$V1[5] <- "CXCL8"
entrezIDs <- mget(rt$V1 ,org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)



entrezID_gene <- na.omit(entrezIDs)


kk <- enrichGO(gene = entrezID_gene,
               OrgDb = org.Hs.eg.db,
               pvalueCutoff =0.05,
               qvalueCutoff = 0.05,
               ont="all",
               readable =T
)

write.table(kk,file="04venn/GO.xls",quote=F,row.names = F,sep = "\t")

pdf(file="04venn/GO_barplot.pdf",width = 9,height = 9)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY",label_format = 30) +
  # scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  facet_grid(ONTOLOGY~., scale='free')
dev.off()


pdf(file="04venn/GO_bubble.pdf",width = 9,height = 9)
dotplot(kk,showCategory = 10,split="ONTOLOGY",label_format = 30) +
  # scale_x_discrete(labels=function(x) str_wrap(x, width=10))+
  facet_grid(ONTOLOGY~., scale='free')
dev.off()



kk <- enrichKEGG(gene = entrezID_gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)


y <- setReadable(kk,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
write.table(y,file="04venn/KEGG.xls",quote=F,row.names = F,sep = "\t")

pdf(file="04venn/KEGG_barplot.pdf",width = 7.5,height = 3.5)
barplot(kk, drop = TRUE, showCategory = 20)
dev.off()


pdf(file="04venn/KEGG_bubble.pdf",width = 7.5,height = 3.5)
dotplot(kk, showCategory = 20)
dev.off()






rm(list = ls())
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")


rt=read.table("02GSVA/GSE70866_EMTScore_diff.xls",sep="\t",header=T,check.names=F,row.names = 1)
entrezIDs <- mget(rownames(rt), org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
out = out[,c(1,7)]
out = cbind(symbol=rownames(out),out)
write.table(out,file="02GSVA/id.xls",quote=F,row.names = F,sep="\t")

rt <- read.table("02GSVA/id.xls",sep="\t",header = T)
entrezID_gene <- na.omit(rt$entrezID)



kk <- enrichGO(gene = entrezID_gene,
               OrgDb = org.Hs.eg.db,
               pvalueCutoff =0.05,
               qvalueCutoff = 0.05,
               ont="all",
               readable =T
)

write.table(kk,file="02GSVA/GO.xls",quote=F,row.names = F,sep = "\t")

pdf(file="02GSVA/GO_barplot.pdf",width = 9,height = 9)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY",label_format = 30) +
  # scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  facet_grid(ONTOLOGY~., scale='free')
dev.off()


pdf(file="02GSVA/GO_bubble.pdf",width = 9,height = 9)
dotplot(kk,showCategory = 10,split="ONTOLOGY",label_format = 30) +
  # scale_x_discrete(labels=function(x) str_wrap(x, width=10))+
  facet_grid(ONTOLOGY~., scale='free')
dev.off()



kk <- enrichKEGG(gene = entrezID_gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)


y <- setReadable(kk,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
write.table(y,file="02GSVA/KEGG.xls",quote=F,row.names = F,sep = "\t")

pdf(file="02GSVA/KEGG_barplot.pdf",width = 7.5,height = 3.5)
barplot(kk, drop = TRUE, showCategory = 20)
dev.off()


pdf(file="02GSVA/KEGG_bubble.pdf",width = 7.5,height = 3.5)
dotplot(kk, showCategory = 20)
dev.off()



####################
# GSEA分析
# GSEA GO富集分析
library("clusterProfiler")
library("org.Hs.eg.db")
DEGs <- read.csv("highlowRisk_all.csv",header = TRUE,row.names = 1)

DEGs <- DEGs[order(DEGs$logFC,decreasing = T),]
entrezIDs <- mget(rownames(DEGs), org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(DEGs,entrezID=entrezIDs)
out = out[,c(1,7)]
out <- cbind(symbol = rownames(out),out)
write.table(out,file="GSEA_id.xls",sep="\t",quote=F,row.names = F)  

GSEA_data <- read.table("GSEA_id.xls",sep="\t",header = T) 
GSEA_data <- na.omit(GSEA_data)
GSEA_gene <- GSEA_data$logFC
names(GSEA_gene) <- GSEA_data$entrezID

# ggsea <-  gseGO(
#   GSEA_gene,
#   ont = "BP",
#   OrgDb = org.Hs.eg.db,
#   keyType = "ENTREZID",
#   pvalueCutoff = 0.05,
#   eps = 0,
#   verbose = TRUE,
#   seed = FALSE
# 
# )
ggsea <- gseGO(geneList     = GSEA_gene,
              OrgDb        = org.Hs.eg.db,
              ont          = "MF",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)



ggsea.symbol <- setReadable(ggsea,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
# ridgeplot(ggsea, 10)
# ggsave("F:/YQD007/cutoff0.5/Result/3.Enrichment/GOBP_GSEAps0_DEG_all.png",width = 8,height = 7,dpi=600)
# ggsave("F:/YQD007/cutoff0.5/Result/3.Enrichment/GOBP_GSEAps0_DEG_all.pdf",width = 8,height = 7)

gseaplot2(ggsea, 1:6)
ggsave("GOMF_GSEA_DEG_all.png",width = 10,height = 8,dpi=600)
ggsave("GOMF_GSEA_DEG_all.pdf",width = 10,height = 8)
write.csv(ggsea.symbol,file = "GOMF_GSEA_DEG_all.csv",row.names = F)

###########GESA KEGG enrichment

# kegggse <-  gseKEGG(
#   GSEA_gene,
#   organism = "hsa",
#   keyType = "kegg",
# 
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   eps = 0,
#   verbose = TRUE,
#   seed = FALSE
# )

kegggse <- gseKEGG(geneList     = GSEA_gene,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 100,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
kegggse.symbol <- setReadable(kegggse,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

gseaplot2(kegggse,c(4,17,48,49,55,94))
ggsave("KEGG_GSEA_DEG_all.png",width = 10,height = 8,dpi=600)
ggsave("KEGG_GSEA_DEG_all.pdf",width = 10,height = 8)
write.csv(kegggse.symbol,file = "KEGG_GSEA_DEG_all.csv",row.names = F)


