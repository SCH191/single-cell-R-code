library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(Seurat)
library(GSVA)
library(pheatmap)
library(patchwork)
library(fgsea)
library(tibble)
library(fanyi)
library(EnhancedVolcano)


#ident1 vs ident2 单个图###########################################
OrgDb1="org.Hs.eg.db"   #鼠：org.Mm.eg.db 人：org.Hs.eg.db
top = 10#选取通路数
nident.1 = '0'
nident.2 = 'Control'
ABC = ABC
#auto run
library(readxl)

markers = FindMarkers(ABC,ident.1 = nident.1,ident.2 = nident.2,only.pos = F,logfc.threshold = 0.25,min.pct = 0.25)
up = subset(markers,subset = avg_log2FC>0.25)
down = subset(markers,subset = avg_log2FC< (-0.25))

up$avg_log2FC = abs(up$avg_log2FC)
down$avg_log2FC = abs(down$avg_log2FC)

x = rownames(up)
tmp=gsub("\\..*","",x)
eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)
genelist <- eg$ENTREZID
go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
up = go[which(go@result$ONTOLOGY!='CC'),c(3,10,8)]

x = rownames(down)
tmp=gsub("\\..*","",x)
eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)
genelist <- eg$ENTREZID
go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
down = go[which(go@result$ONTOLOGY!='CC'),c(3,10,8)]

up$change = c('up')
down$change = c('down')

#
down$N = c('')
up$N = c('')
for (i in 1:length(up$Description)) {
  for (j in 1:length(down$Description)){
    if (rownames(up)[i] == rownames(down)[j]){
      up[i,'N']='p'
      down[j,'N']='p'
    }
    
  }
}
up = subset(up,subset=N!='p')
down = subset(down,subset=N!='p')
#

up = top_n(up,n = top, wt = qvalue)
down = top_n(down,n = top, wt = qvalue)
up$N=NULL
down$N = NULL
bp <- rbind(up, down)
colnames(bp) <- c("Des", "Count", "Qvalue", "Change")

data <- bp %>% mutate(Qvalue2 = ifelse(Change == "up", -log10(Qvalue), log10(Qvalue))) %>% arrange(Change,Qvalue2) 
data <- data %>% mutate(Count = ifelse(Change == "up", Count, -Count))
data$Des = factor(data$Des, levels = unique(data$Des),ordered = T)
head(data)


tmp <- with(data, labeling::extended(-range(Qvalue2)[2], range(Qvalue2)[2], m = 7))

lm <- tmp[c(1,length(tmp))]

# 绘图
a = ggplot(data, aes(x=Des, y=Log10qvalue)) +
  geom_segment( aes(x=Des, xend=Des, y=0, yend=Qvalue2, color=Change), size=5, alpha=0.9) +
  theme_light() +
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 15),
     legend.text = element_text(size = 20),
    # legend.position = c(0.9, 0.5),
    # 不显示图例
    legend.position = "none"
  ) + coord_flip()
  +xlab("") +
  ylab("-log10(Qvalue)")+
  ylim(lm)+
  scale_y_continuous(breaks = tmp,
                     labels = abs(tmp),
                     # 双Y轴
                     sec.axis = sec_axis(~./0.3,
                                         name = 'Number of genes',
                                         breaks = seq(-100,100,10),
                                         labels = abs(seq(-100,100,10))))+
  # 翻转一下
   #+ 
  # 绘制折线图
  #geom_line(aes(Des, Count*0.3, group = 1), size=0.8, color='#71ad46') + geom_point(aes(Des, Count*0.3, group = 1), color='#71ad46')

pdf(paste(nident.1,'_vs_',nident.2,'_GO.pdf',sep = ''),width = 10,height = 6)
a
dev.off()
#ident1 vs ident2 分两张图######################################################
OrgDb1="org.Mm.eg.db"   #鼠：org.Mm.eg.db 人：org.Hs.eg.db
top = 20#选取通路数
nident.1 = '0'
nident.2 = '1'
top_use = 0 #使用topn基因用作富集，0则使用全部符合的基因
ABC = ABC
#auto run
markers = FindMarkers(ABC,ident.1 = nident.1,ident.2 = nident.2,only.pos = F,logfc.threshold = 0.25,min.pct = 0.25)
up = subset(markers,subset = avg_log2FC>0.25&p_val<0.05)
down = subset(markers,subset = avg_log2FC< (-0.25)&p_val<0.05)

up$avg_log2FC = abs(up$avg_log2FC)
down$avg_log2FC = abs(down$avg_log2FC)
if (top_use!=0) {
up = top_n(up,n = top_use, wt = avg_log2FC)
down = top_n(down,n = top_use, wt = avg_log2FC)
}
dir.create(paste(nident.1,'_vs_',nident.2,'_GO',sep = ''))
setwd(paste('./',nident.1,'_vs_',nident.2,'_GO',sep = ''))
write.csv(markers,'DEG.csv')
dir.create(nident.1)
setwd(paste('./',nident.1,sep=''))
x = rownames(down)
tmp=gsub("\\..*","",x)
eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)
genelist <- eg$ENTREZID
go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
write.csv(go@result,'result.csv')
go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
try(ggsave("1.1.1.GO_barplot.pdf",barplot(go, split="ONTOLOGY",showCategory = 15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
try(ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
try(ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
try(ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
try(ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
try(ggsave("1.1.1.GO_barplot.png",barplot(go, split="ONTOLOGY",showCategory=15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
try(ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
try(ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
try(ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
try(ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))

try(ggsave('3.1GO_gene_Connection.png',enrichplot::cnetplot(go,circular=FALSE,colorEdge = TRUE),width=14,height=10))#基因-通路关联网络图
try(ggsave('3.2KEGG_gene_Connection.png',enrichplot::cnetplot(kegg,circular=FALSE,colorEdge = TRUE),width=14,height=10))#circluar为指定是否环化，基因过多时建议设置为FALSE
try(GO2 <- pairwise_termsim(go))
try(KEGG2 <- pairwise_termsim(kegg))
try(ggsave('4.1GO_functions_Connection.png',enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk"),width=25,height=25))
try(ggsave('4.2KEGG_functions_Connection.png',enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk"),width=25,height=25))
setwd('../')

dir.create(nident.2)
setwd(paste('./',nident.2,sep=''))
x = rownames(up)
tmp=gsub("\\..*","",x)
eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)
genelist <- eg$ENTREZID
go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
write.csv(go@result,'result.csv')
go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
try(ggsave("1.1.1.GO_barplot.pdf",barplot(go, split="ONTOLOGY",showCategory = 15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
try(ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
try(ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
try(ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
try(ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
try(ggsave("1.1.1.GO_barplot.png",barplot(go, split="ONTOLOGY",showCategory=15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
try(ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
try(ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
try(ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
try(ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))

try(ggsave('3.1GO_gene_Connection.png',enrichplot::cnetplot(go,circular=FALSE,colorEdge = TRUE),width=14,height=10))#基因-通路关联网络图
try(ggsave('3.2KEGG_gene_Connection.png',enrichplot::cnetplot(kegg,circular=FALSE,colorEdge = TRUE),width=14,height=10))#circluar为指定是否环化，基因过多时建议设置为FALSE
try(GO2 <- pairwise_termsim(go))
try(KEGG2 <- pairwise_termsim(kegg))
try(ggsave('4.1GO_functions_Connection.png',enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk"),width=25,height=25))
try(ggsave('4.2KEGG_functions_Connection.png',enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk"),width=25,height=25))
setwd('../')
tampt = subset(merge , subset = manusubcelltype == nident.1|manusubcelltype == nident.2)
tampt <- NormalizeData(tampt)
tampt <- FindVariableFeatures(tampt, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(tampt)
tampt <- ScaleData(tampt, features = all.genes)
ABC.markers <- FindAllMarkers(tampt, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
ABC.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10
a = DoHeatmap(tampt, features = top10$gene) + NoLegend()
try(ggsave('heatmap.png',a,width=20,height=20))#基因-通路关联网络图
setwd('../')


#ident vs all 批处理######################################################
ABC =  ABC
DefaultAssay(ABC) = 'RNA'
Idents(ABC) = ABC@meta.data$manucelltype
OrgDb1="org.Hs.eg.db"   #鼠：org.Mm.eg.db 人：org.Hs.eg.db
top_use = 0 #使用topn基因用作富集，0则使用全部符合的基因
CN = TRUE                 #是否添加中文翻译结果，TRUE添加，FALSE不添加

#

dir.create('Diff_Group_By_manucelltype')
setwd(paste("./",'Diff_Group_By_manucelltype',sep = ''))
markers <- FindAllMarkers(ABC, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers , "marker.csv")

markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10
Heatmap_15987 = DoHeatmap(ABC, features = top10$gene) + NoLegend()
ggsave('heatmap.png',Heatmap_15987,width=20,height=15)

for (i in 1:length(levels(ABC))) {
  dir.create(as.character(levels(ABC)[i]))
  setwd(paste("./",as.character(levels(ABC)[i]),sep = ''))
  Focus=rownames( markers [intersect(intersect(which(markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=0.5)),which( as.character(markers [,"cluster"])==as.character(levels(ABC)[i]))),])
  Focus_no_only.pos=rownames( markers [intersect(which(markers [,"p_val"]<0.05),which( as.character(markers [,"cluster"])==as.character(levels(ABC)[i]))),])

  x = markers[Focus,'gene']
  if (top_use!=0) {
    x = markers[rownames(top_n(markers[x,],n = top_use, wt = avg_log2FC)),'gene']
  }
  tmp=gsub("\\..*","",x)
  eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  write.csv(go@result,'go_result.csv')
  go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  try(ggsave("1.1.1.GO_barplot.pdf",barplot(go, split="ONTOLOGY",showCategory = 15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
  try(ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
  kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  try(ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
  try(ggsave("1.1.1.GO_barplot.png",barplot(go, split="ONTOLOGY",showCategory=15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
  try(ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
  if(CN){  try(ggsave("1.1.1.GO_barplot_ZH.png",translate_ggplot(barplot(go, split="ONTOLOGY",showCategory = 15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),axis = 'y',from = 'en',to='zh'),width=14,height=10))
    try(ggsave("2.1.KEGG_dotplot_ZH.png",translate_ggplot(dotplot(kegg, showCategory=30,label_format=100),axis = 'y',from = 'en',to='zh'),width=14,height=10))
  }
  try(ggsave('3.1GO_gene_Connection.png',enrichplot::cnetplot(go,circular=FALSE,colorEdge = TRUE),width=14,height=10))#基因-通路关联网络图
  try(ggsave('3.2KEGG_gene_Connection.png',enrichplot::cnetplot(kegg,circular=FALSE,colorEdge = TRUE),width=14,height=10))#circluar为指定是否环化，基因过多时建议设置为FALSE
  try(GO2 <- pairwise_termsim(go))
  try(KEGG2 <- pairwise_termsim(kegg))
  try(ggsave('4.1GO_functions_Connection.png',enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk"),width=25,height=25))
  try(ggsave('4.2KEGG_functions_Connection.png',enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk"),width=25,height=25))
  

  
  setwd('../')
}
setwd('../')


#按细胞类型分别比较实验组对照组差异表达基因富集################################################
dir.create('CellType_Diff')
setwd('./CellType_Diff')
ABC@meta.data$test = paste(ABC@meta.data$Group_MTOR_PIK3_TSC,ABC@meta.data$orig.ident,sep='__')
ABC@meta.data$test = factor(ABC@meta.data$test,levels = sort(levels(as.factor(ABC@meta.data$test))))
tampt = ABC

DefaultAssay(tampt) = 'RNA'
tampt <- NormalizeData(tampt)
tampt <- FindVariableFeatures(tampt, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(tampt)
tampt <- ScaleData(tampt, features = all.genes)


divide_by = 'manucelltyp' #细胞分组划分参考
text_ident = 'Group'       #差异对比参考
disease.type = 'TLE'  #差异组1
control.type = 'Control'   #差异组2
cutoff = 0.5               #富集基因差异认定阈值
OrgDb1="org.Mm.eg.db"      #样本参考集数据，鼠：'org.Mm.eg.db'，人：'org.Hs.eg.db'
CN = TRUE                 #是否添加中文翻译结果，TRUE添加，FALSE不添加
heat_ident = 'Group'  #热图分组依据


#tampt@meta.data$test = paste(tampt@meta.data$Group_MTOR,tampt@meta.data$orig.ident,sep='__')
#tampt@meta.data$test = factor(tampt@meta.data$test,levels = sort(levels(as.factor(tampt@meta.data$test))))

eval(parse(text = paste('tempt_diff = SplitObject(tampt,split.by = "',divide_by,'")',sep = '')))
rm(tampt)
gc()
for (i in 1:length(tempt_diff)){
  gc()
  dir.create(names(tempt_diff[i]))
  setwd(paste('./',names(tempt_diff[i]),sep = ''))
  flag = TRUE
  #Idents(tempt_diff[[i]]) = tempt_diff[[i]]@meta.data$Group
  eval(parse(text = paste('Idents(tempt_diff[[i]]) = tempt_diff[[i]]@meta.data$',text_ident,sep = '')))
  tempt_diff[[i]] <- NormalizeData(tempt_diff[[i]])
  tempt_diff[[i]] <- FindVariableFeatures(tempt_diff[[i]], selection.method = "vst", nfeatures = 4000)
  all.genes <- rownames(tempt_diff[[i]])
  tempt_diff[[i]] <- ScaleData(tempt_diff[[i]], features = all.genes)
  tryCatch({
    markers <- FindAllMarkers(tempt_diff[[i]], only.pos = F, min.pct = 0.1, logfc.threshold = 0.1)
  },warning = function(w){
    flag<<-FALSE
    print("warning")
  },error = function(e){
    flag<<-FALSE
    print("error")
  })
  if(!flag){
    setwd('../')
    flag = TRUE
    next
  }
  df = subset(markers , subset = cluster==disease.type)
  markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10
  #Idents(tempt_diff[[i]]) = tempt_diff[[i]]@meta.data$orig.ident
  eval(parse(text = paste('Idents(tempt_diff[[i]]) = tempt_diff[[i]]@meta.data$',heat_ident,sep = '')))
  Heatmap_15987 = DoHeatmap(tempt_diff[[i]], features = top10$gene) + NoLegend()+scale_fill_gradientn(colors = c("navy","white","firebrick3"))
  ggsave('heatmap.png',Heatmap_15987,width=20,height=15)
  ggsave('Volcano.png',EnhancedVolcano(df, lab = df$gene, x = 'avg_log2FC' , y = 'p_val' , xlim = c(-max(abs(df$avg_log2FC))*1.2,max(abs(df$avg_log2FC))*1.2) , ylim = c(0,1.2*max(-log10(ifelse(df$p_val_adj == 0, 1e-300, df$p_val_adj)))) , pCutoff = 0.000001 , FCcutoff = cutoff),width=20,height=15)
  ggsave('DotPlot.png',DotPlot(tempt_diff[[i]], features =unique(top10$gene),cols = c('navy','red'))  + theme(axis.text.x = element_text(face = "italic",angle=90, hjust=1)),width=15,height=15)
  ggsave('heatmap.pdf',Heatmap_15987,width=20,height=15)
  ggsave('Volcano.pdf',EnhancedVolcano(df, lab = df$gene, x = 'avg_log2FC' , y = 'p_val' , xlim = c(-max(abs(df$avg_log2FC))*1.2,max(abs(df$avg_log2FC))*1.2) , ylim = c(0,1.2*max(-log10(ifelse(df$p_val_adj == 0, 1e-300, df$p_val_adj)))) , pCutoff = 0.000001 , FCcutoff = cutoff),width=20,height=15)
  ggsave('DotPlot.pdf',DotPlot(tempt_diff[[i]], features =unique(top10$gene),cols = c('navy','red'))  + theme(axis.text.x = element_text(face = "italic",angle=90, hjust=1)),width=15,height=15)
  
  write.csv(markers , "diff_gene.csv")
  #markers$gene = rownames(markers)
  Focus=rownames( markers [intersect(intersect(which(markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=cutoff)),which( as.character(markers [,"cluster"])==disease.type)),])
  Control=rownames( markers [intersect(intersect(which(markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=cutoff)),which( as.character(markers [,"cluster"])==control.type)),])
  Focus_no_only.pos=rownames( markers [intersect(which(markers [,"p_val"]<0.05),which( as.character(markers [,"cluster"])==disease.type)),])
  Control_no_only.pos=rownames( markers [intersect(which(markers [,"p_val"]<0.05),which( as.character(markers [,"cluster"])==control.type)),])
  
  
  dir.create(disease.type)
  setwd(paste('./',disease.type,sep=''))
  x = markers[Focus,'gene']
  tmp=gsub("\\..*","",x)
  eg=NULL
  #try(eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1))
  eg = tryCatch({eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)},error = function(e){})
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  if (!is.null(go)){
  write.csv(go@result,'go_result.csv')
  go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  try(ggsave("1.1.1.GO_barplot.pdf",barplot(go, split="ONTOLOGY",showCategory = 15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
  try(ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
  kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  try(ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
  try(ggsave("1.1.1.GO_barplot.png",barplot(go, split="ONTOLOGY",showCategory=15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
  try(ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
  if(CN){  try(ggsave("1.1.1.GO_barplot_ZH.png",translate_ggplot(barplot(go, split="ONTOLOGY",showCategory = 15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),axis = 'y',from = 'en',to='zh'),width=14,height=10))
    try(ggsave("2.1.KEGG_dotplot_ZH.png",translate_ggplot(dotplot(kegg, showCategory=30,label_format=100),axis = 'y',from = 'en',to='zh'),width=14,height=10))
    }  
  try(ggsave('3.1GO_gene_Connection.png',enrichplot::cnetplot(go,circular=FALSE,colorEdge = TRUE),width=14,height=10))#基因-通路关联网络图
  try(ggsave('3.2KEGG_gene_Connection.png',enrichplot::cnetplot(kegg,circular=FALSE,colorEdge = TRUE),width=14,height=10))#circluar为指定是否环化，基因过多时建议设置为FALSE
  try(GO2 <- pairwise_termsim(go))
  try(KEGG2 <- pairwise_termsim(kegg))
  try(ggsave('4.1GO_functions_Connection.png',enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk"),width=25,height=25))
  try(ggsave('4.2KEGG_functions_Connection.png',enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk"),width=25,height=25))
  
  # x = markers[Focus_no_only.pos,]
  # tmp=gsub("\\..*","",x$gene)
  # eg=NULL
  # eg = tryCatch({eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)},error = function(e){})
  # genelist <- eg$ENTREZID
  # go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # if (!is.null(go)){
  # go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # 
  # genedata<-data.frame(ID=tmp,logFC=x$avg_log2FC)
  # BP.row=rownames(go[which(go[,"ONTOLOGY"]=='BP'),])
  # CC.row=rownames(go[which(go[,"ONTOLOGY"]=='CC'),])
  # MF.row=rownames(go[which(go[,"ONTOLOGY"]=='MF'),])
  # 
  # GOplotIn_BP<-go[1:10,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
  # GOplotIn_CC<-go[length(BP.row)+1:min(length(BP.row)+length(CC.row),length(BP.row)+10),c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
  # GOplotIn_MF<-go[length(BP.row)+length(CC.row)+1:min(length(BP.row)+length(CC.row)+length(MF.row),length(BP.row)+length(CC.row)+10),c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
  # GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
  # GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
  # GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
  # names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
  # names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
  # names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
  # GOplotIn_BP$Category = "BP"#分类信息
  # GOplotIn_CC$Category = "CC"
  # GOplotIn_MF$Category = "MF"
  # try(circ_BP<-circle_dat(GOplotIn_BP,genedata) )#GOplot导入数据格式整理
  # try(circ_CC<-circle_dat(GOplotIn_CC,genedata) )
  # try(circ_MF<-circle_dat(GOplotIn_MF,genedata) )
  # try(chord_BP<-chord_dat(data = circ_BP,genes = genedata)) #生成含有选定基因的数据框
  # try(chord_CC<-chord_dat(data = circ_CC,genes = genedata)) 
  # try(chord_MF<-chord_dat(data = circ_MF,genes = genedata)) 
  # try(ggsave('GOChord_BP.png',GOChord(data = chord_BP,#弦图
  #                                     title = 'GO-Biological Process',space = 0.01,
  #                                     limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
  #                                     lfc.col = c('red','white','blue'), 
  #                                     process.label = 8),width=25,height=25))
  # try(ggsave('GOChord_CC.png',GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
  #                                     limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
  #                                     lfc.col = c('red','white','blue'), 
  #                                     process.label = 8) ,width=25,height=25))
  # try(ggsave('GOChord_MF.png',GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
  #                                     limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
  #                                     lfc.col = c('red','white','blue'), 
  #                                     process.label = 8),width=25,height=25))
  # try(ggsave('GOCircle_BP.png',GOCircle(circ_BP),width=15,height=15))
  # try(ggsave('GOCircle_CC.png',GOCircle(circ_CC),width=15,height=15))
  # try(ggsave('GOCircle_MF.png',GOCircle(circ_MF) ,width=15,height=15))
  # }
  }
  setwd('../')
  
  dir.create(control.type)
  setwd(paste('./',control.type,sep=''))
  x = markers[Control,'gene']
  tmp=gsub("\\..*","",x)
  eg=NULL
  eg = tryCatch({eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)},error = function(e){})
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  if (!is.null(go)){
  write.csv(go@result,'go_result.csv')
  go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  try(ggsave("1.1.1.GO_barplot.pdf",barplot(go, split="ONTOLOGY",showCategory = 15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
  try(ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
  kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  try(ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
  try(ggsave("1.1.1.GO_barplot.png",barplot(go, split="ONTOLOGY",showCategory=15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
  try(ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
  try(ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
  if(CN){  try(ggsave("1.1.1.GO_barplot_ZH.png",translate_ggplot(barplot(go, split="ONTOLOGY",showCategory = 15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),axis = 'y',from = 'en',to='zh'),width=14,height=10))
    try(ggsave("2.1.KEGG_dotplot_ZH.png",translate_ggplot(dotplot(kegg, showCategory=30,label_format=100),axis = 'y',from = 'en',to='zh'),width=14,height=10))
      }
  try(ggsave('3.1GO_gene_Connection.png',enrichplot::cnetplot(go,circular=FALSE,colorEdge = TRUE),width=14,height=10))#基因-通路关联网络图
  try(ggsave('3.2KEGG_gene_Connection.png',enrichplot::cnetplot(kegg,circular=FALSE,colorEdge = TRUE),width=14,height=10))#circluar为指定是否环化，基因过多时建议设置为FALSE
  try(GO2 <- pairwise_termsim(go))
  try(KEGG2 <- pairwise_termsim(kegg))
  try(ggsave('4.1GO_functions_Connection.png',enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk"),width=25,height=25))
  try(ggsave('4.2KEGG_functions_Connection.png',enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk"),width=25,height=25))
  
  # x = markers[Control_no_only.pos,]
  # tmp=gsub("\\..*","",x$gene)
  # eg=NULL
  # eg = tryCatch({eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)},error = function(e){})
  # genelist <- eg$ENTREZID
  # go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # if (!is.null(go)){
  # go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  # 
  # genedata<-data.frame(ID=tmp,logFC=x$avg_log2FC)
  # BP.row=rownames(go[which(go[,"ONTOLOGY"]=='BP'),])
  # CC.row=rownames(go[which(go[,"ONTOLOGY"]=='CC'),])
  # MF.row=rownames(go[which(go[,"ONTOLOGY"]=='MF'),])
  # 
  # GOplotIn_BP<-go[1:10,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
  # GOplotIn_CC<-go[length(BP.row)+1:min(length(BP.row)+length(CC.row),length(BP.row)+10),c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
  # GOplotIn_MF<-go[length(BP.row)+length(CC.row)+1:min(length(BP.row)+length(CC.row)+length(MF.row),length(BP.row)+length(CC.row)+10),c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
  # GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
  # GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
  # GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
  # names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
  # names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
  # names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
  # GOplotIn_BP$Category = "BP"#分类信息
  # GOplotIn_CC$Category = "CC"
  # GOplotIn_MF$Category = "MF"
  # try(circ_BP<-circle_dat(GOplotIn_BP,genedata) )#GOplot导入数据格式整理
  # try(circ_CC<-circle_dat(GOplotIn_CC,genedata) )
  # try(circ_MF<-circle_dat(GOplotIn_MF,genedata) )
  # try(chord_BP<-chord_dat(data = circ_BP,genes = genedata)) #生成含有选定基因的数据框
  # try(chord_CC<-chord_dat(data = circ_CC,genes = genedata)) 
  # try(chord_MF<-chord_dat(data = circ_MF,genes = genedata)) 
  # try(ggsave('GOChord_BP.png',GOChord(data = chord_BP,#弦图
  #                                     title = 'GO-Biological Process',space = 0.01,
  #                                     limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
  #                                     lfc.col = c('red','white','blue'), 
  #                                     process.label = 8),width=25,height=25))
  # try(ggsave('GOChord_CC.png',GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
  #                                     limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
  #                                     lfc.col = c('red','white','blue'), 
  #                                     process.label = 8) ,width=25,height=25))
  # try(ggsave('GOChord_MF.png',GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
  #                                     limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
  #                                     lfc.col = c('red','white','blue'), 
  #                                     process.label = 8),width=25,height=25))
  # try(ggsave('GOCircle_BP.png',GOCircle(circ_BP),width=15,height=15))
  # try(ggsave('GOCircle_CC.png',GOCircle(circ_CC),width=15,height=15))
  # try(ggsave('GOCircle_MF.png',GOCircle(circ_MF) ,width=15,height=15))
  # }
  }
  setwd('../')
  
  setwd('../')
}

setwd('../')






#multiGroup GO enrichment#######################
dir.create('MultiGroup_Diff')
setwd('./MultiGroup_Diff')

ABC@meta.data$test = paste(ABC@meta.data$Group,ABC@meta.data$manusubcelltype,sep='__')
ABC@meta.data$test = factor(ABC@meta.data$test,levels = sort(levels(as.factor(ABC@meta.data$test))))
Idents(ABC) = ABC@meta.data$test
df <- FindAllMarkers(ABC, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25)
#write.csv(df , "marker.csv")

df_sig  <- df[df$p_val_adj < 0.05, ]

group <- data.frame(gene=df_sig$gene,
                    group=df_sig$cluster)

Gene_ID <- bitr(df_sig$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

#构建文件并分析
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')

data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

data_GO_sim <- simplify(data_GO, 
                        cutoff=0.7, 
                        by="p.adjust", 
                        select_fun=min)


dotplot(data_GO_sim, showCategory=5,font.size = 8,label_format = 100)
translate_ggplot(dotplot(data_GO_sim, showCategory=5,font.size = 8,label_format = 100),axis = 'y',from = 'en',to='zh')
data_GO_sim_fil <- data_GO_sim@compareClusterResult

#导出数据，在ggplot2中可视化,每组中挑选需要可视化的terms进行可视化；
df_GO = data_GO
#df_GO <- read.csv("df_GO.csv", header = T)
library(forcats)
df_GO$Description <- as.factor(df_GO$Description)
df_GO$Description <- fct_inorder(df_GO$Description)

ggplot(df_GO, aes(Cluster, Description)) +
  geom_point(aes(fill=p.adjust, size=Count), shape=21)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text = element_text(color = 'black', size = 10))+
  scale_fill_gradient(low="purple",high="yellow")+
  labs(x=NULL,y=NULL)+
  coord_flip()

#按细胞类型分别比较每个患者的实验组对照组差异表达基因富集以及整合后泽差异基因top10表达###########################################


dir.create('CellType_Diff_patient_top10')
setwd('./CellType_Diff_patient_top10')
ABC = ABC

DefaultAssay(ABC) = 'RNA'
ABC <- NormalizeData(ABC)
ABC <- FindVariableFeatures(ABC, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(ABC)
ABC <- ScaleData(ABC, features = all.genes)

patient_by = 'Patient_No.'
divide_by = 'manucelltype'
Idents(ABC) = ABC@meta.data$Group
disease.type = 'FCD'
control.type = 'noFCD'
cutoff = 0.5
OrgDb1="org.Hs.eg.db"   #鼠：org.Mm.eg.db 人：org.Hs.eg.db

eval(parse(text = paste('tempt_diff2 = SplitObject(ABC,split.by = "',divide_by,'")',sep = '')))

for (j in 1:length(tempt_diff2)){
  dir.create(names(tempt_diff2[j]))
  setwd(paste('./',names(tempt_diff2[j]),sep = ''))
  eval(parse(text = paste('tempt_diff = SplitObject(tempt_diff2[[j]],split.by = "',patient_by,'")',sep = '')))
  top10marker = NULL
  top10marker = list()
  for (i in 1:length(tempt_diff)) {
  dir.create(names(tempt_diff[i]))
  setwd(paste('./',names(tempt_diff[i]),sep = ''))
    

  flag = TRUE
  Idents(tempt_diff[[i]]) = tempt_diff[[i]]@meta.data$Group
  tempt_diff[[i]] <- NormalizeData(tempt_diff[[i]])
  tempt_diff[[i]] <- FindVariableFeatures(tempt_diff[[i]], selection.method = "vst", nfeatures = 4000)
  all.genes <- rownames(tempt_diff[[i]])
  tempt_diff[[i]] <- ScaleData(tempt_diff[[i]], features = all.genes)
  tryCatch({
    markers <- FindAllMarkers(tempt_diff[[i]], only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  },warning = function(w){
    flag<<-FALSE
    print("warning")
  },error = function(e){
    flag<<-FALSE
    print("error")
  })
  if(!flag){
    setwd('../')
    flag = TRUE
    next
  }
  markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10
  top10marker = c(top10marker,top10$gene)
  Idents(tempt_diff[[i]]) = tempt_diff[[i]]@meta.data$orig.ident
  Heatmap_15987 = DoHeatmap(tempt_diff[[i]], features = top10$gene) + NoLegend()
  ggsave('heatmap.png',Heatmap_15987,width=20,height=15)
  ggsave('DotPlot.png',DotPlot(tempt_diff[[i]], features =unique(top10$gene),cols = c('blue','red'))  + theme(axis.text.x = element_text(face = "italic",angle=90, hjust=1)),width=15,height=15)
  write.csv(markers , "diff_gene.csv")
  #markers$gene = rownames(markers)
  Focus=rownames( markers [intersect(intersect(which(markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=cutoff)),which( as.character(markers [,"cluster"])==disease.type)),])
  Control=rownames( markers [intersect(intersect(which(markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=cutoff)),which( as.character(markers [,"cluster"])==control.type)),])
  Focus_no_only.pos=rownames( markers [intersect(which(markers [,"p_val"]<0.05),which( as.character(markers [,"cluster"])==disease.type)),])
  Control_no_only.pos=rownames( markers [intersect(which(markers [,"p_val"]<0.05),which( as.character(markers [,"cluster"])==control.type)),])
  
  {
  dir.create(disease.type)
  setwd(paste('./',disease.type,sep=''))
  x = markers[Focus,'gene']
  tmp=gsub("\\..*","",x)
  eg=NULL
  #try(eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1))
  eg = tryCatch({eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)},error = function(e){})
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  if (!is.null(go)){
    write.csv(go@result,'go_result.csv')
    go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    try(ggsave("1.1.1.GO_barplot.pdf",barplot(go, split="ONTOLOGY",showCategory = 15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
    try(ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
    try(ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
    try(ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
    kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
    try(ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
    try(ggsave("1.1.1.GO_barplot.png",barplot(go, split="ONTOLOGY",showCategory=15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
    try(ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
    try(ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
    try(ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
    try(ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
    
    try(ggsave('3.1GO_gene_Connection.png',enrichplot::cnetplot(go,circular=FALSE,colorEdge = TRUE),width=14,height=10))#基因-通路关联网络图
    try(ggsave('3.2KEGG_gene_Connection.png',enrichplot::cnetplot(kegg,circular=FALSE,colorEdge = TRUE),width=14,height=10))#circluar为指定是否环化，基因过多时建议设置为FALSE
    try(GO2 <- pairwise_termsim(go))
    try(KEGG2 <- pairwise_termsim(kegg))
    try(ggsave('4.1GO_functions_Connection.png',enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk"),width=25,height=25))
    try(ggsave('4.2KEGG_functions_Connection.png',enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk"),width=25,height=25))
    
    x = markers[Focus_no_only.pos,]
    tmp=gsub("\\..*","",x$gene)
    eg=NULL
    eg = tryCatch({eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)},error = function(e){})
    genelist <- eg$ENTREZID
    go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    if (!is.null(go)){
      go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      
      genedata<-data.frame(ID=tmp,logFC=x$avg_log2FC)
      BP.row=rownames(go[which(go[,"ONTOLOGY"]=='BP'),])
      CC.row=rownames(go[which(go[,"ONTOLOGY"]=='CC'),])
      MF.row=rownames(go[which(go[,"ONTOLOGY"]=='MF'),])
      
      GOplotIn_BP<-go[1:10,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
      GOplotIn_CC<-go[length(BP.row)+1:min(length(BP.row)+length(CC.row),length(BP.row)+10),c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
      GOplotIn_MF<-go[length(BP.row)+length(CC.row)+1:min(length(BP.row)+length(CC.row)+length(MF.row),length(BP.row)+length(CC.row)+10),c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
      GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
      GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
      GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
      names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
      names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
      names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
      GOplotIn_BP$Category = "BP"#分类信息
      GOplotIn_CC$Category = "CC"
      GOplotIn_MF$Category = "MF"
      try(circ_BP<-circle_dat(GOplotIn_BP,genedata) )#GOplot导入数据格式整理
      try(circ_CC<-circle_dat(GOplotIn_CC,genedata) )
      try(circ_MF<-circle_dat(GOplotIn_MF,genedata) )
      try(chord_BP<-chord_dat(data = circ_BP,genes = genedata)) #生成含有选定基因的数据框
      try(chord_CC<-chord_dat(data = circ_CC,genes = genedata)) 
      try(chord_MF<-chord_dat(data = circ_MF,genes = genedata)) 
      try(ggsave('GOChord_BP.png',GOChord(data = chord_BP,#弦图
                                          title = 'GO-Biological Process',space = 0.01,
                                          limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
                                          lfc.col = c('red','white','blue'), 
                                          process.label = 8),width=25,height=25))
      try(ggsave('GOChord_CC.png',GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
                                          limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
                                          lfc.col = c('red','white','blue'), 
                                          process.label = 8) ,width=25,height=25))
      try(ggsave('GOChord_MF.png',GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
                                          limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
                                          lfc.col = c('red','white','blue'), 
                                          process.label = 8),width=25,height=25))
      try(ggsave('GOCircle_BP.png',GOCircle(circ_BP),width=15,height=15))
      try(ggsave('GOCircle_CC.png',GOCircle(circ_CC),width=15,height=15))
      try(ggsave('GOCircle_MF.png',GOCircle(circ_MF) ,width=15,height=15))
    }
  }
  setwd('../')
  }
  {
  dir.create(control.type)
  setwd(paste('./',control.type,sep=''))
  x = markers[Control,'gene']
  tmp=gsub("\\..*","",x)
  eg=NULL
  eg = tryCatch({eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)},error = function(e){})
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  if (!is.null(go)){
    write.csv(go@result,'go_result.csv')
    go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    try(ggsave("1.1.1.GO_barplot.pdf",barplot(go, split="ONTOLOGY",showCategory = 15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
    try(ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
    try(ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
    try(ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
    kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
    try(ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
    try(ggsave("1.1.1.GO_barplot.png",barplot(go, split="ONTOLOGY",showCategory=15,label_format=100,drop=T)+facet_grid(ONTOLOGY~., scale="free"),width=14,height=10))
    try(ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,label_format=100,drop=T),width=14,height=10))
    try(ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,label_format=100,drop=T),width=14,height=10))
    try(ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,label_format=100,drop=T),width=14,height=10))
    try(ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30,label_format=100),width=14,height=10))
    
    try(ggsave('3.1GO_gene_Connection.png',enrichplot::cnetplot(go,circular=FALSE,colorEdge = TRUE),width=14,height=10))#基因-通路关联网络图
    try(ggsave('3.2KEGG_gene_Connection.png',enrichplot::cnetplot(kegg,circular=FALSE,colorEdge = TRUE),width=14,height=10))#circluar为指定是否环化，基因过多时建议设置为FALSE
    try(GO2 <- pairwise_termsim(go))
    try(KEGG2 <- pairwise_termsim(kegg))
    try(ggsave('4.1GO_functions_Connection.png',enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk"),width=25,height=25))
    try(ggsave('4.2KEGG_functions_Connection.png',enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk"),width=25,height=25))
    
    x = markers[Control_no_only.pos,]
    tmp=gsub("\\..*","",x$gene)
    eg=NULL
    eg = tryCatch({eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)},error = function(e){})
    genelist <- eg$ENTREZID
    go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
    if (!is.null(go)){
      go.BP <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      go.CC <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      go.MF <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      Go.BP <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      Go.CC <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      Go.MF <- enrichGO(genelist,readable = T, OrgDb = OrgDb1, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
      
      genedata<-data.frame(ID=tmp,logFC=x$avg_log2FC)
      BP.row=rownames(go[which(go[,"ONTOLOGY"]=='BP'),])
      CC.row=rownames(go[which(go[,"ONTOLOGY"]=='CC'),])
      MF.row=rownames(go[which(go[,"ONTOLOGY"]=='MF'),])
      
      GOplotIn_BP<-go[1:10,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
      GOplotIn_CC<-go[length(BP.row)+1:min(length(BP.row)+length(CC.row),length(BP.row)+10),c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
      GOplotIn_MF<-go[length(BP.row)+length(CC.row)+1:min(length(BP.row)+length(CC.row)+length(MF.row),length(BP.row)+length(CC.row)+10),c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
      GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
      GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
      GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
      names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
      names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
      names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
      GOplotIn_BP$Category = "BP"#分类信息
      GOplotIn_CC$Category = "CC"
      GOplotIn_MF$Category = "MF"
      try(circ_BP<-circle_dat(GOplotIn_BP,genedata) )#GOplot导入数据格式整理
      try(circ_CC<-circle_dat(GOplotIn_CC,genedata) )
      try(circ_MF<-circle_dat(GOplotIn_MF,genedata) )
      try(chord_BP<-chord_dat(data = circ_BP,genes = genedata)) #生成含有选定基因的数据框
      try(chord_CC<-chord_dat(data = circ_CC,genes = genedata)) 
      try(chord_MF<-chord_dat(data = circ_MF,genes = genedata)) 
      try(ggsave('GOChord_BP.png',GOChord(data = chord_BP,#弦图
                                          title = 'GO-Biological Process',space = 0.01,
                                          limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
                                          lfc.col = c('red','white','blue'), 
                                          process.label = 8),width=25,height=25))
      try(ggsave('GOChord_CC.png',GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
                                          limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
                                          lfc.col = c('red','white','blue'), 
                                          process.label = 8) ,width=25,height=25))
      try(ggsave('GOChord_MF.png',GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
                                          limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
                                          lfc.col = c('red','white','blue'), 
                                          process.label = 8),width=25,height=25))
      try(ggsave('GOCircle_BP.png',GOCircle(circ_BP),width=15,height=15))
      try(ggsave('GOCircle_CC.png',GOCircle(circ_CC),width=15,height=15))
      try(ggsave('GOCircle_MF.png',GOCircle(circ_MF) ,width=15,height=15))
    }
  }
  setwd('../')
  }
  setwd('../')
  }
  
  tempt_diff2[[j]] <- NormalizeData(tempt_diff2[[j]])
  tempt_diff2[[j]] <- FindVariableFeatures(tempt_diff2[[j]], selection.method = "vst", nfeatures = 4000)
  all.genes <- rownames(tempt_diff2[[j]])
  tempt_diff2[[j]] <- ScaleData(tempt_diff2[[j]], features = all.genes)
  Idents(tempt_diff2[[j]]) = tempt_diff2[[j]]@meta.data$orig.ident
  Heatmap_15987 = DoHeatmap(tempt_diff2[[j]], features = as.character(top10marker)) + NoLegend()
  ggsave('heatmap.png',Heatmap_15987,width=20,height=length(tempt_diff2)*2)
  ggsave('DotPlot.png',DotPlot(tempt_diff2[[j]], features =unique(as.character(top10marker)),cols = c('blue','yellow'))  + theme(axis.text.x = element_text(face = "italic",angle=90, hjust=1)),width=length(tempt_diff2)*4,height=15)
  
  setwd('../')
}

setwd('../')









##########################################

pdf("1.2.1.GO_BP_barplot.pdf",width=14,height=10)
print(barplot(go.BP,showCategory=20,drop=T))
dev.off()

pdf("1.2.2.GO_BP_dotplot.pdf",width=14,height=10)
print(dotplot(go.BP,showCategory=50))
dev.off()

pdf("1.2.3.GO_BP_plotGOgraph.pdf",width=14,height=10)
print(plotGOgraph(Go.BP))
dev.off()
##############################################

pdf("1.3.1.GO_CC_barplot.pdf",width=14,height=10)
print(barplot(go.CC,showCategory=20,drop=T))
dev.off()

pdf("1.3.2.GO_CC_dotplot.pdf",width=14,height=10)
print(dotplot(go.CC,showCategory=50))
dev.off()

pdf("1.3.3.GO_CC_plotGOgraph.pdf",width=14,height=10)
print(plotGOgraph(Go.CC))
dev.off()
##############################################

pdf("1.4.1.GO_MF_barplot.pdf",width=14,height=10)
print(barplot(go.MF,showCategory=20,drop=T))
dev.off()

pdf("1.4.2.GO_MF_dotplot.pdf",width=14,height=10)
print(dotplot(go.MF,showCategory=50))
dev.off()

pdf("1.4.3.GO_MF_plotGOgraph.pdf",width=14,height=10)
print(plotGOgraph(Go.MF))
dev.off()
##############################################

kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
pdf("2.1.KEGG_dotplot.pdf",width=14,height=10)
print(dotplot(kegg, showCategory=30))
dev.off()

##############################################
degenes = merge(x,eg,by.x="gene",by.y=1)
geneList <- degenes[,c("gene","ENTREZID","ENSEMBL","avg_log2FC")]
geneList.sort <- arrange(geneList, desc(avg_log2FC)) 
gene = geneList.sort$ENTREZID

gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
enrich <- enricher(gene, TERM2GENE=c5) 

glist <- geneList$avg_log2FC     
names(glist) <- as.character(geneList$ENTREZID)
glist <- sort(glist,decreasing = T)

gsea <- GSEA(glist, TERM2GENE=c5, verbose=FALSE, pvalueCutoff = 0.5); head(gsea)


pdf("3.1.GSEA_1_plot.pdf",width=14,height=10)
print(gseaplot2(gsea,1))
dev.off()

pdf("3.2.GSEA_3_plot.pdf",width=14,height=10)
print(gseaplot2(gsea,1:3))
dev.off()
#################################################

#markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25, logfc.threshold = 0)  ##ͬ??ѡȡC5??C0??C3?رȽϣ????????????????????????л?????ǰ?˿?????ΪC5?ϵ???????ΪC5?µ?????
#markers$genes = rownames(markers)
markers.sub = subset(markers,idents = )
cluster.genes<- markers.sub %>% arrange(desc(avg_logFC)) %>% dplyr::select(genes,avg_logFC) #??????logFC????
ranks<- deframe(cluster.genes)

mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")## ???????򼯣?ѡȡC2
fgsea_sets = mdb_c2 [grep("^KEGG",mdb_c2$gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000) #????fgsea

p <-ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES)) +
  coord_flip() +
  labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA") ##????????????ǰ20????Ŀ
pdf('GSEA-fgsea.pdf',width=8,height=5)
print(p)
dev.off()
pdf('fgsea_KEGG_PRIMARY_IMMUNODEFICIENCY.pdf',width=8,height=5)
plotEnrichment(fgsea_sets[["KEGG_PRIMARY_IMMUNODEFICIENCY"]],ranks) + labs(title="KEGG_PRIMARY_IMMUNODEFICIENCY") #??ĳһ?ض?ͨ·????
dev.off()


##################################################
setwd("")
pbmc <-readRDS("pbmc.rds")
expr <- as.data.frame(pbmc@assays$RNA@data)             #????????
meta <- pbmc@meta.data[,c("orig.ident","seurat_clusters")]           #????
m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")           #ѡȡ????????
#
#H??hallmark gene sets---??֢???????򼯺ϣ???50 gene sets??
#C1??positional gene sets---Ⱦɫ??λ?û??򼯺ϣ???299 gene sets??
#C2??curated gene sets---????????֪???ݿ⣬???׵Ȼ???????Ϣ??????5529 gene sets??
#C3??motif gene sets---???ذл??򼯺ϣ???��miRNA-target?Լ?ת¼????-target????ģʽ??3735 gene sets??
#C4??computational gene sets---??????????Ԥ????��?Ļ??򼯺ϣ???Ҫ?ǺͰ?֢???صĻ?????858 gene sets??
#C5??GO gene sets---Gene Ontology??Ӧ?Ļ??򼯺ϣ?10192 gene sets??
#C6??oncogenic signatures---?°????򼯺ϣ??󲿷?��Դ??NCBI GEO????оƬ???ݣ?189 gene sets??
#C7??immunologic signatures---???????ػ??򼯣?4872 gene sets??
#C8??single cell identitiy gene sets, 302 gene sets
#



