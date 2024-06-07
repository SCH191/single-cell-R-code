library(Seurat)
#library(SeuratData)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat) 
library(SingleR)
library(ggplot2)
library(reshape2)
library(showtext)
hpca.se=celldex::HumanPrimaryCellAtlasData() 
Blue.se=BlueprintEncodeData() 
Immune.se=DatabaseImmuneCellExpressionData()
Nover.se=NovershternHematopoieticData()
MonacoIm.se=MonacoImmuneData()
ImmGen.se=ImmGenData() 
Mouse.se=MouseRNAseqData() 

setwd("G:/")
font_add('YaHei Consolas', regular = 'G:\\YaHei_Consolas_Hybrid_1.12.ttf')
showtext_auto()
setwd("./kongjian/data2")

all_data=c("FCD_1","FCD2_1","FCD3_8","FCD2_CON3_1")
i=all_data[2]

ABC.data = Read10X(data.dir=paste("G:/filterd data/4/",i,"/2.2.filtered_feature_bc_matrix",sep=""))
ABC <- CreateSeuratObject(counts = ABC.data, project = i, min.cells = 3, min.features = 200)
img=Seurat::Read10X_Image(image.dir=paste("G:/filterd data/4/",i,"/spatial",sep=""))
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = ABC)]
ABC[['image']] <- img


#plot1 <- VlnPlot(ABC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
#plot2 <- SpatialFeaturePlot(ABC, features = "nCount_Spatial") + theme(legend.position = "right")
#plot3 <- VlnPlot(ABC, features ="nFeature_Spatial",pt.size = 0.1,cols ="tomato") + NoLegend()
#plot4 <- FeatureScatter(ABC, feature1 ="nCount_Spatial",feature2 ="nFeature_Spatial")+ NoLegend()
#plot1/plot3 | plot2/plot4

ABC <- SCTransform(ABC, assay = "RNA", return.only.var.genes = FALSE, verbose = FALSE)

ABC <- RunPCA(ABC, assay = "SCT", verbose = FALSE)
ABC <- FindNeighbors(ABC, reduction = "pca", dims = 1:50)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.5)
ABC <- RunUMAP(ABC, reduction = "pca",neighbors = 10,n.epochs=50,negative.sample.rate = 25,min.dist=0, dims = 1:50)
ABC <- RunTSNE(object=ABC, reduction = "pca" , dims=1:50)

p1 <- DimPlot(ABC, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(ABC, label = TRUE, label.size = 3,pt.size.factor = 2)
p3=TSNEPlot(object=ABC , label=TRUE,repel=T, group.by = c("ident"))
p1/p3|p2

pdf(paste(i,"UMAP&tSNE&image.pdf",sep=""),width=14,height=12)
p1/p3|p2
dev.off()

DefaultAssay(ABC) = 'RNA'
markers <- FindAllMarkers(ABC, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers , paste(i,"seurat_clusters_all_marker.csv",sep='_'))

dir.create(paste(i,"GOenrichment",sep=""))
setwd(paste("./",paste(i,"GOenrichment",sep=""),sep = ''))
OrgDb1="org.Hs.eg.db"   #鼠：org.Mm.eg.db 人：org.Hs.eg.db

for (i in 1:length(levels(ABC))) {
  dir.create(as.character(levels(ABC)[i]))
  setwd(paste("./",as.character(levels(ABC)[i]),sep = ''))
  Focus=rownames( markers [intersect(intersect(which(markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=0.5)),which( as.character(markers [,"cluster"])==as.character(levels(ABC)[i]))),])
  Focus_no_only.pos=rownames( markers [intersect(which(markers [,"p_val"]<0.05),which( as.character(markers [,"cluster"])==as.character(levels(ABC)[i]))),])
  
  x = markers[Focus,'gene']
  tmp=gsub("\\..*","",x)
  eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,readable = T,OrgDb=OrgDb1,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
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
}
setwd('../')

################################################
meta=ABC@meta.data 
ABC_for_SingleR <- GetAssayData(ABC, slot="data")
ABC.hesc <- SingleR(test = ABC_for_SingleR, ref = hpca.se, labels = hpca.se$label.main)
ABC@meta.data$celltype <-ABC.hesc$labels

saveRDS(ABC,paste(i,".rds",sep=""))

j=ABC
name=i

pdf(paste(name,".pdf",sep=""),width=12,height=12)
DimPlot(j, group.by = c("seurat_clusters"),reduction = "umap",label=T,repel=T)+ NoLegend()
DimPlot(j, group.by = c( "celltype"),reduction = "umap",label=T,repel=T)+ NoLegend()
DimPlot(j, group.by = c("seurat_clusters"),reduction = "tsne",label=T,repel=T)+ NoLegend()
DimPlot(j, group.by = c("celltype"),reduction = "tsne",label=T,repel=T)+ NoLegend()
SpatialDimPlot(j, label = TRUE, label.size = 3,pt.size.factor = 2)

marker=c("MOG","OPALIN","KLK6","GPR37","ANLN")
plot(c(1:10),family="YaHei Consolas",main="??ͻ????ϸ?? (Oligo)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("PDGFRA","GPR17","CCND1","TNR","NEU4")
plot(c(1:10),family="YaHei Consolas",main="??ͻ????ϸ??ǰ?? (OPC)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("SLC4A4","NTSR2","GFAP","ACSBG1","ALDOC")         #(KCNJ10?????ο?)
plot(c(1:10),family="YaHei Consolas",main="???ν???ϸ?? (Astrocyte)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("C1QC","CD68","C1QB","TREM2","CCL4")
plot(c(1:10),family="YaHei Consolas",main="С????ϸ?? (Microglia)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("KLRD1","NKG7","GZMA","XCL1","GZMK")
plot(c(1:10),family="YaHei Consolas",main="??Ȼɱ??ϸ?? (NK)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("CD79A","CCR7","CD27","IGLC2","CXCR4")
plot(c(1:10),family="YaHei Consolas",main="Bϸ?? (B cell)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("TRBC2","CD3D","CD3G","GZMK","CD3E")
plot(c(1:10),family="YaHei Consolas",main="Tϸ?? (T cell)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("APOBEC3A","TET2","MEFV","HCK","ITGAX","RBPJ")
plot(c(1:10),family="YaHei Consolas",main="????ϸ?? (Monocytes)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("CD66b","CD11b","CD15","CD16","CXCL8","FCGR3B","MNDA")
plot(c(1:10),family="YaHei Consolas",main="????��ϸ?? (Neutrophlis)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("EPCAM","CD24")
plot(c(1:10),family="YaHei Consolas",main="??Ƥϸ?? (Epithelial cell)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("CLDN5","FLT1","PECAM1","VWF","CLEC14A")
plot(c(1:10),family="YaHei Consolas",main="??Ƥϸ?? (Endothelial)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("KCNMB1","MRVI1","KCNAB1","DMPK")
plot(c(1:10),family="YaHei Consolas",main="ƽ????ϸ?? (Smooth muscle cells)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("COL6A2","VIM","COL1A1","PDGFRB")
plot(c(1:10),family="YaHei Consolas",main="????άϸ?? (Fibroblasts)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("TAGLN","ACTA2","PDGFRB","DES")
plot(c(1:10),family="YaHei Consolas",main="??ϸ?? (Mural)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("STMN2","THY1")
plot(c(1:10),family="YaHei Consolas",main="????Ԫϸ?? (Neuron)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("SLC17A7","SATB2")
plot(c(1:10),family="YaHei Consolas",main="?˷???????Ԫϸ?? (Ex Neuron)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("GAD1","GAD2","SLC32A1","SLC6A1")
plot(c(1:10),family="YaHei Consolas",main="??????????Ԫϸ?? (In Neuron)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("CALB1","CCK","EN1","ETV1","EVX1")
plot(c(1:10),family="YaHei Consolas",main="?м?????Ԫ (inter Neuron)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("GAD1","SLC32A1","PVALB","SLC6A1")
plot(c(1:10),family="YaHei Consolas",main="GABA????Ԫ (GABAergic Neuron)")
VlnPlot(j, features =marker,ncol=1)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

dev.off()

ABC.markers <- FindAllMarkers(ABC, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)

pdf(paste(name,"_space_locate.pdf",sep=""),width=15,height=50)
SpatialDimPlot(ABC, cells.highlight = CellsByIdentities(object = ABC, idents = c(0:(length(levels(ABC))-1))), facet.highlight = TRUE, ncol = 3)
dev.off()

pdf(paste(i,"_HeatMap.pdf",sep=""),width=25,height=35)
ABC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(ABC, features = top10$gene) + NoLegend()
dev.off()



write.csv(ABC.markers , paste(i,"_all_marker.csv",sep=""))
###################################################################################