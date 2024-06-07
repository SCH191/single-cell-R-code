library(Seurat)
#library(SeuratData)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat) 
library(SingleR)
library(ggplot2)
library(reshape2)
library(patchwork)
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
setwd("./danxibao/data2")

all_data=c("FCD_16","normal_16")
i=all_data[2]

ABC.data <- Read10X(data.dir = paste("G:/danxibao/data2/",i,sep=""))
#��????Cellranger count Filtered_gene_bc_matrices?ļ???

#????Seurat??????
ABC <- CreateSeuratObject(counts = ABC.data, project = i, min.cells = 3, min.features = 200)

ABC[["percent.mt"]] <- PercentageFeatureSet(ABC, pattern = "^MT-")
ABC <- subset(ABC, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)


ABC <- NormalizeData(ABC)
ABC <- FindVariableFeatures(ABC, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(ABC)
ABC <- ScaleData(ABC, features = all.genes)

ABC <- SCTransform(ABC, assay = "RNA", return.only.var.genes = FALSE, verbose = FALSE)



ABC <- RunPCA(ABC, assay = "SCT", verbose = FALSE)
ABC <- FindNeighbors(ABC, reduction = "pca", dims = 1:50)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 4)
ABC <- RunUMAP(ABC, reduction = "pca", dims = 1:50)
ABC <- RunTSNE(object=ABC, reduction = "pca" , dims=1:50)

DimPlot(ABC, group.by = c("seurat_clusters", "celltype","orig.ident","manucelltype"),reduction = "umap",label=T,repel = T)

p1 <- DimPlot(ABC, reduction = "umap",group.by = c('seurat_clusters','orig.ident','CellType','manucelltype'), label = TRUE)
p3=TSNEPlot(object=ABC , label=TRUE,repel=T, group.by = c('seurat_clusters','orig.ident','CellType','manucelltype'))
p1/p3

meta=ABC@meta.data 
ABC_for_SingleR <- GetAssayData(ABC, slot="data")
ABC.hesc <- SingleR(test = ABC_for_SingleR, ref = Mouse.se, labels = Mouse.se$label.main)
ABC@meta.data$celltype <-ABC.hesc$labels
DimPlot(ABC, group.by = c("seurat_clusters", "celltype","orig.ident","manucelltype"),reduction = "umap",label=T,repel = T)

saveRDS(ABC,paste(i,".rds",sep=""))

ABC.markers <- FindAllMarkers(ABC, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(ABC.markers , "seurat_clusters_all_marker.csv")

ABC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(top10 , "seurat_clusters_all_marker_top10.csv")
DoHeatmap(ABC, features = top10$gene) + NoLegend()
FeaturePlot(ABC, features =top10$gene,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"),ncol=10,reduction = 'umap')
##################
ABC <- NormalizeData(ABC)
ABC <- FindVariableFeatures(ABC, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(ABC)
ABC <- ScaleData(ABC, features = all.genes)
ABC <- RunPCA(ABC, verbose = FALSE)
ABC <- FindNeighbors(ABC, reduction = "pca", dims = 1:50)

ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.2)

group_by = c('seurat_clusters','orig.ident','Group','manucelltype','manusubcelltype')

ABC <- RunUMAP(ABC, reduction = "pca",n.neighbors = 30L,negative.sample.rate = 5 , dims = 1:50,n.epochs = 500,min.dist = 0.3,spread = 1)
DimPlot(ABC,raster=FALSE, reduction = "umap",group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))


group_by = c('seurat_clusters','orig.ident','Group','Method','Patient_No.','manucelltype','manusubcelltype')

ABC <- RunTSNE(object=ABC, reduction = "pca" , dims=1:50)
#ABC@meta.data$manusubcelltype = ABC@meta.data$manucelltype2
DimPlot(ABC,raster=FALSE, reduction = "umap",group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))/DimPlot(ABC,raster=FALSE, reduction = "tsne",group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))

pix = DimPlot(ABC,raster=FALSE, reduction = "umap",group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))/DimPlot(ABC,raster=FALSE, reduction = "tsne",group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))
ggsave('umap1.png',pix,width=80,height=18,limitsize = FALSE)
##harmony####################
library(harmony)
ABC = RunHarmony(ABC , group.by.vars = 'orig.ident',plot_convergence = T)
ABC <- FindNeighbors(ABC, reduction = "harmony", dims = 1:50)
ABC <- FindClusters(ABC, resolution = 0.2,verbose = FALSE)
ABC <- RunUMAP(ABC, reduction = "harmony",n.neighbors = 30L,negative.sample.rate = 5 , dims = 1:50,n.epochs = 500,min.dist = 0.3,spread = 1)
group_by = c('seurat_clusters','orig.ident','Group','manucelltype','manusubcelltype')
DimPlot(ABC,raster=FALSE, reduction = "umap",group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))
ABC <- RunTSNE(object=ABC, reduction = "harmony" , dims=1:50)
#ABC@meta.data$manusubcelltype = ABC@meta.data$manucelltype2
DimPlot(ABC,raster=FALSE, reduction = "umap",group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))/DimPlot(ABC,raster=FALSE, reduction = "tsne",group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))

##########################################
DefaultAssay(ABC) = 'RNA'
ABC <- NormalizeData(ABC)
ABC <- FindVariableFeatures(ABC, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(ABC)
ABC <- ScaleData(ABC, features = all.genes)
ABC.markers <- FindAllMarkers(ABC, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25)
#write.csv(ABC.markers , "seurat_clusters_all_marker.csv")
ABC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(ABC, features = top10$gene) + NoLegend()

##########################################
current.cluster.ids <- c(1,2,3,4,5,13,15,16,19,21,9,14, 6)
new.cluster.ids <- c('Ex','Ex','Ex','Ex','Ex','Ex','Ex','Ex','Ex','Ex','Ex','Ex','In')
ABC@meta.data$manusubcelltype <- plyr::mapvalues(x = ABC@meta.data[,"manucelltype"], from = current.cluster.ids, to = new.cluster.ids)

DimPlot(ABC, group.by = c("seurat_clusters", "manusubcelltype","orig.ident","manucelltype"),reduction = "umap",label=T,repel = T)

tmp=1
merge@meta.data$manusubcelltype = as.character(merge@meta.data$manusubcelltype)
for (i in 1:length(merge@meta.data$manusubcelltype)) {
  if (rownames(merge@meta.data[i,]) == rownames(ABC@meta.data[tmp,])) {
    merge@meta.data$manusubcelltype[i] = as.character(ABC@meta.data$manusubcelltype[tmp])
    tmp = tmp+1
  }
}
merge@meta.data$manusubcelltype = as.factor(merge@meta.data$manusubcelltype)
DimPlot(merge, reduction = "umap",raster=FALSE,group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))/DimPlot(merge, reduction = "tsne",raster=FALSE,group.by = group_by, label = TRUE ,repel = TRUE,ncol = length(group_by))


#clustree###################################
library(clustree)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.025)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.05)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.075)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.1)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.125)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.15)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.175)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.2)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.225)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.25)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.275)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.3)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.325)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.35)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.375)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.40)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.425)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.45)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.475)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.5)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.525)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.55)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.575)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.6)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.625)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.65)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.675)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.7)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.725)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.75)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.775)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.8)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.825)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.85)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.875)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.9)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.925)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.95)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 0.975)
ABC <- FindClusters(ABC, verbose = FALSE, resolution = 1)

clustree(ABC)

DimPlot(ABC,raster=FALSE, reduction = "umap",group.by = c('integrated_snn_res.0.05','integrated_snn_res.0.1','integrated_snn_res.0.15','integrated_snn_res.0.2','integrated_snn_res.0.25','integrated_snn_res.0.3'), label = TRUE ,repel = TRUE,ncol = 3)

#在umap中看nFeature和UMI分布#####################
yd = DimPlot(ABC,pt.size = 0.1)
data = cbind(ABC@meta.data,yd$data)
head(data)
a1=ggplot(data,aes(umap_1,umap_2))+geom_point(aes(color = nFeature_RNA),size = 0.1)+scale_color_gradient(low = "#FFFF33",high = "#3A0088")+theme_bw()

a2=ggplot(data,aes(umap_1,umap_2))+geom_point(aes(color = nCount_RNA),size = 0.1)+scale_color_gradient(low = "#FFFF33",high = "#3A0088")+theme_bw()

a3=ggplot(data,aes(umap_1,umap_2))+geom_point(aes(color = percent.mt),size = 0.1)+scale_color_gradient(low = "#FFFF33",high = "#3A0088")+theme_bw()
a1+a2+a3
#基因表达皮尔森相关性#############################
library(Seurat)
pbmc = ABC
#table(pbmc$cell_type)
av = AverageExpression(pbmc,
                       group.by = "seurat_clusters",
                       assays = "RNA")
av=av[[1]]
#head(av)
#               Naive CD4 T Memory CD4 T  CD14+ Mono          B      CD8 T FCGR3A+ Mono         NK         DC Platelet
# AL627309.1    0.006007987   0.04854338 0.006065400 0.00000000 0.01995673   0.00000000 0.00000000 0.00000000        0
# AP006222.2    0.000000000   0.01088471 0.008397321 0.00000000 0.01157323   0.00000000 0.00000000 0.00000000        0
# RP11-206L10.2 0.007306336   0.00000000 0.000000000 0.02065031 0.00000000   0.00000000 0.00000000 0.08462847        0
# RP11-206L10.9 0.000000000   0.01050116 0.000000000 0.00000000 0.00000000   0.01200008 0.00000000 0.00000000        0
# LINC00115     0.014872162   0.03753737 0.031095957 0.03888541 0.01892413   0.01469374 0.06302423 0.00000000        0
# NOC2L         0.501822970   0.27253750 0.346874253 0.58653489 0.59129394   0.50026775 0.65705381 0.32861136        0

#选出标准差最大的1000个基因
cg=names(tail(sort(apply(av, 1, sd)),5000))
#查看这1000个基因在各细胞群中的表达矩阵
#View(av[cg,])
#查看细胞群的相关性矩阵
#View(cor(av[cg,],method = 'spearman'))
#pheatmap绘制热图
pheatmap::pheatmap(cor(av[cg,],method = 'spearman')) #默认是Pearson

#Umap调整#########################################################
#n_neighbors：控制模糊搜索区域的半径：更少邻域 到 更多邻域;
#min_dist：低维下允许的行间最小距离：更集中 到 更分散；(0~1)
#metric：选择距离的测度方法：欧氏距离、曼哈顿距离等；
#n_epochs：优化步骤的迭代次数。
ABC = RunUMAP(ABC,neighbors = 50,n.epochs=500,negative.sample.rate = 10,min.dist=0.1,dims = 1:50) 
DimPlot(ABC, group.by = c("seurat_clusters", "manusubcelltype","orig.ident","manucelltype"),reduction = "umap",label=T,repel = T)

#聚类树##############################################################################################

ABC <- BuildClusterTree(ABC,assay="RNA",dims = 1:50,reduction = "pca")
library(ggtree)
myPhyTree <- Tool(object=ABC, slot = "BuildClusterTree")
ggtree(myPhyTree)+geom_tiplab()+theme_tree()

#手动框选细胞##################################################################################################
ABC = ABC

# CellSelector
p1=DimPlot(ABC,reduction = 'umap',group.by = 'manusubcelltype',raster=FALSE); p1
# Follow instructions in the terminal to select points
cells.located <- CellSelector(plot = p1)
cells.located #返回细胞id
# [1] "ACGTCGCTCCTGAA-1" "ACTTCAACAAGCAA-1" "CACCGGGACTTCTA-1" "CGATACGACAGGAG-1" "CGCAGGTGGGAACG-1"
# [6] "GCCCATACAGCGTT-1" "GCCTCAACTCTTTG-1" "GTATTAGAAACAGA-1"
#
#pbmc2 = pbmc
#ABC@meta.data$manusubcelltype = NA
Idents(ABC) = ABC@meta.data$manusubcelltype
# Automatically set the identity class of selected cells and return a new Seurat object
ABC <- CellSelector(plot = p1, object = ABC, ident = 'TBD')
ABC@meta.data$manusubcelltype = Idents(ABC)
DimPlot(ABC,group.by = c('manucelltype','manusubcelltype'),label = T,repel = T,reduction = 'umap',raster=FALSE)
#> table(pbmc2@meta.data$seurat_clusters) #没变
#0   1   2   3   4   5   6   7   8 
#684 481 476 344 291 162 155  32  13

#> table(pbmc2@active.ident) #新增一类
#SelectedCells             0             1             2             3             4             5 
#36           681           480           476           344           291           162 
#6             8 
#155            13

#> levels(pbmc2) #新增一类
#[1] "SelectedCells" "0"             "1"             "2"             "3"             "4"            
#[7] "5"             "6"             "8" 


#标注分群细胞区域#################################################

library(scRNAtoolVis)
clusterCornerAxes(object = ABC,reduction = 'umap',
                  noSplit = T,
                  clusterCol = 'seurat_clusters',
                  cornerTextSize = 3.5,
                  themebg = 'bwCorner',
                  addCircle = TRUE,
                  cicAlpha = 0.2,
                  nbin = 200)

#decontX去除游离细胞污染#################################



library(Seurat)
library(decontX)

# 如果你已经有一个处理好的seurat对象，但不想转化数据类型，可以直接使用这个方法
# 读取seurat对象
ABC <- 

# 提取矩阵(genes x cells)
# 需要注意的是，你应该查看你的矩阵是否已经注释了行名和列明，以及如果你的矩阵是scanpy对象中提取的矩阵应该行列转置
counts <- ABC@assays$RNA@counts

# 需要把结果储存在新变量
decontX_results <- decontX(counts) 

# 你可以使用str()查看结果的形式

# RNA污染的计算结果储存在：decontX_results$contamination
# 他是一个数值型向量，长度与细胞数量一致，顺序也与矩阵的colnames一致
# 我们可以直接把他写入metadata
ABC$Contamination =decontX_results$contamination
head(ABC@meta.data) # 可以查看一下结果，这里就不展示了

# 我们对他进行可视化
library(ggplot2)
FeaturePlot(ABC, 
            features = 'Contamination', 
            raster=T       # 细胞过多时候需要加这个参数
) + 
  scale_color_viridis_c()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab('scVI_UMAP_1')+
  ylab('scVI_UMAP_2')

#现在我们根据Contamination进行筛选，decontX官方并没有给予明确的建议，这里参考的了soupX的官方文档，并基于自己摸索，根据这个数据的情况，选择 Contamination < 0.2是比较好的选择。
low_con_scobj = ABC[,ABC$Contamination < 0.2] #保留小于等于0.2的细胞



#marker gene输出展示pdf###########################################################################################################


j=ABC
name=i

pdf(paste(name,".pdf",sep=""),width=12,height=12)
DimPlot(j, group.by = c("seurat_clusters"),reduction = "umap",label=T,repel=T)+ NoLegend()
DimPlot(j, group.by = c( "celltype"),reduction = "umap",label=T,repel=T)+ NoLegend()
DimPlot(j, group.by = c("seurat_clusters"),reduction = "tsne",label=T,repel=T)+ NoLegend()
DimPlot(j, group.by = c("celltype"),reduction = "tsne",label=T,repel=T)+ NoLegend()
SpatialDimPlot(j, label = TRUE, label.size = 3,pt.size.factor = 2)

marker=c("MOG","OPALIN","KLK6","GPR37","ANLN")
plot(c(1:10),family="YaHei Consolas",main="少突胶质细胞 (Oligo)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("PDGFRA","GPR17","CCND1","TNR","NEU4")
plot(c(1:10),family="YaHei Consolas",main="少突胶质细胞前体 (OPC)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("SLC4A4","NTSR2","GFAP","ACSBG1","ALDOC")         #(KCNJ10仅作参考)
plot(c(1:10),family="YaHei Consolas",main="星形胶质细胞 (Astrocyte)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("C1QC","CD68","C1QB","TREM2","CCL4")
plot(c(1:10),family="YaHei Consolas",main="小胶质细胞 (Microglia)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("KLRD1","NKG7","GZMA","XCL1","GZMK")
plot(c(1:10),family="YaHei Consolas",main="自然杀伤细胞 (NK)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("CD79A","CCR7","CD27","IGLC2","CXCR4")
plot(c(1:10),family="YaHei Consolas",main="B细胞 (B cell)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("TRBC2","CD3D","CD3G","GZMK","CD3E")
plot(c(1:10),family="YaHei Consolas",main="T细胞 (T cell)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("APOBEC3A","TET2","MEFV","HCK","ITGAX","RBPJ")
plot(c(1:10),family="YaHei Consolas",main="单核细胞 (Monocytes)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("CD66b","CD11b","CD15","CD16","CXCL8","FCGR3B","MNDA")
plot(c(1:10),family="YaHei Consolas",main="中性粒细胞 (Neutrophlis)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("HCRTR2","PDGFD","NEUN","RYR3","SDK2","TRPC3","SYT10","STK32B")
plot(c(1:10),family="YaHei Consolas",main="树突状细胞 (DC)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("EPCAM","CD24")
plot(c(1:10),family="YaHei Consolas",main="上皮细胞 (Epithelial cell)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("CLDN5","FLT1","PECAM1","VWF","CLEC14A")
plot(c(1:10),family="YaHei Consolas",main="内皮细胞 (Endothelial)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("KCNMB1","MRVI1","KCNAB1","DMPK")
plot(c(1:10),family="YaHei Consolas",main="平滑肌细胞 (Smooth muscle cells)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("COL6A2","VIM","COL1A1","PDGFRB")
plot(c(1:10),family="YaHei Consolas",main="成纤维细胞 (Fibroblasts)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("TAGLN","ACTA2","PDGFRB","DES")
plot(c(1:10),family="YaHei Consolas",main="周细胞 (Mural)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("STMN2","THY1")
plot(c(1:10),family="YaHei Consolas",main="神经元细胞 (Neuron)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("SLC17A7","SATB2")
plot(c(1:10),family="YaHei Consolas",main="兴奋性神经元细胞 (Ex Neuron)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("GAD1","GAD2","SLC32A1","SLC6A1")
plot(c(1:10),family="YaHei Consolas",main="抑制性神经元细胞 (In Neuron)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("CALB1","CCK","EN1","ETV1","EVX1")
plot(c(1:10),family="YaHei Consolas",main="中间神经元 (inter Neuron)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c("GAD1","SLC32A1","PVALB","SLC6A1")
plot(c(1:10),family="YaHei Consolas",main="GABA神经元 (GABAergic Neuron)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")


dev.off()


ABC.markers <- FindAllMarkers(ABC, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ABC.markers , paste(i,"_all_markers.csv",sep=""))

pdf(paste(i,"_HeatMap.pdf",sep=""),width=25,height=35)
ABC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(ABC, features = top10$gene) + NoLegend()
dev.off()


#Seurat::DotPlot####################################
library(dplyr);library(ggtree);library(cowplot)
df <- DotPlot(ABC,features = marker)$data

mat <- df %>% dplyr::select(-pct.exp, -avg.exp.scaled) %>%  tidyr::pivot_wider(names_from = features.plot,values_from = avg.exp ) %>% data.frame() 
row.names(mat) <- mat$id
mat <- mat[,-1]

clust <- hclust(dist(mat %>% as.matrix()))
ddgram <- as.dendrogram(clust) 

ggtree_plot <- ggtree::ggtree(ddgram,branch.length="none")+geom_tippoint(color="#FDAC4F", shape=8, size=3)

df$id <-  factor(df$id , levels = clust$labels[clust$order])
dotplot <- ggplot(df,aes(x=features.plot,y = id,size = pct.exp, color = avg.exp.scaled))+
  geom_point() + 
  scale_size("% detected", range = c(0,6)) +
  scale_y_discrete(position = "right") +
  scale_color_gradientn(colours = viridis::viridis(20),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nexpression") +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("Markers") + theme_bw() +
  theme(
    axis.text.x = element_text(size=10, angle=0, hjust=0.5, color="black"),
    axis.text.y = element_text(size=12, color="black"),
    axis.title = element_text(size=14)
  )
ggtree_plot_yset <- ggtree_plot + aplot::ylim2(dotplot)
p <- plot_grid(ggtree_plot_yset,NULL,dotplot, nrow = 1, 
               rel_widths = c(0.5,-0.025, 2), 
               align = 'h')

pdf("dotplot.pdf",height=8,width=8)
print(p)
dev.off()

#HoverLocator() 鼠标悬停可见细胞id及其它 meta.data########################################################################
pbmc = ABC

p1=DimPlot(pbmc,reduction = 'tsne',raster=FALSE); p1
HoverLocator(p1)
HoverLocator(plot = p1, information = FetchData(object = pbmc, vars = 'percent.mt'))  #细胞的其他信息
HoverLocator(plot = p1, information = FetchData(object = pbmc, vars = c('manucelltype', "manusubcelltype") )) #更多信息


#高亮显示一部分细胞###########################################################


DimPlot(pbmc, label=T)
DimPlot(pbmc, raster=FALSE,cells.highlight=5) #一个细胞
DimPlot(pbmc, cells.highlight=WhichCells(pbmc, idents=5)) #一群细胞，method1
DimPlot(pbmc, cells.highlight=WhichCells(pbmc, idents=5),
        cols.highlight = "blue" ) #指定高亮细胞的颜色
DimPlot(pbmc, cells.highlight=WhichCells(pbmc, idents=5),
        cols.highlight = "blue",
        sizes.highlight = 2) #指定高亮细胞的大小

DimPlot(pbmc, cells.highlight=WhichCells(pbmc, expression = CD79A>3)) #一群细胞，method2
# 取子集后，取细胞名字
DimPlot(pbmc, #group.by = "isExpHigh", #该参数在高亮显示时不影响显示
        cells.highlight=Cells(subset(pbmc, subset=isExpHigh==T)) ) #一群细胞，method3

# 多群细胞
DimPlot(pbmc,
        cells.highlight=list(
          c5=WhichCells(pbmc, idents=5),
          c8=WhichCells(pbmc, idents=8),
          c2=WhichCells(pbmc, idents=2)
        ),
        sizes.highlight = c(4.5, 1.5, 0.5), #这就很奇怪了，尺寸和出现的对应
        #cols.highlight=list(c5="red", c8="orange")
        cols.highlight=c("orange", "blue", "purple") #颜色对应关系总是不对
        # 貌似它总是按照order(listNames)的顺序染色
)


#批处理marker gene输出展示pdf########################################################################################
for (i in 1:length(brain)) {
  brain[[i]][["percent.mt"]] <- PercentageFeatureSet(brain[[i]], pattern = "^MT-")
  brain[[i]] <- subset(brain[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
  brain[[i]] <- SCTransform(brain[[i]], assay = "RNA", return.only.var.genes = FALSE, verbose = FALSE,)
  brain[[i]] <- RunPCA(brain[[i]], assay = "RNA", verbose = FALSE)
  brain[[i]] <- FindNeighbors(brain[[i]], reduction = "pca", dims = 1:20)
  brain[[i]] <- FindClusters(brain[[i]], verbose = FALSE, resolution = 20)
  brain[[i]] <- RunUMAP(brain[[i]], reduction = "pca", dims = 1:20)
  brain[[i]] <- RunTSNE(object=brain[[i]], reduction = "pca" , dims=1:20)
  
  pdf(paste(brain[[i]]@meta.data$orig.ident[1],".pdf",sep=""),width=12,height=12)
  DimPlot(brain[[i]], group.by = c("seurat_clusters"),reduction = "umap",label=T,repel=T)+ NoLegend()
  DimPlot(brain[[i]], group.by = c("seurat_clusters"),reduction = "tsne",label=T,repel=T)+ NoLegend()

  marker=c("MOG","OPALIN","KLK6","GPR37","ANLN")
  plot(c(1:10),family="YaHei Consolas",main="少突胶质细胞 (Oligo)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("PDGFRA","GPR17","CCND1","TNR","NEU4")
  plot(c(1:10),family="YaHei Consolas",main="少突胶质细胞前体 (OPC)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("SLC4A4","NTSR2","GFAP","ACSBG1","ALDOC")         #(KCNJ10仅作参考)
  plot(c(1:10),family="YaHei Consolas",main="星形胶质细胞 (Astrocyte)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("C1QC","CD68","C1QB","TREM2","CCL4")
  plot(c(1:10),family="YaHei Consolas",main="小胶质细胞 (Microglia)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("KLRD1","NKG7","GZMA","XCL1","GZMK")
  plot(c(1:10),family="YaHei Consolas",main="自然杀伤细胞 (NK)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("CD79A","CCR7","CD27","IGLC2","CXCR4")
  plot(c(1:10),family="YaHei Consolas",main="B细胞 (B cell)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("TRBC2","CD3D","CD3G","GZMK","CD3E")
  plot(c(1:10),family="YaHei Consolas",main="T细胞 (T cell)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("APOBEC3A","TET2","MEFV","HCK","ITGAX","RBPJ")
  plot(c(1:10),family="YaHei Consolas",main="单核细胞 (Monocytes)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("CD66b","CD11b","CD15","CD16","CXCL8","FCGR3B","MNDA")
  plot(c(1:10),family="YaHei Consolas",main="中性粒细胞 (Neutrophlis)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("HCRTR2","PDGFD","NEUN","RYR3","SDK2","TRPC3","SYT10","STK32B")
  plot(c(1:10),family="YaHei Consolas",main="树突状细胞 (DC)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("EPCAM","CD24")
  plot(c(1:10),family="YaHei Consolas",main="上皮细胞 (Epithelial cell)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("CLDN5","FLT1","PECAM1","VWF","CLEC14A")
  plot(c(1:10),family="YaHei Consolas",main="内皮细胞 (Endothelial)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("KCNMB1","MRVI1","KCNAB1","DMPK")
  plot(c(1:10),family="YaHei Consolas",main="平滑肌细胞 (Smooth muscle cells)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("COL6A2","VIM","COL1A1","PDGFRB")
  plot(c(1:10),family="YaHei Consolas",main="成纤维细胞 (Fibroblasts)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("TAGLN","ACTA2","PDGFRB","DES")
  plot(c(1:10),family="YaHei Consolas",main="周细胞 (Mural)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("STMN2","THY1")
  plot(c(1:10),family="YaHei Consolas",main="神经元细胞 (Neuron)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("SLC17A7","SATB2")
  plot(c(1:10),family="YaHei Consolas",main="兴奋性神经元细胞 (Ex Neuron)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("GAD1","GAD2","SLC32A1","SLC6A1")
  plot(c(1:10),family="YaHei Consolas",main="抑制性神经元细胞 (In Neuron)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("CALB1","CCK","EN1","ETV1","EVX1")
  plot(c(1:10),family="YaHei Consolas",main="中间神经元 (inter Neuron)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  marker=c("GAD1","SLC32A1","PVALB","SLC6A1")
  plot(c(1:10),family="YaHei Consolas",main="GABA神经元 (GABAergic Neuron)")
  VlnPlot(brain[[i]], features =marker,ncol=1)
  FeaturePlot(brain[[i]], features =marker)
  FeaturePlot(brain[[i]], features =marker,min.cutoff="q10",max.cutoff="q95")
  
  
  dev.off()
  
  
  ABC.markers <- FindAllMarkers(brain[[i]], only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  ABC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
  write.csv(ABC.markers , paste(brain[[i]]@meta.data$orig.ident[1],"_all_markers.csv",sep=""))
  write.csv(top10 , paste(brain[[i]]@meta.data$orig.ident[1],"_top10_markers.csv",sep=""))
  

}
# 啊#######################################################

i=1
ggsave(paste(i,"_.png"),DimPlot(brain[[i]], reduction = "umap", group.by = c("seurat_clusters", "orig.ident",'manucelltype'),label = T,repel = T)/DimPlot(brain[[i]], reduction = "tsne", group.by = c("seurat_clusters", "orig.ident"),label = T,repel = T),width=14,height=10)
j=brain[[i]]
name=i

pdf(paste(name,".pdf",sep=""),width=12,height=12)

marker=c('PRPH','TUBB3','RBFOX3','NEFM','DISP2','NRGN','ELAVL3')
plot(c(1:10),family="YaHei Consolas",main="神经元 Neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('ROBO1','NEUROD1','NTRK1','SLC1A1','NOTCH3','NOTCH1','LHX6','NES')
plot(c(1:10),family="YaHei Consolas",main="未成熟神经元 Immature neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('DDC','PNMT','DBH','SLC18A','NPFF','SLC12A7','SYT1','TH')
plot(c(1:10),family="YaHei Consolas",main="肾上腺素能神经元 Adrenergic neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('CHAT','SLC18A3','SLC5A7','TAC1','ACLY','ACHE','BRCA1')
plot(c(1:10),family="YaHei Consolas",main="胆碱能神经元 Cholinergic neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('TH','SLC6A3','PITX3','SMAD3','NEUROD6','SLC18A3','SLC18A2','DDC')
plot(c(1:10),family="YaHei Consolas",main="多巴胺能神经元Dopaminergic neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('GAD1','GAD2','PAX2','SLC32A1','SLC6A1','VIP','PVALB')
plot(c(1:10),family="YaHei Consolas",main="GABA能神经元 GABAergic neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('SLC17A7','SLC17A6','MEIS2','SLC17A8','SLC1A1','SLC1A2','SLC1A6')
plot(c(1:10),family="YaHei Consolas",main="谷氨酰胺能神经元 Glutaminergic neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('SLC6A9','SLC32A1')
plot(c(1:10),family="YaHei Consolas",main="甘氨酸能神经元 Glycinergic neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('PVALB','SST','VIP','EOMES','KCNC2','CALB1','CCK','VSX2')
plot(c(1:10),family="YaHei Consolas",main="中间神经元 Interneurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('CHAT','MNX1','ISL2','VSX2','EN1','EVX1','EVX2','FGF1')
plot(c(1:10),family="YaHei Consolas",main="运动神经元 Motor neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('TH','DDC','DBH','SLC18A2','SLC18A3','SLC6A2','SLC39A11','SLC9B2')
plot(c(1:10),family="YaHei Consolas",main="去甲肾上腺素能神经元 Noradrenergic neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('TSHZ1','CALB1','ADGRB1','HCN1','SLC24A2','GABRA1','CLSTN3')
plot(c(1:10),family="YaHei Consolas",main="浦肯野神经元 Purkinje neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('TPH1','SLC6A4','FEV','TPH2','DDC','SLC18A2','ESM1','SLC22A3')
plot(c(1:10),family="YaHei Consolas",main="5-羟色胺能神经元 Serotonergic neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('FGF13','CREBRF','KCNAB2','ATL1','HSPB8','FKBP1B','TLN2','SYNM')
plot(c(1:10),family="YaHei Consolas",main="三叉神经元 Trigeminal neurons")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('GFAP','ACSBG1','SLC1A2','GJA1','SLC1A3','S100B','SOX9','S1PR1')
plot(c(1:10),family="YaHei Consolas",main="星形胶质细胞Astrocytes")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('ITGAM','CMTM6','AIF1','ADRB2','CSF1R','P2RY12','CD40','TMEM119')
plot(c(1:10),family="YaHei Consolas",main="小胶质细胞 Microglia")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('MBP','MOG','MAG','NIPAL4','CLDN14','KLK6','PEX5L','SEC14L5')
plot(c(1:10),family="YaHei Consolas",main="少突胶质细胞 Oligodendrocytes")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('PDGFRA','NEU4','NKX6-2','FYN','TNR','ALDOC','PCDH15','OLIG1')
plot(c(1:10),family="YaHei Consolas",main="少突胶质祖细胞 Oligodendrocyte progenitor cells")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('MPZ','CRYAB','NF2','ALDOC','APOD','SOX10','CMTM5','BCHE','EGFL8')
plot(c(1:10),family="YaHei Consolas",main="雪旺细胞 Schwann cells")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")

marker=c('RABL2','FOXJ1','PIFO','RSPH1','TEKT4','CALML4','CRYGN','HDC','FAM183B')
plot(c(1:10),family="YaHei Consolas",main="室管膜细胞 Ependymal cells")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")


dev.off()


























