library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)


path <- r"(G:\filterd data\67\output-XETG00208__0022108__TLE_pre__20240126__053708)"
# Load the Xenium data
xenium.obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)


#接下来，我们使用 绘制组织上泛抑制性神经元标记 Gad1、抑制性神经元亚型标记 Pvalb 和 Sst 以及星形胶质细胞标记 Gfap 的位置ImageDimPlot()。
ImageDimPlot(xenium.obj, fov = "fov", molecules = c("Gad1", "Sst", "Pvalb", "Gfap"), nmols = 20000)
#ImageFeaturePlot()在这里，我们使用类似于可视化 2D 嵌入表达的函数，在每个细胞水平上可视化一些关键层标记基因的表达水平FeaturePlot()。
#我们手动将每个基因的 大约调整到其计数分布的max.cutoff第 90 个百分位数（可以用 指定）以提高对比度。max.cutoff='q90'
ImageFeaturePlot(xenium.obj, features = c("Cux2", "Rorb", "Bcl11b", "Foxp2"), max.cutoff = c(25,35, 12, 10), size = 0.75, cols = c("white", "red"))

#我们可以使用该功能放大选定的区域Crop()。放大后，我们可以可视化细胞分割边界以及单个分子。
cropped.coords <- Crop(xenium.obj[["fov"]], x = c(5000, 7500), y = c(1500, 3000), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",coord.fixed = FALSE, molecules = c("Gad1", "Sst", "Npy2r", "Pvalb", "Nrn1"), nmols = 10000)
ImageDimPlot(xenium.obj,fov='zoom', cols = "polychrome", size = 0.75, axes = TRUE ,border.color = "white",border.size = 0.1)
#接下来，我们使用 SCTransform 进行归一化，然后进行标准降维和聚类。此步骤从开始到结束大约需要 5 分钟。
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)

#DimPlot()然后，我们可以根据每个细胞在 UMAP 空间中的聚类对每个细胞进行着色，或者用 覆盖在图像上，从而可视化聚类的结果ImageDimPlot()。
DimPlot(xenium.obj)

FeaturePlot(xenium.obj, features = c("Cux2", "Bcl11b", "Foxp2", "Gad1", "Sst", "Gfap"))

#现在，我们可以使用ImageDimPlot()对上一步中确定的簇标签着色的单元格位置进行着色。
DefaultBoundary(xenium.obj[["fov"]]) <- "centroids"
DefaultBoundary(xenium.obj[["fov"]]) <- "segmentation"
ImageDimPlot(xenium.obj, cols = "polychrome", size = 0.75, axes = TRUE ,border.color = "white",border.size = 0.1)

##########从 Seurat 查询和参考对象中提取计数、聚类和点信息，以构建ReferenceRCTDSpatialRNA用于注释的对象。然后将注释的输出添加到 Seurat 对象。

library(spacexr)

query.counts <- GetAssayData(xenium.obj, assay = "Xenium", slot = "counts")[, Cells(xenium.obj[["fov"]])]
coords <- GetTissueCoordinates(xenium.obj[["fov"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# allen.corted.ref can be downloaded here:
# https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1
allen.cortex.ref <- readRDS(r"(G:\work\UNSC reference\mouse Allen\allen_cortex.rds)")
allen.cortex.ref <- UpdateSeuratObject(allen.cortex.ref)

Idents(allen.cortex.ref) <- "subclass"
# remove CR cells because there aren't enough of them for annotation
allen.cortex.ref <- subset(allen.cortex.ref, subset = subclass != "CR")
counts <- GetAssayData(allen.cortex.ref, assay = "RNA", slot = "counts")
cluster <- as.factor(allen.cortex.ref$subclass)
names(cluster) <- colnames(allen.cortex.ref)
nUMI <- allen.cortex.ref$nCount_RNA
names(nUMI) <- colnames(allen.cortex.ref)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
xenium.obj$predicted.celltype <- annotations
keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
xenium.obj <- subset(xenium.obj, cells = keep.cells)

#
xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "predicted.celltype",niches.k = 5, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", size = 1.5, cols = "polychrome",
                              dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
celltype.plot | niche.plot

table(xenium.obj$predicted.celltype, xenium.obj$niches)

##########

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)

new.cluster.ids <- c('Oligodendrocyte',
                     'Microglia',
                     'EC',
                     'Astrocyte',
                     'Fibroblasts',
                     'Neuron',
                     'Neuron',
                     'Astrocyte',
                     'Neuron',#lamp5,cux2
                     'Neuron',
                     'OPC',
                     'Macrophage|Neuron',
                     'In Neuron',
                     'Neuron',
                     'In Neuron',
                     'Neuron',
                     'Neuron',
                     'Oligodendrocyte',
                     'Fibro|Neuron',
                     'EC',
                     'Neuron',
                     'EC|Neuron|Oligo',
                     'Fibro|EC|Neuron',
                     'Neutrophils|B cells',
                     'In Neuron',
                     'EC|Neuron')
ABC@meta.data$manucelltype <- plyr::mapvalues(x = ABC@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)


setwd(r"(G:\work\24.2\Xenium)")
marker = read.csv('Xenium_mBrain_v1.1_metadata.csv')
merge = readRDS('merge4.rds')
genelist = split(marker,marker$Annotation)

#################################
xenium.obj = merge 
genelist = split(marker,marker$Annotation)
dir.create('Xenium_gene')
setwd('./Xenium_gene')


for(j in 1:length(genelist)) {
  tamptlist = genelist[[j]]
  tamptgenelist = tamptlist$Genes
  dirnames = names(genelist[j])
  dirnames=gsub("/","_",dirnames)
  dir.create(dirnames)
  setwd(paste('./',dirnames,sep = ''))
  for(i in 1:length(tamptgenelist)){
    marker = tamptgenelist[i]
    DefaultBoundary(xenium.obj[["fov"]]) <- "centroids"
    DefaultBoundary(xenium.obj[["fov.2"]]) <- "centroids"
    DefaultBoundary(xenium.obj[["fov.3"]]) <- "centroids"
    DefaultBoundary(xenium.obj[["fov.4"]]) <- "centroids"
    
    TLE_post_centroids = ImageDimPlot(xenium.obj, fov = "fov", molecules =marker , nmols = 20000, cols = "polychrome")
    TLE_pre_centroids = ImageDimPlot(xenium.obj, fov = "fov.2", molecules =marker , nmols = 20000, cols = "polychrome")
    Ctri_post_centroids = ImageDimPlot(xenium.obj, fov = "fov.3", molecules = marker, nmols = 20000, cols = "polychrome")
    Ctri_pre_centroids = ImageDimPlot(xenium.obj, fov = "fov.4", molecules = marker, nmols = 20000, cols = "polychrome")
    
    DefaultBoundary(xenium.obj[["fov"]]) <- "segmentation"
    DefaultBoundary(xenium.obj[["fov.2"]]) <- "segmentation"
    DefaultBoundary(xenium.obj[["fov.3"]]) <- "segmentation"
    DefaultBoundary(xenium.obj[["fov.4"]]) <- "segmentation"
    
    TLE_post_segmentation = ImageDimPlot(xenium.obj, fov = "fov", molecules = marker, nmols = 20000, cols = "polychrome")
    TLE_pre_segmentation = ImageDimPlot(xenium.obj, fov = "fov.2", molecules = marker, nmols = 20000, cols = "polychrome")
    Ctri_post_segmentation = ImageDimPlot(xenium.obj, fov = "fov.3", molecules = marker, nmols = 20000, cols = "polychrome")
    Ctri_pre_segmentation = ImageDimPlot(xenium.obj, fov = "fov.4", molecules = marker, nmols = 20000, cols = "polychrome")
    
    ggsave(paste(marker,'TLE_post_centroids.png',sep='_'),TLE_post_centroids ,width=25,height=16.5)
    ggsave(paste(marker,'TLE_pre_centroids.png',sep='_'), TLE_pre_centroids,width=25,height=16.5)
    ggsave(paste(marker,'Ctri_post_centroids.png',sep='_'),Ctri_post_centroids ,width=25,height=16.5)
    ggsave(paste(marker,'Ctri_pre_centroids.png',sep='_'), Ctri_pre_centroids,width=25,height=16.5)
    ggsave(paste(marker,'TLE_post_segmentation.png',sep='_'), TLE_post_segmentation,width=25,height=16.5)
    ggsave(paste(marker,'TLE_pre_segmentation.png',sep='_'),TLE_pre_segmentation ,width=25,height=16.5)
    ggsave(paste(marker,'Ctri_post_segmentation.png',sep='_'),Ctri_post_segmentation ,width=25,height=16.5)
    ggsave(paste(marker,'Ctri_pre_segmentation.png',sep='_'),Ctri_pre_segmentation ,width=25,height=16.5)
    
    
  }
  
  setwd('../')
}

setwd('../')










for(i in 1:length(markers)){
  marker = markers[i]
  DefaultBoundary(xenium.obj[["E_R"]]) <- "centroids"
  DefaultBoundary(xenium.obj[["fov.2"]]) <- "centroids"
  DefaultBoundary(xenium.obj[["C_L"]]) <- "centroids"
  DefaultBoundary(xenium.obj[["fov.4"]]) <- "centroids"
  
  TLE_post_centroids = ImageDimPlot(xenium.obj, fov = "E_R", molecules =marker , nmols = 20000, cols = "polychrome",group.by = 'seurat_clusters')
  TLE_pre_centroids = ImageDimPlot(xenium.obj, fov = "fov.2", molecules =marker , nmols = 20000, cols = "polychrome",group.by = 'seurat_clusters')
  Ctri_post_centroids = ImageDimPlot(xenium.obj, fov = "C_L", molecules =marker , nmols = 20000, cols = "polychrome",group.by = 'seurat_clusters')
  Ctri_pre_centroids = ImageDimPlot(xenium.obj, fov = "fov.4", molecules =marker , nmols = 20000, cols = "polychrome",group.by = 'seurat_clusters')
  
  DefaultBoundary(xenium.obj[["E_R"]]) <- "segmentation"
  DefaultBoundary(xenium.obj[["fov.2"]]) <- "segmentation"
  DefaultBoundary(xenium.obj[["C_L"]]) <- "segmentation"
  DefaultBoundary(xenium.obj[["fov.4"]]) <- "segmentation"
  
  TLE_post_segmentation = ImageFeaturePlot(xenium.obj, fov = "E_R", features = marker, cols = c("white", "red"), nmols = 20000)
  TLE_pre_segmentation = ImageFeaturePlot(xenium.obj, fov = "fov.2", features = marker, cols = c("white", "red"), nmols = 20000)
  Ctri_post_segmentation = ImageFeaturePlot(xenium.obj, fov = "C_L", features = marker, cols = c("white", "red"), nmols = 20000)
  Ctri_pre_segmentation = ImageFeaturePlot(xenium.obj, fov = "fov.4", features = marker, cols = c("white", "red"), nmols = 20000)
  
  ggsave(paste(marker,'E_R_centroids.png',sep='_'),TLE_post_centroids ,width=15,height=9.9)
  ggsave(paste(marker,'E_A_centroids.png',sep='_'), TLE_pre_centroids,width=25,height=16.5)
  ggsave(paste(marker,'C_L_centroids.png',sep='_'),Ctri_post_centroids ,width=15,height=9.9)
  ggsave(paste(marker,'C_A_centroids.png',sep='_'), Ctri_pre_centroids,width=25,height=16.5)
  ggsave(paste(marker,'E_R_segmentation.png',sep='_'), TLE_post_segmentation,width=25,height=16.5)
  ggsave(paste(marker,'E_A_segmentation.png',sep='_'),TLE_pre_segmentation ,width=25,height=16.5)
  ggsave(paste(marker,'C_L_segmentation.png',sep='_'),Ctri_post_segmentation ,width=25,height=16.5)
  ggsave(paste(marker,'C_A_segmentation.png',sep='_'),Ctri_pre_segmentation ,width=25,height=16.5)
  
  ggsave(paste(marker,'E_R_centroids.pdf',sep='_'),TLE_post_centroids ,width=15,height=9.9)
  ggsave(paste(marker,'E_A_centroids.pdf',sep='_'), TLE_pre_centroids,width=25,height=16.5)
  ggsave(paste(marker,'C_L_centroids.pdf',sep='_'),Ctri_post_centroids ,width=15,height=9.9)
  ggsave(paste(marker,'C_A_centroids.pdf',sep='_'), Ctri_pre_centroids,width=25,height=16.5)
  ggsave(paste(marker,'E_R_segmentation.pdf',sep='_'), TLE_post_segmentation,width=25,height=16.5)
  ggsave(paste(marker,'E_A_segmentation.pdf',sep='_'),TLE_pre_segmentation ,width=25,height=16.5)
  ggsave(paste(marker,'C_L_segmentation.pdf',sep='_'),Ctri_post_segmentation ,width=25,height=16.5)
  ggsave(paste(marker,'C_A_segmentation.pdf',sep='_'),Ctri_pre_segmentation ,width=25,height=16.5)
  
}

setwd('../')







