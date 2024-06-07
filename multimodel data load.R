library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)

setwd("G:/work/ATAC")


#############我们将根据基因表达数据创建一个 Seurat 对象，然后添加 ATAC-seq 数据作为第二次检测#####

# the 10x hdf5 file contains both data types.
path = R"(G:\filterd data\5\FCD)"
inputdata.10x <- Read10X(data.dir = paste(path,"/filtered_feature_bc_matrix",sep=""))
#or
#inputdata.10x <- Read10X_h5(paste(path,"/filtered_feature_bc_matrix.h5",sep=""))
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file = paste(path,"/atac_fragments.tsv.gz",sep="")
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay



#####我们根据每种模式检测到的分子数量以及线粒体百分比执行基本 QC
VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,log = TRUE, pt.size = 0) + NoLegend()
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)


#####接下来，我们使用 RNA 和 ATAC-seq 数据的标准方法独立地对两种分析进行预处理和降维
# RNA analysis
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")



#####我们计算了一个 WNN 图，表示 RNA 和 ATAC-seq 模态的加权组合。我们将此图用于 UMAP 可视化和聚类
pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)


p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


#####我们注释下面的集群
# perform sub-clustering on cluster 6 to find additional structure
pbmc <- FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(pbmc) <- "sub.cluster"

# add annotations
pbmc <- RenameIdents(pbmc, '19' = 'pDC','20' = 'HSPC','15' = 'cDC')
pbmc <- RenameIdents(pbmc, '0' = 'CD14 Mono', '9' ='CD14 Mono', '5' = 'CD16 Mono')
pbmc <- RenameIdents(pbmc, '10' = 'Naive B', '11' = 'Intermediate B', '17' = 'Memory B', '21' = 'Plasma')
pbmc <- RenameIdents(pbmc, '7' = 'NK')
pbmc <- RenameIdents(pbmc, '4' = 'CD4 TCM', '13'= "CD4 TEM", '3' = "CD4 TCM", '16' ="Treg", '1' ="CD4 Naive", '14' = "CD4 Naive")
pbmc <- RenameIdents(pbmc, '2' = 'CD8 Naive', '8'= "CD8 Naive", '12' = 'CD8 TEM_1', '6_0' = 'CD8 TEM_2', '6_1' ='CD8 TEM_2', '6_4' ='CD8 TEM_2')
pbmc <- RenameIdents(pbmc, '18' = 'MAIT')
pbmc <- RenameIdents(pbmc, '6_2' ='gdT', '6_3' = 'gdT')
pbmc$celltype <- Idents(pbmc)



#####我们可以可视化基于基因表达、ATAC-seq 或 WNN 分析的聚类。
#####差异比之前的分析更细微（您可以探索权重，它们比我们的 CITE-seq 示例中的分配更均匀）
#####但我们发现 WNN 提供了最清晰的细胞状态分离
p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


#####例如，ATAC-seq 数据有助于分离 CD4 和 CD8 T 细胞状态。
#####这是由于存在多个基因座，这些基因座在不同的 T 细胞亚型之间表现出不同的可及性。
#####我们可以使用Signac 可视化小插图中的工具，将 CD8A 基因座的“假体”轨迹与基因表达水平的小提琴图一起可视化
## to make the visualization easier, subset T cell clusters
celltype.names <- levels(pbmc)
tcell.names <- grep("CD4|CD8|Treg", celltype.names,value = TRUE)
tcells <- subset(pbmc, idents = tcell.names)
CoveragePlot(tcells, region = 'CD8A', features = 'CD8A', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)





setwd("C:/Users/11723/Desktop/ATAC")
FCD=readRDS("FCD.rds")
noFCD=readRDS("noFCD.rds")
tsc=readRDS("tsc.rds")
notsc=readRDS("notsc.rds")


DimPlot(ABC, reduction = "wnn.umap", group.by = c("seurat_clusters", "manusubcelltype2","orig.ident","manucelltype"), label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN")

gene="MTOR"
CoveragePlot(FCD,region = gene, features = gene, assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)/CoveragePlot(noFCD,region = gene, features = gene, assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)




pbmc = tsc

Idents(pbmc) = pbmc@meta.data$manucelltype
DefaultAssay(pbmc) <- "ATAC"
# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)
# link peaks to genes
pbmc <- LinkPeaks(
object = pbmc,
peak.assay = "ATAC",
expression.assay = "SCT",
genes.use = c("CCL4","MTOR")
)

p1 <- CoveragePlot(
object = pbmc,
region = "CCL4",
features = "CCL4",
expression.assay = "SCT",
extend.upstream = 500,
extend.downstream = 10000
)
p2 <- CoveragePlot(
object = pbmc,
region = "MTOR",
features = "MTOR",
expression.assay = "SCT",
extend.upstream = 500,
extend.downstream = 10000
)
patchwork::wrap_plots(p1, p2, ncol = 1)




# extract position frequency matrices for the motifs
pbmc = FCD
DefaultAssay(pbmc)="ATAC"
genelist = c("CCL4","MTOR")
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
pbmc <- AddMotifs(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

# gather the footprinting information for sets of motifs
pbmc <- Footprint(
  object = pbmc,
  motif.name = genelist,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(pbmc, features = genelist)

p2 + patchwork::plot_layout(ncol = 1)
























###################
setwd('G:/')
pdf(paste('notsc',".pdf",sep=""),width=12,height=12)
DefaultAssay(pbmc) = 'SCT'
DimPlot(pbmc, group.by = c("seurat_clusters"),reduction = "wnn.umap",label=T,repel=T)+ NoLegend()
DimPlot(pbmc, group.by = c("seurat_clusters"),reduction = "umap.rna",label=T,repel=T)+ NoLegend()
DimPlot(pbmc, group.by = c("seurat_clusters"),reduction = "umap.atac",label=T,repel=T)+ NoLegend()

marker=c("MOG","OPALIN","KLK6","GPR37","ANLN")
plot(c(1:10),family="YaHei Consolas",main="少突胶质细胞 (Oligo)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("PDGFRA","GPR17","CCND1","TNR","NEU4")
plot(c(1:10),family="YaHei Consolas",main="少突胶质细胞前体 (OPC)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("SLC4A4","NTSR2","GFAP","ACSBG1","ALDOC")         #(KCNJ10仅作参考)
plot(c(1:10),family="YaHei Consolas",main="星形胶质细胞 (Astrocyte)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("C1QC","CD68","C1QB","TREM2","CCL4")
plot(c(1:10),family="YaHei Consolas",main="小胶质细胞 (Microglia)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("KLRD1","NKG7","GZMA","XCL1","GZMK")
plot(c(1:10),family="YaHei Consolas",main="自然杀伤细胞 (NK)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("CD79A","CCR7","CD27","IGLC2","CXCR4")
plot(c(1:10),family="YaHei Consolas",main="B细胞 (B cell)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("TRBC2","CD3D","CD3G","GZMK","CD3E")
plot(c(1:10),family="YaHei Consolas",main="T细胞 (T cell)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("APOBEC3A","TET2","MEFV","HCK","ITGAX","RBPJ")
plot(c(1:10),family="YaHei Consolas",main="单核细胞 (Monocytes)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("CD66b","CD11b","CD15","CD16","CXCL8","FCGR3B","MNDA")
plot(c(1:10),family="YaHei Consolas",main="中性粒细胞 (Neutrophlis)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("HCRTR2","PDGFD","NEUN","RYR3","SDK2","TRPC3","SYT10","STK32B")
plot(c(1:10),family="YaHei Consolas",main="树突状细胞 (DC)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("EPCAM","CD24")
plot(c(1:10),family="YaHei Consolas",main="上皮细胞 (Epithelial cell)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("CLDN5","FLT1","PECAM1","VWF","CLEC14A")
plot(c(1:10),family="YaHei Consolas",main="内皮细胞 (Endothelial)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("KCNMB1","MRVI1","KCNAB1","DMPK")
plot(c(1:10),family="YaHei Consolas",main="平滑肌细胞 (Smooth muscle cells)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("COL6A2","VIM","COL1A1","PDGFRB")
plot(c(1:10),family="YaHei Consolas",main="成纤维细胞 (Fibroblasts)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("TAGLN","ACTA2","PDGFRB","DES")
plot(c(1:10),family="YaHei Consolas",main="周细胞 (Mural)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("STMN2","THY1")
plot(c(1:10),family="YaHei Consolas",main="神经元细胞 (Neuron)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("SLC17A7","SATB2")
plot(c(1:10),family="YaHei Consolas",main="兴奋性神经元细胞 (Ex Neuron)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("GAD1","GAD2","SLC32A1","SLC6A1")
plot(c(1:10),family="YaHei Consolas",main="抑制性神经元细胞 (In Neuron)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("CALB1","CCK","EN1","ETV1","EVX1")
plot(c(1:10),family="YaHei Consolas",main="中间神经元 (inter Neuron)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")

marker=c("GAD1","SLC32A1","PVALB","SLC6A1")
plot(c(1:10),family="YaHei Consolas",main="GABA神经元 (GABAergic Neuron)")
VlnPlot(pbmc, features =marker,ncol=1)
FeaturePlot(pbmc, features =marker,reduction = "wnn.umap")
FeaturePlot(pbmc, features =marker,min.cutoff="q10",max.cutoff="q95",reduction = "wnn.umap")


dev.off()
#new load##########################
path = R"(G:\filterd data\5\noFCD)"


counts <- Read10X(data.dir = paste(path,"/filtered_feature_bc_matrix",sep=""))
fragpath = paste(path,"/atac_fragments.tsv.gz",sep="")



# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

pbmc



DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
#pbmc <- TSSEnrichment(pbmc)

DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)


VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)


# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc


DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:50)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.2)
DefaultAssay(pbmc) <- "ATAC"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)

# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()



####
DefaultAssay(pbmc) <- "ATAC"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

#
idents.plot = pbmc@meta.data$seurat_clusters

marker = 'CX3CL1'

# link peaks to genes
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = marker
)

CoveragePlot(
  object = pbmc,
  region = marker,
  features = marker,
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)


CoverageBrowser(pbmc,marker)

