#0. 载入
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)


##############################################################################
#1. 环境
setwd("G:/danxibao/data2/")
FCD_16=readRDS("./FCD_16.rds")
normal_16=readRDS("./normal_16.rds")
merge = readRDS("./merge.rds")

ABC=FCD_16
DimPlot(ABC, group.by = c("seurat_clusters", "celltype","orig.ident","manucelltype"),reduction = "umap",label=T,repel = T)

##########################################################################
#2. 数据录入
#pbmc <- readRDS("pbmc.rds")
pbmc = ABC
color_by = "manusubcelltype"

data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)



###########################################################################
#####3. 预处理

#RNA-seq是使用PCA，如果是处理ATAC-seq的数据用Latent Semantic Indexing

#⚠️preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
plot_cells(cds)

colnames(colData(cds))

#以之前的Seurat分群来添加颜色，和原有的Seurat分群对比
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by=color_by) + ggtitle('cds.umap')

##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by=color_by) + ggtitle('int.umap')
p = p1|p2
p

#可视化指定基因
ciliated_genes <- c("SLC4A4","GFAP","PLP1","PDGFRA","PTPRC","CLDN5")
plot_cells(cds,genes=ciliated_genes,label_cell_groups=FALSE,show_trajectory_graph=FALSE)


#################################################################################
#######4. 构建细胞轨迹
cds <- cluster_cells(cds)
#cds <- cluster_cells(cds,cluster_method = "louvain")
cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = color_by, label_groups_by_cluster=FALSE,label_leaves=FALSE, label_branch_points=FALSE,group_label_size = 3.5,graph_label_size = 4)
p

plot_cells(cds, color_cells_by = color_by, label_groups_by_cluster=FALSE,label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=3.5)
#黑色的线显示的是graph的结构。数字带白色圆圈表示不同的结局，也就是叶子。
#数字带黑色圆圈代表分叉点，从这个点开始，细胞可以有多个结局。
#这些数字可以通过label_leaves和label_branch_points参数设置。

######4.1使用 3D 轨迹
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")
cds_3d_plot_obj

########4.2 细胞按拟时排序
cds <- order_cells(cds)
#运行上面的代码后，会跳出这个窗口。可以手动在图上选择一个位置，然后点击Done。可以选择多个位置。

p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
p

################################################################################
#5. 差异表达分析

Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
#write.csv(Track_genes, "Trajectory_genes.csv", row.names = F)

Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%pull(gene_short_name) %>% as.character()


#基因表达趋势图
p <- plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="seurat_clusters", min_expr=0.5, ncol = 2)
p


#FeaturePlot图
p <- plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)
p$facet$params$ncol <- 5
p


#寻找共表达基因模块
#Track_genes <- read.csv("Trajectory_genes.csv")
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-1, cores = 6)
#write.csv(gene_module, "Genes_Module.csv", row.names = F)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$seurat_clusters)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
p

#提取拟时分析结果返回seurat对象

pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(pbmc@meta.data)]
pbmc$pseudotime <- pseudotime

pbmc


###########################################
time_b=cds@principal_graph_aux@listData$UMAP$pseudotime
time_b = as.data.frame(time_b)

c=c('TAFA3','LDLR','CST7','CX3CL1','GRN','PTPRC','NR1D1','SYT11')
negative = apply(cds@assays@data@listData$counts[c,],2,mean)
negative = as.data.frame(negative)
time_b[,2] = negative

c=c('CCL3','TAFA3','TTBK1','LRRK2','CTSC')
positive = apply(cds@assays@data@listData$counts[c,],2,mean)
positive = as.data.frame(positive)
time_b[,3] = positive

c=c('CX3CR1','IL33','GRN','TYROBP','TREM2')
immu = apply(cds@assays@data@listData$counts[c,],2,mean)
immu = as.data.frame(immu)
time_b[,4] = immu

c=c('CX3CR1','CCL3','P2RX4','CX3CL1','P2RY12')
positive_migration = apply(cds@assays@data@listData$counts[c,],2,mean)
positive_migration = as.data.frame(positive_migration)
time_b[,5] = positive_migration

pdf('negative regulation of microglial cell activation.pdf',width = 15,height = 8)
ggplot(time_b , aes(x=time_b, y=negative)) +
geom_point(size=0.1) +
geom_smooth(aes(x=time_b, y=negative), se = F, method = 'loess') +ggtitle('negative regulation of microglial cell activation')
dev.off()
#########debug

trace('calculateLW', edit = T, where = asNamespace("monocle3"))