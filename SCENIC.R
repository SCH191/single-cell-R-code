library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

setwd("G:/SCENIC")
#########################           导入Seurat文件
RDSpathway = "G:/danxibao/data2/normal_16.rds"
pbmc <-readRDS(RDSpathway) 
#or
#pbmc = FCD_16
memory.limit(102400)
########################         准备细胞meta信息

cellInfo <- data.frame(pbmc@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="nFeature_RNA")] <- "nGene"
colnames(cellInfo)[which(colnames(cellInfo)=="nCount_RNA")] <- "nUMI"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="manucelltype")] <- "manucelltype"
cellInfo <- cellInfo[,c("sample","nGene","nUMI","cluster","manucelltype")]
saveRDS(cellInfo, file="./int/cellInfo.Rds")
#设置celltype颜色，也可不设置，使用默认颜色
#colVars <- list(celltype=c("Naive CD4 T"="forestgreen",  "CD14+ Mono"="darkorange",  "Memory CD4 T"="magenta4", "B"="hotpink",  "CD8 T"="red3",  "FCGR3A+ Mono"="skyblue",   "NK"="darkblue",  "DC"="yellow",  "Platelet"="grey"))
#colVars$celltype <- colVars$celltype[intersect(names(colVars$celltype), cellInfo$celltype)]
#saveRDS(colVars, file="int/colVars.Rds")

##################################################    初始化 SCENIC 设置,设置分析环境

org <- "hgnc" #mgi 代表鼠标，hgnc 代表人，dmel 代表苍蝇 ，我们这儿pbmc3k为人的，选择hgnc
dbDir <- "cisTarget_databases" # RcisTarget databases location

myDatasetTitle <- "FCD scRNA-seq(normal_16)" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=8)  #nCores=10，代表开启10个线程计算
scenicOptions@inputDatasetInfo$cellInfo <- "./int/cellInfo.Rds"

#若未设定颜色文件则跳过下面一行代码
#scenicOptions@inputDatasetInfo$colVars <- "./int/colVars.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

###################################################       构建Co-expression network

exprMat <- as.matrix(pbmc@assays$RNA@counts) #准备表达矩阵，为了节省计算资源，随机只抽取部分细胞，计算共表达网络

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=3*.01*ncol(exprMat),minSamples=ncol(exprMat)*.01) #基因过滤/选择，去除最有可能是噪音的基因
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions) ##计算相关性矩阵，1.2_corrMat.Rds：基因之间的相关性矩阵

###############################################GENIE3，GENIE3 的输入通常是一个表达式矩阵和一个候选调节器列表
exprMat_filtered <- log2(exprMat_filtered+1)  #完整矩阵

#导出用于GRNboost(in python)的输入文件
#write.csv(t(exprMat_filtered),"exprMat_filtered.csv")
exportsForArboreto(exprMat, scenicOptions, dir = ".")

runGenie3(exprMat_filtered, scenicOptions, nParts = 10,resumePreviousRun=TRUE) #nParts参数，是把表达矩阵分成n份分开计算, !!!!!!这一步时间极其长！！！！！！！
#1.4，GENIE3_linkList.Rds：GENIE3最终结果，是把“1.3_”开头的文件合并在一起
#exportsForArboreto(exprMat, scenicOptions, dir = "int")
#GRNBoost



#导入GRNboost(in python)的输出文件，并转为R包可识别的RDS文件，修改weight值的值域
GRNBoost_output <- read.delim("E:/Linux share/output.tsv", header=TRUE)
colnames(GRNBoost_output) <- c("TF","Target","weight")
class(GRNBoost_output[,"weight"])
GRNBoost_output[,"weight"] = GRNBoost_output[,"weight"]/max(GRNBoost_output[,"weight"])
#GRNBoost_output[,"weight"] = factor(GRNBoost_output[,"weight"])

class(GRNBoost_output[,"weight"])
#class(GRNBoost_output[,"weight"]) = numeric
saveRDS(GRNBoost_output, file="G:/SCENIC/int/1.4_GENIE3_linkList.Rds")

head(GRNBoost_output)

#############################################
#运行到这儿，才进入SCENIC关键分析步骤
#1. 获取共表达模块
#2. 获取调节子（使用 RcisTarget）：TF 基序分析）
#3. 对细胞中的 GRN评分（使用 AUCell）
#4. 根据 GRN 活动对细胞进行聚类

exprMat_log <- log2(exprMat+1) #log标准话原始矩阵
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 8
scenicOptions@settings$seed <- 123
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) #1. 获取共表达模块

scenicOptions@settings$nCores = 1
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)  #2. 获取regulons
scenicOptions@settings$nCores = 8
#运行完runSCENIC_2_createRegulons后产生：
#output/Step2_MotifEnrichment.tsv
#output/Step2_MotifEnrichment_preview.html
#output/Step2_regulonTargetsInfo.tsv 文件





#可使用参数coexMethod=c("top5perTarget")，coexMethod可供选择，作者尝试了多种策略过滤低相关性TF-Target,并建议是6种过滤标准都用,其实也并没有耗时多少。默认都计算。
#w001：以每个TF为核心保留weight>0.001的基因形成共表达模块；#w005，#top50
#top5perTarget：每个基因保留weight值top5的TF得到精简的TF-Target关联表，然后把基因分配给TF构建共表达模块；#top10perTarget；#top50perTarget
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log) #3. 对细胞中的 GRN（调节子）评分
#运行完runSCENIC_3_scoreCells后，产生regulonAUC矩阵，可用热图展示，t-SNE图分别展示每个regulon的活性分布情况。


saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_log)
#runSCENIC_4_aucell_binarize，是进行 二进制转换及衍生分析 ，这步可不做





#####################################################################可视化
library(doParallel)
setwd("G:/SCENIC")
nPcs <- c(5,15,50)
scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
pdf("int/AUC_tsne.pdf",width=10,height=10)
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=TRUE, varName="manucelltype", cex=.5)
dev.off()

pbmc <-readRDS(RDSpathway)  #导入pbmc3k数据
#current.cluster.ids <- c(0:8)
#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
#pbmc@meta.data$celltype <- plyr::mapvalues(x = pbmc@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)
head(pbmc@meta.data)

##导入原始regulonAUC矩阵，也就是运行完runSCENIC_3_scoreCells后产生的AUC矩阵，可以查看每个GRN在每个细胞中的AUC活性打分

AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- data.frame(t(AUCmatrix@assays@data@listData$AUC), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC) #把(替换成_
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC) #把)去掉，最后如把 KLF3_extended (79g) 替换成KLF3_extended_79g
colnames(AUCmatrix) <- RegulonName_AUC
pbmcauc <- AddMetaData(pbmc, AUCmatrix) #把AUC矩阵添加到pbmc的metadata信息中
pbmcauc@assays$integrated <- NULL
saveRDS(pbmcauc,'pbmcauc.rds')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
pbmcbin <- AddMetaData(pbmc, BINmatrix)
pbmcbin@assays$integrated <- NULL
saveRDS(pbmcbin, 'pbmcbin.rds')

##利用Seurat可视化AUC
dir.create('scenic_seurat')
setwd("./scenic_seurat")
dir.create('featurePlot')
dir.create('vlnplot')
setwd("../")
#FeaturePlot
GRNs <-intersect(colnames(AUCmatrix),colnames(BINmatrix))
for(i in 1:length(GRNs)){
  p1 = FeaturePlot(pbmcauc, features=GRNs[i], label=T, reduction = 'umap')
  p2 = FeaturePlot(pbmcbin, features=GRNs[i], label=T, reduction = 'umap')
  p3 = DimPlot(pbmc, reduction = 'umap', group.by = "manucelltype", label=T)
  plotc = p1|p2|p3
  ggsave(paste("scenic_seurat/featurePlot/",GRNs[i],".png",sep=""), plotc, width=14 ,height=4)
}

#RidgePlot&VlnPlot ,和单细胞展示基因表达一样，这儿可以利用小提琴图展示GRNs的活性得分
for(i in 1:length(GRNs)){
  p1 = RidgePlot(pbmcauc, features =GRNs[i], group.by="celltype") + theme(legend.position='none')
  p2 = VlnPlot(pbmcauc, features =GRNs[i], pt.size = 0, group.by="celltype") + theme(legend.position='none')
  plotc = p1 + p2
  ggsave(paste("scenic_seurat/vlnplot/","Ridge-Vln_",GRNs[i],".png",sep=""),plotc, width=10, height=8)
}

#pySCENIC############################################################################
#0.loom 文件制备
   
library(Seurat)
sce = readRDS('G:/work/GMH/merge.rds')
table(Idents(sce))
p1=DimPlot(sce,label = T) 
p1

write.csv(t(as.matrix(sce@assays$RNA@counts)), file = "fibo_1000.csv")


#1.pySCENIC run in docker ........


#2.提取 out_SCENIC.loom 信息
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(SCENIC)

loom <- open_loom('G:/work/SCENIC/sample_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)


#3.加载SeuratData(可跳过)##########
library(SeuratData) #加载seurat数据集  
data("pbmc3k")  
sce <- pbmc3k.final   
table(sce$seurat_clusters)
table(Idents(sce))
sce$celltype = Idents(sce)

library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                   'CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB',  # mac
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'FCGR3A',
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'KLRB1','NCR1', # NK 
                   'FGF7','MME', 'ACTA2', ## fibo 
                   'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo 
                   'MKI67' , 'TOP2A', 
                   'PECAM1', 'VWF',  ## endo 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check
p <- DotPlot(sce , features = unique(genes_to_check),assay='RNA'  )  + coord_flip() +   theme(axis.text.x=element_text(angle=45,hjust = 1))

p 
ggsave('check_last_markers.pdf',height = 11,width = 11)
DimPlot(sce,reduction = "umap",label=T ) 
sce$sub_celltype =  sce$celltype
DimPlot(sce,reduction = "umap",label=T,group.by = "sub_celltype" )
ggsave('umap-by-sub_celltype.pdf')
Idents(sce) <- sce$sub_celltype

sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.8)
table(sce@meta.data$RNA_snn_res.0.8)  
DimPlot(sce,reduction = "umap",label=T ) 
ggsave('umap-by-sub_RNA_snn_res.0.8.pdf')


#4.四种可视化
sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
dim(sub_regulonAUC)
sce 
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(sce))
#[1] TRUE

cellClusters <- data.frame(row.names = colnames(sce), seurat_clusters = as.character(sce$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(sce), celltype = sce$sub_celltype)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 
save(sub_regulonAUC,cellTypes,cellClusters,sce,file = 'for_rss_and_visual.Rdata')

#尴尬的是TCF4(+)我这个数据里面没有，换了个PAX5(+)和REL(+)
regulonsToPlot = c('PAX5(+)','REL(+)')
regulonsToPlot
sce@meta.data = cbind(sce@meta.data ,t(assay(   sub_regulonAUC[regulonsToPlot,])))
Idents(sce) <- sce$sub_celltype
table(Idents(sce))

DotPlot(sce, features = unique(regulonsToPlot)) + RotatedAxis()
RidgePlot(sce, features = regulonsToPlot , ncol = 1)
VlnPlot(sce, features = regulonsToPlot,pt.size = 0 ) 
FeaturePlot(sce, features = regulonsToPlot )


#5.热图查看TF分布
pheatmap(regulonActivity_byGroup_Scaled)


#6.rss查看特异TF
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)


#7.其他查看TF方式
library(dplyr) 
rss=regulonActivity_byGroup_Scaled
head(rss)
library(dplyr) 
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 
#rownames(rowcn) = rownames(n)
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T)


#4.四种可视化
sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
dim(sub_regulonAUC)
sce 
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(sce))
#[1] TRUE

cellClusters <- data.frame(row.names = colnames(sce), seurat_clusters = as.character(sce$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(sce), celltype = as.character(Idents(sce)))
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 
#save(sub_regulonAUC,cellTypes,cellClusters,sce,file = 'for_rss_and_visual.Rdata')
dim(sub_regulonAUC)

#5.TF
# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), cellTypes[,selectedResolution]) 

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)
#[1]  220 2638 #似乎没啥区别

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,function(cells) rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression. 
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
# 参考：https://www.jianshu.com/p/115d07af3029
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)
#[1] 220   9
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)



#热图查看TF分布
pheatmap(regulonActivity_byGroup_Scaled)


#rss查看特异TF
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)


#其他查看TF方式
library(dplyr) 
rss=regulonActivity_byGroup_Scaled
head(rss)
library(dplyr) 
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  #median、sd、mean
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 
#rownames(rowcn) = rownames(n)
pheatmap(n,annotation_row = rowcn,show_rownames = T)



