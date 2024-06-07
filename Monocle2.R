library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)

setwd("G:/work/test/新建文件夹/micro 2")
ABC = readRDS('Microglia.rds')
DefaultAssay(ABC) = 'RNA'
pbmc = ABC #导入注释好的seurat对象（已注释）

##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix = as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data = pbmc@meta.data 
p_data$celltype = pbmc@active.ident  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
##提取基因信息 如生物类型、gc含量等
f_data = data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)

#构建CDS对象
pd = new('AnnotatedDataFrame', data = p_data) 
fd = new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds = newCellDataSet(expr_matrix,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())


cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

#######################################
##使用seurat选择的高变基因⚠️
express_genes = VariableFeatures(pbmc)
cds = setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

##使用clusters差异表达基因
deg.cluster = FindAllMarkers(pbmc)
express_genes = subset(deg.cluster,p_val_adj<0.05)$gene
cds = setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

##使用monocle选择的高变基因⚠️
disp_table = dispersionTable(cds)
disp.genes = subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds = setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)





#理想情况下，我们希望尽可能少地使用正在研究的系统生物学的先验知识。我们希望从数据中发现重要的排序基因，而不是依赖于文献和教科书，因为这可能会在排序中引入偏见。我们将从一种更简单的方法开始，但是我们通常推荐一种更复杂的方法，称为“dpFeature”。
#这一步输入的expressed_genes来自于步骤4。
#⚠️⚠️后续分析使用的是该方法
expressed_genes = VariableFeatures(pbmc) 
diff = differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~celltype",cores=1) 
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)

##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg = subset(diff, qval < 0.01) #选出2829个基因
deg = deg[order(deg$qval,decreasing=F),]
head(deg)

##差异基因的结果文件保存
#write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)

## 轨迹构建基因可视化
ordergene = rownames(deg) 
cds = setOrderingFilter(cds, ordergene)  
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
#pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
#dev.off()
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)



######################################

cds <- reduceDimension(cds, 
                       max_components = 2, 
                       num_dim = 6, 
                       reduction_method = 'DDRTree', 
                       residualModelFormulaStr = "~orig.ident", #去除样本影响
                       verbose = F)
cds= orderCells(cds)
cds= orderCells(cds, root_state = 5)
#⚠️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds = orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点

####################################

plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 

plot_cell_trajectory(cds,color_by="manusubcelltype2", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "manusubcelltype2") + facet_wrap("~manusubcelltype2", nrow = 1)

plot_cell_trajectory(cds,color_by="State", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)

#
p1 = plot_cell_trajectory(cds, x = 1, y = 2, color_by = "manusubcelltype2") + 
  theme(legend.position='none',panel.border = element_blank())# + #去掉第一个的legend
  #scale_color_manual(values = colour) 
p2 = plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                      color_by = "manusubcelltype2")+
  #scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) 
p1|p2

#
library(ggpubr)
df = pData(cds) 
## pData(cds)取出的是cds对象中cds@phenoData@data的内容
View(df)
ggplot(df, aes(Pseudotime, colour = manusubcelltype2, fill=manusubcelltype2)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()

#
#指定基因
s.genes = c('P2RY12','CCL4','JUN','CD83','IL1B','APOE','APOC1')
p1 = plot_genes_jitter(cds[s.genes,], grouping = "manusubcelltype2", color_by = "manusubcelltype2")
p2 = plot_genes_violin(cds[s.genes,], grouping = "manusubcelltype2", color_by = "manusubcelltype2")
p3 = plot_genes_in_pseudotime(cds[s.genes,], color_by = "manusubcelltype2")
plotc = p1|p2|p3
plotc



