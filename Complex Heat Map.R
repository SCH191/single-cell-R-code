library(ComplexHeatmap)
library(dendextend)
library(magick)
library(dplyr)


#future::plan("multiprocess", workers = 6);options(future.globals.maxSize = 100000 * 1024^5) #设置任务多线程
##Data process >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#参数##################################################################################

Heatmap(matrix, col, name, 
        #matrix：数字或字符型矩阵（可以是离散或连续型数值）
        #col:定义热图颜色，对离散型数据，col可以是一个向量；对连续型数据，col可以是一个函数，也可以用colorRamp2 函数生成
        #name:热图图例名称
        na_col = 'grey',
        #ComplexHeatmap允许数据中含有NA,需要通过参数na_col来控制NA的颜色
        color_space = 'LAB',
        #当矩阵是数值型矩阵，col是一个向量时，控制内插颜色
        rect_gp = gpar(col = NA),
        #rect_gp:热图体区矩形的参数，如设置举行边框为白色
        cell_fun = NULL,
        #cell_fun：自定义在cell中增加绘图项的函数。7个参数：i(row index,矩阵中的行index）, j(column index，矩阵中的列index), x,y(热图体区中中间点的坐标）,width,height(cell的宽度和高度）,fill(cell的填充颜色）
        row_title = character(0),
        #row_title：行标题
        row_title_side = c('left', 'right'),
        #row_title_side：行标题位置，左('left')，右('right')
        row_title_gp = gpar(fontsize = 14),
        #row_title_gp:设置行标题的文本属性，此处为字体大小为14
        row_title_rot = switch(row_title_side[1], 'left' = 90, 'right' = 270),
        #row_title_rot:行标题的旋转角度，可选为0,90,270
        column_title = character(0),
        #column_title：列标题
        column_title_side = c('top', 'bottom'),
        #column_title_side：列标题位置，上('top')，下('bottom')
        column_title_gp = gpar(fontsize = 14),
        #column_title_gp：设置列标题的文本属性
        column_title_rot = 0,
        #column_title_rot：列标题的旋转角度，可选为0,90,270
        cluster_rows = TRUE,
        #cluster_rows:是否行聚类
        clustering_distance_rows = 'euclidean',
        #clustering_distance_rows：行聚类的距离方法，默认为'euclidean'，也可以为自定义函数
        clustering_method_rows = 'complete',
        #clustering_method_rows:行聚类的方法，默认为'complete'，可参考hclust
        row_dend_side = c('left', 'right'),
        #row_dend_side：行聚类树位置，左('left')，右('right')
        row_dend_width = unit(10, 'mm'),
        #row_dend_width：行聚类树的宽度，unit对象
        show_row_dend = TRUE,
        #show_row_dend：是否展示行聚类树
        row_dend_reorder = TRUE,
        #row_dend_reorder:对行重新排序，该值可以是逻辑值或包含用于重新排序行的权重的向量
        row_dend_gp = gpar(),
        #row_dend_gp：绘图线的图形参数。如果已经提供了带有边渲染的树形图对象，则该参数将被忽略。
        row_hclust_side = row_dend_side,
        #row_hclust_side：已弃用
        row_hclust_width = row_dend_width,
        #row_hclust_width：已弃用
        show_row_hclust = show_row_dend,
        #show_row_hclust：已弃用
        row_hclust_reorder = row_dend_reorder,
        #row_hclust_reorder：已弃用
        row_hclust_gp = row_dend_gp,
        #row_hclust_gp：已弃用
        cluster_columns = TRUE,
        #cluster_columns：是否列聚类
        clustering_distance_columns = 'euclidean',
        #clustering_distance_columns：列聚类的距离方法，也可以为自定义函数
        clustering_method_columns = 'complete',
        #clustering_method_columns：列聚类方法，可参考hclust
        column_dend_side = c('top', 'bottom'),
        #column_dend_side：列聚类树位置，上('top')，下('bottom')
        column_dend_height = unit(10, 'mm'),
        #column_dend_height:行聚类树的高度，unit对象
        show_column_dend = TRUE,
        #show_column_dend:是否展示列聚类树
        column_dend_gp = gpar(),
        #column_dend_gp:绘图线的图形参数。如果已经提供了带有边渲染的树形图对象，则该参数将被忽略。
        column_dend_reorder = TRUE,
        #column_dend_reorder：对列重新排序，该值可以是逻辑值或包含用于重新排序列的权重的向量
        column_hclust_side = column_dend_side,
        #column_hclust_side：已弃用
        column_hclust_height = column_dend_height,
        #column_hclust_height：已弃用
        show_column_hclust = show_column_dend,
        #show_column_hclust：已弃用
        column_hclust_gp = column_dend_gp,
        #column_hclust_gp：已弃用
        column_hclust_reorder = column_dend_reorder,
        #column_hclust_reorder：已弃用
        row_order = NULL,
        #row_order：行的顺序。如果选择此热图作为主热图，则可以轻松调整热图列表的行顺序。手动设置行顺序应关闭群集
        column_order = NULL,
        #column_order：列的顺序。它可以轻松调整矩阵和列注释的列顺序
        row_names_side = c('right', 'left'),
        #row_names_side：行名称位置。
        show_row_names = TRUE,
        #show_row_names：是否展示行名称
        row_names_max_width = default_row_names_max_width(),
        #row_names_max_width：行名称的最大宽度。因为某些时候行名称可能很长，所以显示它们都是不合理的。
        row_names_gp = gpar(fontsize = 12),
        #row_names_gp：行名称文本属性
        column_names_side = c('bottom', 'top'),
        #column_names_side：列名称位置
        show_column_names = TRUE,
        #show_column_names：是否展示列名称
        column_names_max_height = default_column_names_max_height(),
        #column_names_max_height：行名称的最大宽度。
        column_names_gp = gpar(fontsize = 12),
        #column_names_gp：列名称文本属性
        top_annotation = new('HeatmapAnnotation'),
        #top_annotation：用HeatmapAnnotation函数构建的注释对象，在顶部添加注释信息
        top_annotation_height = top_annotation@size,
        #top_annotation_height：顶部注释信息展示的总高度
        bottom_annotation = new('HeatmapAnnotation'),
        #bottom_annotation：用HeatmapAnnotation函数构建的底部注释对象
        bottom_annotation_height = bottom_annotation@size,
        #bottom_annotation_height：底部注释信息展示的总高度
        km = 1,
        #km：对行做k-means聚类的类数，若k>1,热图会根据k-means聚类对行进行分裂,对每个cluster,进行层次聚类
        km_title = 'cluster%i',
        #km_title：设置km时每个cluster的行标题。它必须是格式为'。*％i。*'的文本，其中'％i'由cluster的索引替换
        split = NULL,
        #split:行按照split定义的向量或者数据框进行分裂。但是，如果cluster_rows是聚类对象，则split可以是单个数字，表示将根据树上的拆分来拆分行
        gap = unit(1, 'mm'),
        #gap:如果热图按行分割，则行切片之间的间隙应为单位对象。如果是矢量，则热图中的顺序对应于从上到下
        combined_name_fun = function(x) paste(x, collapse = '/'),
        #combined_name_fun:如果热图按行分割，如何为每个切片创建组合行标题？ 此函数的输入参数是一个向量，它包含split中每列下的级别名称。
        width = NULL,
        #width:单个热图的宽度应该是固定的单位对象。 当热图附加到热图列表时，它用于布局。
        show_heatmap_legend = TRUE,
        #show_heatmap_legend:是否展示图例
        heatmap_legend_param = list(title = name),
        #heatmap_legend_param：热图图例设置（标题，位置，方向，高度等）参数列表，详情可见color_mapping_legend，ColorMapping-method。例如：heatmap_legend_param = list(title= 'legend', title_position ='topcenter',
        legend_height=unit(8,'cm'),legend_direction='vertical',
use_raster = FALSE,
#use_raster：是否将热图图像体渲染为光栅图像。当矩阵很大时，它有助于减小文件大小。如果设置了cell_fun，则强制use_raster为FALSE
raster_device = c('png', 'jpeg', 'tiff', 'CairoPNG', 'CairoJPEG', 'CairoTIFF'),
#raster_device：用于生成光栅图像的图形设备
raster_quality = 2,
#raster_quality：设置为大于1的值将改善光栅图像的质量。
raster_device_param = list()
#raster_device_param：所选图形设备的其他参数列表。
)


ABC@meta.data$manucelltype = ordered(ABC@meta.data$manucelltype,levels = c('Astrocyte','OPC','Oligo','Microglia','Fibroblasts','Epithelial_cell','Endothelial_cell','Smooth_muscle_cells','NK','Neutrophlis','Monocytes'))
############################
scRNA = ABC  #加载Seurat对象数据
### 计算seurat_cluster间的 差异基因
markers <- FindAllMarkers(object = scRNA, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
### 挑选每个cluster top5 的基因画表达热图
select.features <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) 
### 随机为 细胞样本分配一个模拟标签
#scRNA@meta.data$sample <- sample(c("species.1","species.2","species.3","species.4"),size = ncol(scRNA),replace = T) 
scRNA@meta.data$sample = paste(scRNA@meta.data$Group,'_',scRNA@meta.data$manucelltype,sep = '')
### 提取Seurat绘制热图的矩阵数据
data <- Seurat::DoHeatmap(scRNA, features = select.features$gene, group.by = "seurat_clusters", group.bar = T, size = 4)$data 
### 提取细胞标识

cell.meta <- scRNA@meta.data %>% tibble::rownames_to_column(var="Cell") %>% dplyr::select(Cell,seurat_clusters,sample) %>% arrange(seurat_clusters) %>% mutate(Identity = seurat_clusters)

### 长数据转换宽数据
#counts <- data %>% dplyr::select(!Identity) %>% tidyr::pivot_wider(names_from =  Cell, values_from = Expression) %>% data.frame() %>% tibble::column_to_rownames(var = "Feature") %>% dplyr::select(cell.meta$Cell)
#mat <- as.matrix(counts) #格式转换为矩阵
tempt <- ScaleData(scRNA, features = select.features$gene)
mat = GetAssayData(tempt,slot = 'scale.data')


colormaps <- c(RColorBrewer::brewer.pal(name="Dark2", n = 8),RColorBrewer::brewer.pal(name="Pastel1", n = 9),RColorBrewer::brewer.pal(name="Paired", n = 12),RColorBrewer::brewer.pal(name="Set1", n = 9),RColorBrewer::brewer.pal(name="Pastel2", n = 8),RColorBrewer::brewer.pal(name="Set3", n = 12))
scales::show_col(colormaps)  


#ComPlexHeatMap############################################################################
# Dg_tree
### 创建表达矩阵的样本组间聚类树【计算组内样本均值进行建树】
dend1 = cluster_between_groups(mat, cell.meta$Identity)
### 你希望这些样本被聚类成几簇 【按树枝颜色区分】
dend1 = color_branches(dend1, k = 8)  
### 树样式调整
dend1 = dend1 %>% set("branches_lwd", 1.5) # 聚类树树枝线条 厚度
### dend1 = dend1 %>% raise.dendrogram (3) #聚类树底端线条厚度
### dend1 = dend1 %>% highlight_branches_col(viridis::viridis(100)) #聚类树颜色调整
### dend1 = dend1 %>% highlight_branches_col(rev(viridis::magma(1000)))  #聚类树颜色调整


# 注释1 empty
ha_top_1 <- HeatmapAnnotation(
  empty = anno_empty(border = FALSE,height = unit(0.1, "cm")), #添加空的注释块
  annotation_name_side = "left",which = "column" 
)



# 注释2 Group
### 获取树聚类后的矩阵样本的排列顺序
HM <- Heatmap(mat,cluster_columns = dend1)
HM = draw(HM)
### 根据树聚类的样本排列顺序 来排列细胞信息表cell.meta
group.data <- cell.meta[column_order(HM),] 
### 提取按树聚类排布的样本簇标签顺序
group_order_label <- unique(group.data$Identity,fromLast = F) %>% as.vector()

### 创建样本簇标识 色板
color.cl <- colormaps[seq(length(unique(cell.meta$Identity)))] 
### 创建簇标识色块注释对象
ha_top_2 <- HeatmapAnnotation( 
  Group = anno_block(gp = gpar(fill = color.cl,col = 0),
                     labels = group_order_label, #块注释标签
                     labels_gp = gpar(col = "black",fontsize = 16, fontface = "bold") , #注释文本样式
                     show_name = TRUE , #显示注释对象名
                     height = unit(0.5,"cm") # 注释对象的整体高度
                     # weight = unit(10,"cm") # 注释对象的整体宽度
  ),
  annotation_name_side = "left",#注释对象名显示方向
  which = "column"
)
### anno_block 块注释的图例构建
lgd_Group <- Legend(title = "Group", labels = group_order_label,legend_gp = gpar(fill = color.cl)) 




# 注释3 Batch
ha_top_3 <- HeatmapAnnotation(Batch = cell.meta$sample, 
                              annotation_legend_param = list(Batch = list(title = "Batch",ncol=1)), #注释图例参数调整
                              annotation_name_side = "left",which = "column")




### 统计每个基因的表达量
sum_Normexpr <- scRNA@assays$RNA@data[rownames(mat),] %>% Matrix::rowSums()
ha_rig = rowAnnotation(sum_Normexpr = anno_barplot(sum_Normexpr, bar_width = 1,gp = gpar(fill = "yellow",col="red"),
                                                   border=F, #行注释对象外侧边框
                                                   width = unit(2,"cm"), # 行注释的宽度
                                                   axis_param =(list(side = "top",gp=gpar(fontsize=5,col="red"))) # 坐标轴参数
),
show_annotation_name = FALSE, #不显示注释对象标题
annotation_name_side = "bottom",# 注释标题旋转位置
annotation_name_gp= gpar(fontsize = 8), #注释标题大小
annotation_name_rot = 0 #注释标题旋转
)
### anno_block注释图例对象创建
lgd_sumExpr <- Legend(title = "sum_Norm_expr",at = "",legend_gp = gpar(fill = "yellow")) 





Heatmap(mat,
        cluster_columns = dend1, #列方向添加 簇级 树聚类
        column_split = length(unique(cell.meta$Identity)), #热图列方向按簇拆分
        #热图主体
        column_dend_height = unit(5, "cm"), #树的高度
        clustering_method_columns = "spearson", #树的聚类方法
        column_title = "_OH_MY_Doheatmap_", #列方向大标题
        column_title_side = "bottom",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"), #列方向大标题样式
        
        name = "Expr", #热图名称，表达量图例名
        cluster_rows = FALSE, #关闭行方向聚类
        show_column_names = FALSE, #关闭显示列名
        show_row_names = TRUE, #打开显示行名
        col = viridis::viridis(200), #表达量梯度颜色设置
        na_col = "black", #空值单元格的颜色
        row_title = "cluster_between_groups", #行方向大标题
        row_title_gp = grid::gpar(fontsize = 20,fontface="bold"), #行方向大标题样式
        row_names_side = "left", #行名显示方向
        row_names_gp = grid::gpar(fontsize = 6,fontface="bold"), #行名大小调整
        border = TRUE, #热图图像外边框显示
        
        # 表达量图例 样式设置
        heatmap_legend_param = list(
          title = "Exp",
          border = "red",
          direction = "vertical",
          title_position = "topleft"
          # legend_height = unit(12, "cm") # 热图表达量图例大小
        ),
        # 顶部注释
        top_annotation = c(ha_top_1,ha_top_2,ha_top_3), # 合并多个注释对象
        # 右注释
        right_annotation = ha_rig,
        
        ##图像 光栅化转换
        use_raster = TRUE, raster_quality = 5
) %>% draw(merge_legend = TRUE,padding = unit(c(1, 1, 2, 1), "cm"), # panding:图像编剧下-左-上-右
           annotation_legend_list = list(lgd_Group,lgd_sumExpr) # 添加 自己创建的 legend 对象
)
decorate_column_dend("Expr", {grid.yaxis()}) # 树聚类 修饰

####
png("Complex Heat Map.png",width = 50000,height = 5000,res = 216)
Heatmap(mat,
        cluster_columns = dend1, #列方向添加 簇级 树聚类
        column_split = length(unique(cell.meta$Identity)), #热图列方向按簇拆分
        #热图主体
        column_dend_height = unit(10, "cm"), #树的高度
        clustering_method_columns = "spearson", #树的聚类方法
        column_title = "_OH_MY_Doheatmap_", #列方向大标题
        column_title_side = "bottom",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"), #列方向大标题样式
        
        name = "Expr", #热图名称，表达量图例名
        cluster_rows = FALSE, #关闭行方向聚类
        show_column_names = FALSE, #关闭显示列名
        show_row_names = TRUE, #打开显示行名
        col = viridis::viridis(200), #表达量梯度颜色设置
        na_col = "black", #空值单元格的颜色
        row_title = "cluster_between_groups", #行方向大标题
        row_title_gp = grid::gpar(fontsize = 20,fontface="bold"), #行方向大标题样式
        row_names_side = "left", #行名显示方向
        row_names_gp = grid::gpar(fontsize = 6,fontface="bold"), #行名大小调整
        border = TRUE, #热图图像外边框显示
        
        # 表达量图例 样式设置
        heatmap_legend_param = list(
          title = "Exp",
          border = "red",
          direction = "vertical",
          title_position = "topleft"
          # legend_height = unit(12, "cm") # 热图表达量图例大小
        ),
        # 顶部注释
        top_annotation = c(ha_top_1,ha_top_2,ha_top_3), # 合并多个注释对象
        # 右注释
        right_annotation = ha_rig,
        
        ##图像 光栅化转换
        use_raster = TRUE, raster_quality = 5
) %>% draw(merge_legend = TRUE,padding = unit(c(1, 1, 2, 1), "cm"), # panding:图像编剧下-左-上-右
           annotation_legend_list = list(lgd_Group,lgd_sumExpr) # 添加 自己创建的 legend 对象
)
decorate_column_dend("Expr", {grid.yaxis()}) # 树聚类 修饰
dev.off()

write.csv(markers,'All_diff.csv')
write.csv(select.features,'Top5_all_diff.csv')


####################
library(gplots)
library(pheatmap)

markers <- FindAllMarkers(object = ABC, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
tempt <- ScaleData(ABC, features = top10$gene)
mat = GetAssayData(tempt,slot = 'scale.data')

annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(heatdata)

heatmap.2(mat,trace = "none")

pheatmap(mat, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = F,# 显示行名
         show_colnames = F,# 显示列名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("blue", "white","red"))(100))
####################
markers <- FindAllMarkers(object = ABC, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(ABC, features = top10$gene)





