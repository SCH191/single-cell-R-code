library(reticulate)
use_python("C:/Users/11723/AppData/Local/Programs/Python/Python39/python.exe")


library(CytoTRACE) #加载包

#results <- CytoTRACE(marrow_10x_expr)#使用自带骨髓 10x scRNA-seq 数据集运行 CytoTRACE

# 或者使用“ncores”（默认值 = 1）进行多线程
#或使用“subsamplingsize”（默认值 = 1,000 个单元格）指示子采样大小
#比如：使用 8 个内核和 1,000 个子样本在快速模式下运行以下数据集
#Windows不可用
#results <- CytoTRACE(marrow_10x_expr, ncores = 8, subsamplesize = 1000)

#在多个 scRNA-seq 批次/数据集上运行
#datasets <- list(marrow_10x_expr, marrow_plate_expr)
#results <- iCytoTRACE(datasets)

ABC

results <- CytoTRACE(as.data.frame(ABC@assays$RNA@counts))

#绘制 CytoTRACE 和 iCytoTRACE 结果

#可视化 CytoTRACE 结果
#plotCytoTRACE(results, phenotype = marrow_10x_pheno)

types_manu = as.character(ABC@meta.data$manusubcelltype)
names(types_manu) = as.character(rownames(ABC@meta.data))
plotCytoTRACE(results,emb = ABC@reductions$umap@cell.embeddings,phenotype = types_manu)

#可视化与 CytoTRACE 相关的基因
plotCytoGenes(results, numOfGenes = 10)