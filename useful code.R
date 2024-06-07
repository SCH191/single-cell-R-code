library(Seurat)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(patchwork)
library(CellChat)
library(ggalluvial)
library(svglite)
library(monocle3)
library(tidyverse)
library(NMF)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(GSVA)
library(pheatmap)
library(fgsea)
library(tibble)
library(SeuratWrappers)
library(gridExtra)
library(Matrix)
library(ggsci)
library(GOplot)
library(scRNAtoolVis)
library(DoubletFinder)
library('fanyi')
library('beepr')
R.utils::setOption( "clusterProfiler.download.method",'auto' )
set_translate_option('20240323002002720','xUMChTWH064IoTV16nke',source = 'baidu')

###########

for (i in length(a)) {
  for (j in height) {
    
  }
}


for (i in 1:length(brain)) {
test = brain[[i]]
markers <- FindAllMarkers(object = test, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
data = rbind(data,top10)
write.csv(top10,paste(names(brain[i]),'.csv',sep = ''))
}


merge = RenameIdents(merge,'HS1.1'='','HS1.2'='','HS2.1'='','HS2.2'='','HS3.1'='','HS3.2'='','HS4.1'='','HS4.2'='','HS5.1'='','HS5.2'='','HS6.1'='','HS6.2'='','HS7.1'='','HS7.2'='','HS8.1'='','HS8.2'='')


j = merge
marker = c('RBFOX3','SLC17A7','GAD1','GAD2','CUX2','RORB','THEMIS','FEZF2','PVALB','SST','VIP','ID2')#神经元亚型
marker = c('ASCL1','NEUROD2','GAD1','PROX1','MEIS2','LHX6','NR2F2','AQP4','OLIG2','MBP','PTPRC','SPARC')#海马部位分型
marker = c('ASC','CASP1','CASP4','CASP5','CASP11','GSDMD','IL1B','IL18','IL33','AIM2')#细胞焦亡与凋亡对比
marker = c('CASP2','CASP7','CASP10','PARP1')#细胞焦亡与凋亡对比
marker = c('GPX4','GCL','GSS','FSP1','ACSL4','LPCAT3','ALOX15')#铁死亡
marker = c('PARP1')#Parthanatos
marker = c('SLC1A2','SLC1A3','GLUL')#谷氨酸转运蛋白、谷氨酰胺合酶
marker = c('ID2','RELN','NRG1','VIP','CBLN1','CALB2','SST','HPGD','TRHDE','PVALB','MEPE','SULF1','PAX6','NEFH','NEFL','NEFM')
marker = c('RBFOX3','GAD1','GAD2','ID2','LAMP5','RELN','PVALB','SULF1','LHX6','SST','SEMA6A','NOS1','VIP','TAC3','CALB2')
#In
marker = c('CUX2','PDGFD','LAMA2','RORB','COBLL1','PLCH1','THEMIS','INPP4B','NR4A2','LRRK1','SEMA3E','IL26','TLE4','NEFH','NEFL','NEFM')
#Ex
marker = c('Pik3ca','Pik3cb','Pik3cd','Akt','Bcl2','Casp8','Casp10','Jnk','Gzmb','Map3k5')#mouse凋亡
marker = c('SOX10','OLIG1','OLIG2','PLP1','MYRF','MOG','MAG','PTPRZ1','PDGFRA','VCAN','BMP4','GPR17','BCAS1','RRAS2','PROM1','THBS3','SLC9A3R2','FOS','KLK6','HOPX','PTGDS','IL33','SERPINA3N','C4B','CDKN1A','TNFRSF12A','IRF9','IFIT1')
#OL系细胞亚型区分（mouse参考）

marker = c('Ear2','Gngt2','Car2','Hba-a1','Ms4a7','C1qa','Apoe','H2-Aa','Lgals3','Ly6c2','Ms4a6c','Ctss','Gata2','Cpa3','Ngp','Camp','S100a9','S100a8','Gata3','Siglech','Klrd1','Nkg7','Il7r','Cd4','Cd8b1','Cd3g','Cd79b','Cd79a','Ifitm1','Angpt1','Adgrl4','Cd34','Kit')
#sham
marker = c('AQP4','C1QL1','BCAS1','MOG','GAD2','RGS4','PROX1','FLT1','DCN','C1QA','LYVE1','RELN','CFAP126','SLC1A3','HMGB2','ASCL1','HMGB2','SOX4','RFC4','NNAT','NPY1R')
#海马 “https://www.nature.com/articles/s41593-022-01073-x#Sec12”

####
list = 'p1'
for (i in 1:length(marker)){
  eval(parse(text = paste("p",i,"=VlnPlot(ABC, features =marker[i],ncol=1,pt.size = 0.0)+geom_boxplot(width = 0.1, outlier.shape = NA,col = 'black', fill='white')+ NoLegend()",sep = '')))
  if(i != 1){
    list = paste(list,'|p',i,sep='')
  }
}
eval(parse(list))

marker = c()

plots_violins = list()
for (i in 1:length(marker)){
  plots_violins[[i]] = VlnPlot(ABC, features =marker[i],cols = c('#327eba','#e06663'),ncol=1,pt.size = 0.0)+geom_boxplot(width = 0.1, outlier.shape = NA,col = 'black', fill='white')+ NoLegend()
}
CombinePlots(plots_violins,ncol = 3)
####
VlnPlot(ABC, features =marker,ncol=min(length(marker),3),pt.size = 0.0,raster=FALSE) #+geom_boxplot(width = 0.2,col = 'black', fill='white')+ NoLegend()
FeaturePlot(ABC, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "#e06663"),ncol=min(length(marker),5),reduction = 'umap',raster=FALSE)
DotPlot(ABC, features =marker,cols = c('#327eba','#e06663'))  + theme(axis.text.x = element_text(face = "italic",angle=90, hjust=1)) #+ coord_flip()
RidgePlot(ABC, features = marker, ncol = 4)
DoHeatmap(ABC, features = marker) + NoLegend()

ABC@meta.data$test = paste(ABC@meta.data$manusubcelltype,ABC@meta.data$Group,sep='__')
ABC@meta.data$test = factor(ABC@meta.data$test,levels = sort(levels(as.factor(ABC@meta.data$test))))
Idents(ABC) = ABC@meta.data$test
DefaultAssay(ABC) = 'RNA'

top = 20
i=3
j = brain[[i]]
marker = c('Pik3cd','Bcl2','Casp8','Map3k5')#diaowang
tiff(paste(names(brain[i]),'diao wang_gene_VlnPlot.tiff',sep = ''),width = 6000,height = 200*top,units = 'px',res = 450)
VlnPlot(j, features =marker,ncol=2,pt.size = 0.0)
dev.off()
tiff(paste(names(brain[i]),'diao wang_gene_DotPlot.tiff',sep = ''),width = 6000,height = 200*top,units = 'px',res = 450)
DotPlot(j, features =marker,dot.scale = 15,cols = c('blue','red'))
dev.off()
tiff(paste(names(brain[i]),'diao wang_gene_RiderPlot.tiff',sep = ''),width = 6000,height = 200*top,units = 'px',res = 450)
RidgePlot(j, features = marker, ncol = 4)
dev.off()
marker = c('Asc','Casp1','Casp4','Casp5','Casp11','Gsdmd')#细胞焦亡与凋亡对比
tiff(paste(names(brain[i]),'jiao wang_gene_VlnPlot.tiff',sep = ''),width = 6000,height = 200*top,units = 'px',res = 450)
VlnPlot(j, features =marker,ncol=2,pt.size = 0.0)
dev.off()
tiff(paste(names(brain[i]),'jiao wang_gene_DotPlot.tiff',sep = ''),width = 6000,height = 200*top,units = 'px',res = 450)
DotPlot(j, features =marker,dot.scale = 15,cols = c('blue','red'))
dev.off()
tiff(paste(names(brain[i]),'jiao wang_gene_RiderPlot.tiff',sep = ''),width = 6000,height = 200*top,units = 'px',res = 450)
RidgePlot(j, features = marker, ncol = 4)
dev.off()



ggsave('1.PNG',RidgePlot(j, features = marker, ncol = 2),width = 20,height = 40)

pdf(paste("out",name,".pdf",sep=""),width=18,height=12)
marker=c('NEFH','NEFL','NEFM','MAP1B','MTOR','RHEB','PIK3R3')
plot(c(1:10),main="Dysplastic Neuron,DN")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('VIM','NES','BCL2','PROM1','MAP2','C3','CLU','EFEMP1','SPARC','SERPINA3')
plot(c(1:10),main="Balloon Cell,BC")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
dev.off()



pdf(paste("out",name,".pdf",sep=""),width=18,height=12)

marker=c('RBFOX3','NF1','SYN1','MAP2')
plot(c(1:10),main="Neuron")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('SLC17A7','SATB2')
plot(c(1:10),main="Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('LAMP5','CUX2')
plot(c(1:10),main="L2_3_CUX2")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('GABRG1','Rorb')
plot(c(1:10),main="L4_Rorb")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('SEMA3E','NR4A2','CDH10','SLIT3','FEZF2','KCNIP1')
plot(c(1:10),main="L5_6_Fezf2")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('NR4A2','GRIK3','NTNG2','THEMIS')
plot(c(1:10),main="L5_6_Themis")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('GAD1','GAD2','CALB2','LHX6','SLC32A1','NXPH1')
plot(c(1:10),main="Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('PVALB','MEPE')
plot(c(1:10),main="PVALB")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('SST','RELN','GRIK1','TAC3','NPY')
plot(c(1:10),main="SST")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('VIP','VALB2','TAC3')
plot(c(1:10),main="VIP")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('NEFH','NEFL','NEFM','MAP1B')
plot(c(1:10),main="Dysplastic Neuron,DN")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('VIM','NES','BCL2','PROM1','MAP2')
plot(c(1:10),main="Balloon Cell,BC")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('GAFP','AQP4','SOX9','ALDH1L1','ATP13A4','FGFR3')
plot(c(1:10),main="Astrocyte")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('AIF1','C1QB','CTSS','CCL3','CCL4')
plot(c(1:10),main="Microglia")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('CSPG4','PDGFRA','VCAN')
plot(c(1:10),main="OPC")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('MOG','MAG','OLIG2','OLIG1','PLP1')
plot(c(1:10),main="Oligodendrocyte")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('VWF','A2M','APOLD1','FLT1','CLDN5','SLC2A1')
plot(c(1:10),main="Endothelial")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)
marker=c('ACTA2','CSPG4')
plot(c(1:10),main="Pericyte")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

dev.off()

#########################################
library(scRNAtoolVis)
clusterCornerAxes(object = merge,reduction = 'umap',
noSplit = T,
clusterCol = 'CellType',
cornerTextSize = 3.5,
themebg = 'bwCorner',
addCircle = TRUE,
cicAlpha = 0.2,
nbin = 200)







for (i in 1:length(ABC@meta.data$orig.ident)){
  if (ABC@meta.data$manucelltype[i] == '46') {
    print(ABC@meta.data$CellType[i])
    print(ABC@meta.data$manucelltype[i])
    ABC@meta.data$manucelltype[i] = ABC@meta.data$CellType[i]
    print(ABC@meta.data$manucelltype[i])
      }
}



####  0.library using codepack  ####
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(CellChat)
library(ggalluvial)
library(svglite)
library(monocle3)
library(tidyverse)
library(NMF)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(GSVA)
library(pheatmap)
library(fgsea)
library(tibble)
library(SeuratWrappers)
library(gridExtra)
library(Matrix)
library(ggsci)
library(GOplot)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
setwd("C:/Users/11723/Desktop/HS")
merge = readRDS('HS_merge.rds')
#项目名称，如H02001。。。。
project.name = '220323_HS2'
#样本疾病类型，如FCD、HS。。。。
disease.type = 'HS'
#输出文件地址（工作地址）
output.address = r"(C:\Users\11723\Desktop\project)"
#分类参考（该列名称，如：CellType\CellGruop\manucelltype）
divide_by = 'manucelltype'
#
ill.data.Read.10X = gsub("\\\\","/",ill.data.Read.10X)
output.address = gsub("\\\\","/",output.address)
####  1.load data  ####
health = readRDS('C:/Users/11723/Desktop/project/health.rds')
merge@meta.data$Group = 'HS'
bp_HS = SplitObject(merge,split.by = 'orig.ident')
head(bp_HS[[1]])









FeaturePlot(merge, features =c('PTPRC','SALL1',"SPARC"),min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))

FeaturePlot(merge, features =c('CST7','GM1673','ST8SIA6','ICAM1','LPL','CLEC7A','PSAT1','CD63'),min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))

FeaturePlot(merge, features =c('APOE','MS4A7','MS4A6C','LYZ2','TGFBI'),min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))

FeaturePlot(merge, features =c('P2RX7','MRC1','FOLR2','NRP1','CD63','CIITA','CLEC12A','PTPRC'),min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))

FeaturePlot(merge, features =c('APOE','MS4A7','LYZ2','CD63','MRC1','CD163','GAS6','NRP1','IGF1','CLEC10A','CD209F','FOLR2','CLEC12A','CD74','H2-AA','KLRA2','LYVE1','P2RX7','CCL8','PLA2G2D','CCR2','LILRA5','TTR','CST7'),min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))




#神经元按皮层分类
FeaturePlot(ABC, features =c('LAMP5','PRSS12','CUX2','RORB','GRIN3A','PCP4','HTR2C','TLE4','GRIK3','OPRK1','NR4A2'),min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))

#兴奋性神经元分层，四类，每类3个marker
FeaturePlot(merge, features =c('CUX2','PDGFD','LAMP5','RORB','COBLL1','GABRG1','THEMIS','NTNG2','NR4A2','FEZF2','SEMA3E','IL26'),min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"),ncol=3)



######################################################################

ABC.data <- Read10X(data.dir = "G:/raw data/result-SHJS-2022-02871/result-SHJS-2022-02871/normal/220323_HS1/2.Basic_analysis/2.2.filtered_feature_bc_matrix")
HS2.1 <- CreateSeuratObject(counts = ABC.data, project = "HS2.1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = "G:/raw data/result-SHJS-2022-02871/result-SHJS-2022-02871/normal/220323_HS2/2.Basic_analysis/2.2.filtered_feature_bc_matrix")
HS2.2 <- CreateSeuratObject(counts = ABC.data, project = "HS2.2", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = "G:/raw data/SHJS-2022-02606/result-SHJS-2022-02606/normal/HS1_220218/2.Basic_analysis/2.2.filtered_feature_bc_matrix")
HS1.1 <- CreateSeuratObject(counts = ABC.data, project = "HS1.1", min.cells = 3, min.features = 200)

ABC.data <- Read10X(data.dir = "G:/raw data/SHJS-2022-02606/result-SHJS-2022-02606/normal/HS2_220218/2.Basic_analysis/2.2.filtered_feature_bc_matrix")
HS1.2 <- CreateSeuratObject(counts = ABC.data, project = "HS1.2", min.cells = 3, min.features = 200)


brain = list()
brain[[1]]=FCD
brain[[2]]=noFCD
brain[[3]]=HS1
brain[[4]]=HS2


for (i in 1:length(brain)) {
  brain[[i]][["percent.mt"]] <- PercentageFeatureSet(brain[[i]], pattern = "^MT-")
  brain[[i]] <- subset(brain[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
  brain[[i]] <- SCTransform(brain[[i]], assay = "RNA",variable.features.n = 2000, return.only.var.genes = FALSE, verbose = FALSE)
  brain[[i]] <- RunPCA(brain[[i]], assay = "SCT", verbose = FALSE)
  brain[[i]] <- FindNeighbors(brain[[i]], reduction = "pca", dims = 1:50)
  brain[[i]] <- FindClusters(brain[[i]], verbose = FALSE, resolution = 4)
  brain[[i]] <- RunUMAP(brain[[i]], reduction = "pca", dims = 1:50)
  brain[[i]] <- RunTSNE(object=brain[[i]], reduction = "pca" , dims=1:50)
  }


FCD = brain[[1]]
noFCD = brain[[2]]
HS1 = brain[[3]]
HS2 = brain[[4]]




setwd("C:/Users/11723/Desktop/sn")










FeaturePlot(merge, features =c("RBFOX3","GFAP","TLR4","NEFH","NEFL","NEFM","MAP1B","VIM","NES"),min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))













Idents(merge) = merge@meta.data$orig.ident
merge = RenameIdents(merge,'H02001'='H','H0000101'='H','L02001'='L','L0000102'='L')
brain = SplitObject(merge,split.by = 'manucelltype')

setwd("C:/Users/11723/Desktop/test/manucelltype_diff")

for (i in 1:length(brain)) {
  ABC = brain[[i]]
  dir.create(as.character(ABC@meta.data$manucelltype[[1]]))
  setwd(paste('./',as.character(ABC@meta.data$manucelltype[1]),sep=''))
  #ABC = brain[[i]]
  diff1 = FindMarkers(brain[[i]],ident.1 = 'H',ident.2 = 'L',min.pct = 0.25,only.pos = FALSE)
  #diff2 = FindMarkers(brain[[i]],ident.1 = 'H0000101',ident.2 = 'L0000102',min.pct = 0.25,only.pos = FALSE)
  write.csv(diff1,'H&L.csv')
  #write.csv(diff2,"H0000101&L0000102.csv")
  data = diff1
  ggsave('Vocanol.png',EnhancedVolcano(data, lab = rownames(data), x = 'avg_log2FC' , y = 'p_val_adj' , xlim = c(-10,10) , ylim = c(0,125) , pCutoff = 0.001 , FCcutoff = 0.5),width = 10,height = 10)
  #data = diff2
  #ggsave('Vocanol2.png',EnhancedVolcano(data, lab = rownames(data), x = 'avg_log2FC' , y = 'p_val_adj' , xlim = c(-10,10) , ylim = c(0,125) , pCutoff = 0.001 , FCcutoff = 0.5),width = 10,height = 10)
  #dir.create('H02001&L02001')
  #setwd('./H02001&L02001')
  markers = diff1
  markers$gene = rownames(markers)
  up=rownames( markers [intersect(which( markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=0.5)),])
  down=rownames( markers [intersect(which( markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]<=-0.5)),])
  dir.create('up')
  setwd('./up')
  x = markers[intersect(up,p2ry12),]
  eg=bitr(x$gene,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb="org.Hs.eg.db")
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #go.BP <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #go.CC <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #go.MF <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #Go.BP <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #Go.CC <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #Go.MF <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  ggsave("1.1.1.GO_barplot.pdf",barplot(go,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30),width=14,height=10)
  ggsave("1.1.1.GO_barplot.png",barplot(go,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30),width=14,height=10)
  setwd('../')
  dir.create('down')
  setwd('./down')
  x = markers[intersect(down,p2ry12),]
  eg=bitr(x$gene,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb="org.Hs.eg.db")
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #go.BP <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #go.CC <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #go.MF <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #Go.BP <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #Go.CC <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  #Go.MF <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  ggsave("1.1.1.GO_barplot.pdf",barplot(go,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30),width=14,height=10)
  ggsave("1.1.1.GO_barplot.png",barplot(go,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  #ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30),width=14,height=10)
  setwd('../')
  setwd('../')
  }



  dir.create('H0000101&L0000102')
  setwd('./H0000101&L0000102')
  markers = diff2
  markers$gene = rownames(markers)
  up=rownames( markers [intersect(which( markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=0.5)),])
  down=rownames( markers [intersect(which( markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]<=-0.5)),])
  dir.create('up')
  setwd('./up')
  x = markers[up,]
  eg=bitr(x$gene,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb="org.Hs.eg.db")
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.BP <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.CC <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.MF <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.BP <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.CC <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.MF <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  ggsave("1.1.1.GO_barplot.pdf",barplot(go,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30),width=14,height=10)
  ggsave("1.1.1.GO_barplot.png",barplot(go,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30),width=14,height=10)
  setwd('../')
  dir.create('down')
  setwd('./down')
  x=markers[down,]
  eg=bitr(x$gene,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb="org.Hs.eg.db")
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.BP <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.CC <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.MF <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.BP <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.CC <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.MF <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  ggsave("1.1.1.GO_barplot.pdf",barplot(go,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30),width=14,height=10)
  ggsave("1.1.1.GO_barplot.png",barplot(go,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30),width=14,height=10)
  setwd('../')
  setwd('../')
  setwd('../')
}



############################################################################


test = ABC.markers






a=c('celltype','FC>5','FC>2','FC>0.5','FC<-0.5','FC<-2','FC<-5')
frame = data.frame(a)

for (i in 1:length(brain)) {
ABC = brain[[i]]
#ABC = brain[[i]]
diff1 = FindMarkers(brain[[i]],ident.1 = 'H',ident.2 = 'L',min.pct = 0.25,only.pos = FALSE)
#diff2 = FindMarkers(brain[[i]],ident.1 = 'H0000101',ident.2 = 'L0000102',min.pct = 0.25,only.pos = FALSE)
#write.csv(diff2,"H0000101&L0000102.csv")
data = diff1
#data = diff2
#ggsave('Vocanol2.png',EnhancedVolcano(data, lab = rownames(data), x = 'avg_log2FC' , y = 'p_val_adj' , xlim = c(-10,10) , ylim = c(0,125) , pCutoff = 0.001 , FCcutoff = 0.5),width = 10,height = 10)
#dir.create('H02001&L02001')
#setwd('./H02001&L02001')
markers = data
up5=rownames( markers [intersect(which( markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=5)),])
up2=rownames( markers [intersect(which( markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=2)),])
up0.5=rownames( markers [intersect(which( markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=0.5)),])
down0.5=rownames( markers [intersect(which( markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]<=-0.5)),])
down2=rownames( markers [intersect(which( markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]<=-2)),])
down5=rownames( markers [intersect(which( markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]<=-5)),])
u5=length(up5)
u2=length(up2)
u0.5=length(up0.5)
d0.5=length(down0.5)
d2=length(down2)
d5=length(down5)
list = c(as.character(ABC@meta.data$manucelltype[1]),u5,u2,u0.5,d0.5,d2,d5)
frame[i+1]=list

}









rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)


sce.markers <- read.csv('amarkers.csv',row.names = 1)
table(sce.markers$cluster)
#dPVL  iCAFs  imPVL myCAFs 
#   120     22    129     46

H02001=sce.markers[sce.markers$cluster=='H02001',]$gene
L02001=sce.markers[sce.markers$cluster=='L02001',]$gene
H0000101=sce.markers[sce.markers$cluster=='H0000101',]$gene
L0000102=sce.markers[sce.markers$cluster=='L0000102',]$gene
total <- list(H02001=H02001,L02001=L02001,H0000101=H0000101,L0000102=L0000102)

# 转ID
library(org.Hs.eg.db)
i=1
for (i in 1:4) {
  #i=1
  ## 把SYMBOL改为ENTREZID
  total[[i]]=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,keys = total[[i]],columns = 'ENTREZID',keytype = 'SYMBOL')[,2]))
  }
lapply(total, head)

xx <- compareCluster(total, fun="enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=0.99)
table(xx@compareClusterResult$Cluster)
head(as.data.frame(xx))

dotplot(xx) #气泡图
df_go_diff <- as.data.frame(xx)












setwd("C:/Users/11723/Desktop/test")
merge = readRDS("./merge.rds")
Idents(merge) = merge@meta.data$orig.ident
ABC = RenameIdents(merge,'H02001'='H','H0000101'='H','L02001'='L','L0000102'='L')
ABC@meta.data$HL = Idents(ABC)
table(ABC@meta.data$HL,ABC@meta.data$manucelltype)
brain = SplitObject(ABC,split.by = 'HL')
brain




for (i in length(ABC@meta.data$manusubcelltype2)) {
  if (as.character(ABC@meta.data$manusubcelltype2[i]) == '1'){
    ABC@meta.data$manusubcelltype2[i] = 'Microglia_JUN+'
  }else if (as.character(ABC@meta.data$manusubcelltype2[i]) == '2'){
    ABC@meta.data$manusubcelltype2[i] = 'Microglia_APOE+'
  }else if (as.character(ABC@meta.data$manusubcelltype2[i]) == '3'){
    ABC@meta.data$manusubcelltype2[i] = 'Microglia_CCL4+'
  }else if (as.character(ABC@meta.data$manusubcelltype2[i]) == '4'){
    ABC@meta.data$manusubcelltype2[i] = 'Microglia_P2RY12+'
  }else if (as.character(ABC@meta.data$manusubcelltype2[i]) == '5'){
    ABC@meta.data$manusubcelltype2[i] = 'Microglia_JUN+'
  }else if (as.character(ABC@meta.data$manusubcelltype2[i]) == '6'){
    ABC@meta.data$manusubcelltype2[i] = 'Microglia_APOC+'
  }else {
    ABC@meta.data$manusubcelltype2[i] = ABC@meta.data$manucelltype[i]
  }
}

##############################################################
ABC = HS2
Idents(ABC) = ABC@meta.data$CellType
ABC = RenameIdents(ABC,'Ex.L2_L3' = 'Ex','Ex.L3_L5' = 'Ex','Ex.L5' = 'Ex','Ex.L5_L6' = 'Ex','Ex.L5b_UMN_CT' = 'Ex','Ex.L5b_UMN_PT' = 'Ex','Ex.L6b' = 'Ex','In.5HT3aR_VIPneg' = 'In','In.5HT3aR_VIPpos' = 'In','In.Basket_PVALB' = 'In','In.Chandelier_PVALB' = 'In','In.Rosehip' = 'In','In.SST_NPYneg' = 'In','In.SST_NPYpos' = 'In')
ABC@meta.data$CellType = Idents(ABC)
HS2 = ABC


###############################################################


for (i in 1:length(brain)){
  dir.create(names(brain[i]))
  setwd(paste('./',names(brain[i]),sep = ''))
  markers = FindAllMarkers(brain[[i]],min.pct = 0.25,only.pos = TRUE)
  up=rownames( markers[intersect(which( markers[,"p_val"]<0.05),which( markers[,"cluster"]=='PN')),])
  PN = markers[up,]
  up=rownames( markers[intersect(which( markers[,"p_val"]<0.05),which( markers[,"cluster"]=='FCD')),])
  FCD = markers[up,]
  up=rownames( markers[intersect(which( markers[,"p_val"]<0.05),which( markers[,"cluster"]=='HS')),])
  HS = markers[up,]
  
  dir.create('PN')
  setwd('./PN')
  x = PN
  eg=bitr(x$gene,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb="org.Hs.eg.db")
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.BP <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.CC <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.MF <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.BP <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.CC <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.MF <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  ggsave("1.1.1.GO_barplot.pdf",barplot(go,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30),width=14,height=10)
  ggsave("1.1.1.GO_barplot.png",barplot(go,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30),width=14,height=10)
  setwd('../')
  
  dir.create('FCD')
  setwd('./FCD')
  x = FCD
  eg=bitr(x$gene,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb="org.Hs.eg.db")
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.BP <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.CC <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.MF <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.BP <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.CC <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.MF <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  ggsave("1.1.1.GO_barplot.pdf",barplot(go,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30),width=14,height=10)
  ggsave("1.1.1.GO_barplot.png",barplot(go,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30),width=14,height=10)
  setwd('../')
  
  dir.create('HS')
  setwd('./HS')
  x = HS
  eg=bitr(x$gene,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb="org.Hs.eg.db")
  genelist <- eg$ENTREZID
  go <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.BP <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.CC <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  go.MF <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.BP <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.CC <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  Go.MF <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
  ggsave("1.1.1.GO_barplot.pdf",barplot(go,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.2.1.GO_BP_barplot.pdf",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.3.1.GO_CC_barplot.pdf",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.4.1.GO_MF_barplot.pdf",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  ggsave("2.1.KEGG_dotplot.pdf",dotplot(kegg, showCategory=30),width=14,height=10)
  ggsave("1.1.1.GO_barplot.png",barplot(go,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.2.1.GO_BP_barplot.png",barplot(go.BP,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.3.1.GO_CC_barplot.png",barplot(go.CC,showCategory=20,drop=T),width=14,height=10)
  ggsave("1.4.1.GO_MF_barplot.png",barplot(go.MF,showCategory=20,drop=T),width=14,height=10)
  ggsave("2.1.KEGG_dotplot.png",dotplot(kegg, showCategory=30),width=14,height=10)
  setwd('../')
  setwd('../')
}






marker = c('CD74',
           'PTPRC',
           'CD46',
           'RIPK2',
           'SYK',
           'HLA-DRB1',
           'LGALS9',
           'CD81',
           'PYCARD',
           'HLA-DMB',
           'RPS3',
           'HLA-DPA1',
           'HLA-DPB1',
           'HMGB1',
           'TMIGD2',
           'SELENOK',
           'B2M',
           'LILRB4',
           'NCKAP1L',
           'AIF1',
           'VSIR',
           'CORO1A',
           'TGFBR2',
           'TNFRSF13C',
           'NCK2',
           'HLA-DRA',
           'XBP1',
           'IGF1',
           'NLRP3',
           'CD4',
           'NFKBID',
           'HLA-DRB5',
           'HLA-DMA',
           'PELI1',
           'HLA-DQB1',
           'DOCK8',
           'RHOA',
           'SIRPA',
           'CSK',
           'SOX4'
)

