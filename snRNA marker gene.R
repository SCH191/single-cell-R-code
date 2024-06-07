setwd("G:/work/HS")
merge = readRDS('HS_merge.rds')

setwd("G:/work/HS8")
merge = readRDS('merge8.rds')

setwd("G:/work/sn/3.HS")
merge = readRDS("merge.rds")

setwd("G:/work/mouse")
merge = readRDS('merge.rds')


setwd("G:/work/AD online")

setwd("G:/work/sn/1.FCD")

FCD=readRDS("FCD.rds")
noFCD=readRDS("noFCD.rds")
merge = readRDS("merge.rds")
HS1=readRDS("HS1.rds")
HS2=readRDS("HS2.rds")


brain = list()
brain[[1]]=FCD
brain[[2]]=noFCD
brain[[3]]=HS1
brain[[4]]=HS2




i=1

dir.create(paste("number",i,"sample"))
setwd(paste("number",i,"sample"))

j=brain[[i]]
name=i
reduciton = 'umap' #'umap','tsne'

ggsave(paste(i,"_umap.png"),DimPlot(j, reduction = "umap", group.by = c("seurat_clusters", "orig.ident",'manucelltype'),label = T,repel = T)/DimPlot(j, reduction = "tsne", group.by = c("seurat_clusters", "orig.ident"),label = T,repel = T),width=14,height=10)

#pdf marker gene Hs##############################################
pdf(paste("out_",name,"_marker_gene.pdf",sep=""),width=18,height=12)

marker=c("RBFOX3","GAD1","SLC17A7","SLC1A3")
plot(c(1:10),main="In & Ex $ non_Neuro")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, reduction = reduciton, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('SLC17A7','SATB2','RORB','CUX2','TLE4','NR4A2','SEMA3C')
plot(c(1:10),main="Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, reduction = reduciton,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('LAMP5','CUX2','PDGFD','FREM3','PRSS12','COL5A2')
plot(c(1:10),main="L2_3_Cux2 Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, reduction = reduciton,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('RORB','COBLL1','SCHLAP1','MME','PLCH1')
plot(c(1:10),main="L4_Rorb Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, reduction = reduciton,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('FEZF2','TLE4','LRRK1','ADRA1A','RORB')
plot(c(1:10),main="L5_6_Fezf2 Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j,reduction = reduciton, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('THEMIS','NR4A2','PXDN','OPRK1')
plot(c(1:10),main="L5_6_Themis Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j,reduction = reduciton, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('GAD1','GAD2','SOX6','PVALB','SST','VIP','LHX6','CALB2','SULF1')
plot(c(1:10),main="Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, reduction = reduciton,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('ID2','LAMP5','CALB2','PAX6','RELN','SV2C','COL5A2','FAM19A1','SPOCK1','MYO5B')
plot(c(1:10),main="Id2 Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, reduction = reduciton,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('PVALB','NOS1','SULF1','LHX6','KCNS3','CRH','PLEKHH2','LGR5')
plot(c(1:10),main="Pvalb Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, reduction = reduciton,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('SST','NOS1','SEMA6A','FAM89A','LHX6','GRIK1')
plot(c(1:10),main="Sst Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j,reduction = reduciton, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('VIP','TAC3','CALB2','LAMA3','FAM19A1','NPR3')
plot(c(1:10),main="Vip Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, reduction = reduciton,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('AQP4','FGFR3',"SLC4A4","GFAP",'SLC1A2','GJA1','SLC1A3','S100B','SOX9')
plot(c(1:10),main="Astrocytes")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, reduction = reduciton,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('CTSS','C1QB',"C1QC","CD68","TREM2","CCL4",'ITGAM','CMTM6','AIF1','ADRB2','CSF1R','P2RY12','CD40','TMEM119')
plot(c(1:10),main="Microglia")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j,reduction = reduciton, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('CSPG4','PDGFRA','VCAN',"GPR17","CCND1","TNR","NEU4",'NKX6-2','FYN','TNR','ALDOC','PCDH15','OLIG1')
plot(c(1:10),main="Oligodendrocyte.Precursors")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j,reduction = reduciton, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('MOG','MAG',"OPALIN","KLK6","GPR37","ANLN",'MBP','NIPAL4','CLDN14','PEX5L','SEC14L5')
plot(c(1:10),main="Oligodendrocytes")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j,reduction = reduciton, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c("CLDN5","FLT1","PECAM1","VWF","CLEC14A")
plot(c(1:10),main="Endothelial")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, reduction = reduciton,features =marker,min.cutoff="q10",max.cutoff="q95")
DotPlot(j, features =marker)

marker=c('DCN','PTGDS','ATP1A2','ITIH5','FLT1','RABL2','FOXJ1','PIFO','RSPH1','TEKT4','CALML4','CRYGN','HDC','FAM183B')
plot(c(1:10),main="Vascular")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j,reduction = reduciton, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)



dev.off()


#png marker gene Hs#############################
dir.create(paste("out_",name,"_marker_gene_dir",sep=""))
setwd(paste("./out_",name,"_marker_gene_dir",sep=""))


#ggsave(paste(".1_",main,"_VlnPlot",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0.1)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
#ggsave(paste(".2_",main,"_FeaturePlot",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
#ggsave(paste(".3_",main,"_DotPlot",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
#ggsave(paste(".4_",main,"_RidgePlot",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)


marker=c("RBFOX3","GAD1","SLC17A7","SLC1A3")
main="In & Ex $ non_Neuro"
ggsave(paste("1.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("1.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("1.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("1.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 2),width = 10,height = 10)

marker=c('SLC17A7','SATB2','RORB','CUX2','TLE4','NR4A2','SEMA3C')
main="Excitatory"
ggsave(paste("2.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("2.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, reduction = reduciton,raster=FALSE,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("2.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("2.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('LAMP5','CUX2','PDGFD','FREM3','PRSS12','COL5A2')
main="L2_3_Cux2 Excitatory"
ggsave(paste("2.11_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("2.21_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, reduction = reduciton,raster=FALSE,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("2.31_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("2.41_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('RORB','COBLL1','SCHLAP1','MME','PLCH1')
main="L4_Rorb Excitatory"
ggsave(paste("2.112_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("2.212_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, reduction = reduciton,raster=FALSE,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("2.312_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("2.412_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('FEZF2','TLE4','LRRK1','ADRA1A','RORB')
main="L5_6_Fezf2 Excitatory"
ggsave(paste("2.1113_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("2.2113_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, reduction = reduciton,raster=FALSE,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("2.3113_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("2.4113_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('THEMIS','NR4A2','PXDN','OPRK1')
main="L5_6_Themis Excitatory"
ggsave(paste("2.11114_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("2.21114_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, reduction = reduciton,raster=FALSE,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("2.31114_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("2.41114_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)


marker=c('GAD1','GAD2','SOX6','PVALB','SST','VIP','LHX6','CALB2','SULF1')
main="Inhibitory"
ggsave(paste("3.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 20)
ggsave(paste("3.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("3.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("3.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('ID2','LAMP5','CALB2','PAX6','RELN','SV2C','COL5A2','FAM19A1','SPOCK1','MYO5B')
main="Id2 Inhibitory"
ggsave(paste("4.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 20)
ggsave(paste("4.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("4.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("4.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('PVALB','NOS1','SULF1','LHX6','KCNS3','CRH','PLEKHH2','LGR5')
main="Pvalb Inhibitory"
ggsave(paste("5.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 18)
ggsave(paste("5.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("5.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("5.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('SST','NOS1','SEMA6A','FAM89A','LHX6','GRIK1')
main="Sst Inhibitory"
ggsave(paste("6.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("6.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("6.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("6.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('VIP','TAC3','CALB2','LAMA3','FAM19A1','NPR3')
main="Vip Inhibitory"
ggsave(paste("7.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("7.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("7.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("7.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('AQP4','FGFR3',"SLC4A4","GFAP",'SLC1A2','GJA1','SLC1A3','S100B','SOX9')
main="Astrocytes"
ggsave(paste("8.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 20)
ggsave(paste("8.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("8.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("8.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('CTSS','C1QB',"C1QC","CD68","TREM2","CCL4",'ITGAM','CMTM6','AIF1','ADRB2','CSF1R','P2RY12','CD40','TMEM119')
main="Microglia"
ggsave(paste("9.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 20)
ggsave(paste("9.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("9.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("9.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('CSPG4','PDGFRA','VCAN',"GPR17","CCND1","TNR","NEU4",'NKX6-2','FYN','ALDOC','PCDH15','OLIG1')
main="Oligodendrocyte.Precursors"
ggsave(paste("10.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 20)
ggsave(paste("10.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("10.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("10.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('MOG','MAG',"OPALIN","KLK6","GPR37","ANLN",'MBP','NIPAL4','CLDN14','PEX5L','SEC14L5')
main="Oligodendrocytes"
ggsave(paste("11.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 20)
ggsave(paste("11.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("11.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("11.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("PECAM1","VWF","CLDN5","FLT1","CLEC14A")
main="Endothelial"
ggsave(paste("12.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("12.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("12.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("12.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("EPCAM","CD24",'KRT19','PROM1','ALDH1A1')
main="Epithelial_cell"
ggsave(paste("13.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("13.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("13.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("13.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("KCNMB1","MRVI1","KCNAB1","DMPK")
main="Smooth_muscle_cells"
ggsave(paste("14.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("14.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("14.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("14.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("COL6A2","VIM","COL1A1","PDGFRB")
main="Fibroblasts"
ggsave(paste("15.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("15.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("15.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("15.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('DCN','PTGDS','ATP1A2','ITIH5','FLT1','RABL2','FOXJ1','PIFO','RSPH1','TEKT4','CALML4','CRYGN','HDC','FAM183B')
main="Vascular"
ggsave(paste("16.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 25)
ggsave(paste("16.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("16.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("16.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('PTPRC')
main="immu"
ggsave(paste("17.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 6)
ggsave(paste("17.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("17.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("17.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('FGFBP2','FCG3RA','CXCR1',"KLRD1","NKG7","GZMA","XCL1","GZMK")
main="NK"
ggsave(paste("18.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 18)
ggsave(paste("18.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton, raster=FALSE,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("18.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("18.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('CD19','MS4A1',"CD79A","CCR7","CD27","IGLC2","CXCR4")
main="B_cell"
ggsave(paste("19.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 18)
ggsave(paste("19.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton, raster=FALSE,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("19.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("19.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('IGHG1','MZB1','SDC1','CD79A')
main="Plasma cells"
ggsave(paste("20.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("20.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton, raster=FALSE,features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("20.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("20.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("CD3D","CD3E","CD8A","TRBC2","CD3G","GZMK")
main="T_cell"
ggsave(paste("21.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("21.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("21.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("21.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("APOBEC3A","TET2","MEFV","HCK","ITGAX","RBPJ")
main="Monocytes"
ggsave(paste("22.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 15)
ggsave(paste("22.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("22.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("22.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("CD163",'CD68','CD14')
main="Macrophages"
ggsave(paste("23.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("23.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("23.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("23.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("CEACAM8","ITGAM","FUT4","FCGR3A","CXCL8","FCGR3B","MNDA")
main="Neutrophlis"
ggsave(paste("24.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 18)
ggsave(paste("24.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("24.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("24.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('ITGAX','CX3CR1','CD86','CD83')
main="DC"
ggsave(paste("25.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("25.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j,reduction = reduciton,raster=FALSE, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("25.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("25.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)


setwd('../')

#png marker gene Mm######################################
dir.create(paste("out_",name,"_marker_gene_dir",sep=""))
setwd(paste("./out_",name,"_marker_gene_dir",sep=""))


#ggsave(paste(".1_",main,"_VlnPlot",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0.1)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
#ggsave(paste(".2_",main,"_FeaturePlot",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
#ggsave(paste(".3_",main,"_DotPlot",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
#ggsave(paste(".4_",main,"_RidgePlot",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)


marker=c("Rbfox3","Gad1","Slc17a7","Slc1a3")
main="In & Ex $ non_Neuro"
ggsave(paste("1.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("1.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("1.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("1.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Slc17a7','Satb2','Rorb','Cux2','Tle4','Nr4a2','Sema3c')
main="Excitatory"
ggsave(paste("2.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("2.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("2.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("2.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Gad1','Gad2','Sox6','Pvalb','Sst','Vip','Lhx6','Calb2','Sulf1')
main="Inhibitory"
ggsave(paste("3.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("3.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("3.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("3.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Id2','Lamp5','Calb2','Pax6','Reln','Sv2c','Col5a2','Fam19a1','Spock1','Myo5b')
main="Id2 Inhibitory"
ggsave(paste("4.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("4.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("4.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("4.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Pvalb','Nos1','Sulf1','Lhx6','Kcns3','Crh','Plekhh2','Lgr5')
main="Pvalb Inhibitory"
ggsave(paste("5.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("5.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("5.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("5.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Sst','Nos1','Sema6a','Fam89a','Lhx6','Grik1')
main="Sst Inhibitory"
ggsave(paste("6.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("6.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("6.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("6.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Vip','Tac3','Calb2','Lama3','Fam19a1','Npr3')
main="Vip Inhibitory"
ggsave(paste("7.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("7.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("7.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("7.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Aqp4','Fgfr3',"Slc4a4","Gfap",'Slc1a2','Gja1','Slc1a3','S100b','Sox9')
main="Astrocytes"
ggsave(paste("8.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("8.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("8.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("8.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('C1qb',"C1qc","Cd68","Trem2","Ccl4",'Itgam','Csf1r','P2ry12','Tmem119')
main="Microglia"
ggsave(paste("9.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("9.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("9.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("9.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Cspg4','Pdgfra',"Ccnd1","Tnr","Neu4",'Aldoc','Pcdh15','Olig1')
main="Oligodendrocyte.Precursors"
ggsave(paste("10.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("10.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("10.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("10.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Mog','Mag',"Opalin","Klk6","Gpr37",'Mbp','Cldn14')
main="Oligodendrocytes"
ggsave(paste("11.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("11.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("11.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("11.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Pecam1","Vwf","Cldn5","Flt1","Clec14a")
main="Endothelial"
ggsave(paste("12.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("12.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("12.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("12.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Epcam","Cd24",'Krt19','Prom1','Aldh1a1')
main="Epithelial_cell"
ggsave(paste("13.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("13.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("13.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("13.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Kcnmb1","Mrvi1","Kcnab1","Dmpk")
main="Smooth_muscle_cells"
ggsave(paste("14.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("14.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("14.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("14.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Col6a2","Vim","Col1a1","Pdgfrb")
main="Fibroblasts"
ggsave(paste("15.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("15.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("15.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("15.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Dcn','Ptgds','Atp1a2','Itih5','Flt1','Rabl2','Tekt4')
main="Vascular"
ggsave(paste("16.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("16.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("16.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("16.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Cd45','Ptprc')
main="immu"
ggsave(paste("17.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("17.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("17.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("17.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)


marker=c('Fgfbp2','Fcg3ra','Cxcr1',"Klrd1","Nkg7","Gzma","Xcl1","Gzmk")
main="NK"
ggsave(paste("18.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("18.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("18.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("18.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Cd19','Ms4a1',"Cd79a","Ccr7","Cd27","Iglc2","Cxcr4")
main="B_cell"
ggsave(paste("19.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("19.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("19.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("19.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Ighg1','Mzb1','Sdc1','Cd79a')
main="Plasma cells"
ggsave(paste("20.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("20.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("20.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("20.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Cd3d","Cd3e","Cd8a","Trbc2","Cd3g","Gzmk")
main="T_cell"
ggsave(paste("21.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("21.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("21.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("21.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Apobec3A","Tet2","Mefv","Hck","Itgax","Rbpj")
main="Monocytes"
ggsave(paste("22.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("22.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("22.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("22.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Cd163",'Cd68','Cd14')
main="Macrophages"
ggsave(paste("23.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("23.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("23.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("23.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Ceacam8","Itgam","Fut4","Fcgr3a","Cxcl8","Fcgr3b","Mnda")
main="Neutrophlis"
ggsave(paste("24.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("24.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("24.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("24.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Itgax','Cx3cr1','Cd86','Cd83')
main="DC"
ggsave(paste("25.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("25.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("25.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("25.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)


setwd('../')

#png marker gene Mm online part#############################################
dir.create(paste("out_",name,"_marker_gene_Mm_online_dir",sep=""))
setwd(paste("./out_",name,"_marker_gene_Mm_online_dir",sep=""))


#ggsave(paste(".1_",main,"_VlnPlot",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0.1)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
#ggsave(paste(".2_",main,"_FeaturePlot",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
#ggsave(paste(".3_",main,"_DotPlot",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
#ggsave(paste(".4_",main,"_RidgePlot",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)


marker=c("Cx3cr1", "P2ry12", "Tmem119")
main="Microglia and non-parenchymal macrophages"
ggsave(paste("1.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("1.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("1.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("1.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Cd163")
main="onparenchymal macrophage"
ggsave(paste("2.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("2.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("2.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("2.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Pdgfra", "Susd5","Cspg4")
main="oligodendrocyte lineage cells"
ggsave(paste("3.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("3.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("3.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("3.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Mog", "Mag", "Gjc2")
main="mature myelinating oligodendrocytes"
ggsave(paste("4.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("4.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("4.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("4.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Pecam1", "Cldn5")
main="pan-endothelial cells"
ggsave(paste("5.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("5.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("5.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("5.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Slco1c1", "Ocln")
main="blood-brain barrier"
ggsave(paste("6.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("6.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("6.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("6.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Aldh1l1", "Slc1a3", "Aqp4")
main="astrocyte"
ggsave(paste("7.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("7.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("7.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("7.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Gdf10", "Vim", "Nbl1", "A2m")
main="Bergman glia"
ggsave(paste("8.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("8.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("8.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("8.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Slc17a7", "Neurod6", "Mab21l1")
main="excitatory neurons"
ggsave(paste("9.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("9.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("9.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("9.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Gad1', "Reln", "Calb1")
main="inhibitory neurons"
ggsave(paste("10.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("10.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("10.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("10.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Dcx", "Dlx2", "Ascl1", "Hes5")
main="neuroprogenitor cells"
ggsave(paste("11.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("11.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("11.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("11.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Des", "Mcam", "Pdgfrb")
main="pericytes"
ggsave(paste("12.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("12.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("12.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("12.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c('Bmx', "Jag1", "Efnb2")
main="arterioles"
ggsave(paste("13.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("13.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("13.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("13.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Nr2f2", 'Flt4', "Vwf")
main="post-capillary venules"
ggsave(paste("14.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("14.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("14.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("14.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Tfrc", "Car4", "Slc16a1")
main="capillaries"
ggsave(paste("15.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("15.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("15.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("15.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Jag1", "Notch1", "Hey1", "Vegfc", "Edn1", "Tmem100")
main="pro-angiogenic arterioles"
ggsave(paste("16.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("16.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("16.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("16.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c(" Vcam1", "Icam1", "Lcn2","Hif1a", "Vwf", "Csf1")
main="post-capillary venules and capillaries enriched with inflammatory genes"
ggsave(paste("17.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("17.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("17.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("17.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

marker=c("Ptprc", "H2-Eb1","Cd14")
main="immune cells"
ggsave(paste("18.1_",main,"_VlnPlot.png",sep=""),VlnPlot(j, features =marker,ncol=1,pt.size = 0)+geom_boxplot(width = 0.2,col = 'black', fill='white'),width = 10,height = 10)
ggsave(paste("18.2_",main,"_FeaturePlot.png",sep=""),FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred")),width = 10,height = 10)
ggsave(paste("18.3_",main,"_DotPlot.png",sep=""),DotPlot(j, features =marker),width = 10,height = 10)
ggsave(paste("18.4_",main,"_RidgePlot.png",sep=""),RidgePlot(j, features = marker, ncol = 3),width = 10,height = 10)

setwd('../')

###############################################################################
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





######################################
i=1
j=brain[[i]]
name=i
pdf(paste("out",name,".pdf",sep=""),width=12,height=12)

marker=c("RBFOX3","GAD1","SLC17A7","SLC1A3")
plot(c(1:10),main="In & Ex $ non_Neuro")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('ID3','GFAP','SERPINI2','WDR49','MT1F','LOC105373158','SLC16A9','ID4','AQP4','TJP2')
plot(c(1:10),main="Astro_L1_2_FGFR3_GFAP&MT1F&ID3")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('LOC105376917','RANBP3L','FGFR3','GLI3','RNF219_AS1','SLC1A3','COL5A3')
plot(c(1:10),main="Astro_L1_6_FGFR3_SLC14A1&LOC105376917&SLC1A3")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('EPAS1','ANXA3','EMCN','ITIH5','COLEC12','EGFL7','ABCB1','FLT1','LOC105377862')
plot(c(1:10),main="Endo_L2_6_NOSTRIN&EMCN")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('CD74','CX3CR1','FYB','SYK','CSF1R','APBB1IP','LPAR6','ADAM28','DOCK8','INPP5D')
plot(c(1:10),main="Micro_L1_3_TYROBP&CSF1R&APBB1IP")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('PRRX1','PDGFRA','VCAN','LOC105371624','SEMA5A','PCDH15','PTPRZ1','LHFPL3')
plot(c(1:10),main="OPC_L1_6_PDGFRA&PCDH15&PDGFRA")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('CNDP1','MOG','LOC101927459','CARNS1','CERCAM','ENPP2','PLP1','ST18','MOBP','TMEM144')
plot(c(1:10),main="Oligo_L1_6_OPALIN&ENPP2&CERCAM")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c("ADARB2","LHX6")
plot(c(1:10),main="抑制性神经元亚型，尾神经节隆起(CGE) 和内侧神经节隆起 (MGE)")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

dev.off()




################################################################################
i=1
j=brain[[i]]
name=i


pdf(paste("out",name,".pdf",sep=""),width=12,height=12)

marker=c('SLC17A7','SATB2','RORB','CUX2','TLE4','NR4A2','SEMA3C')
plot(c(1:10),main="Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('LAMP5','CUX2','PDGFD','FREM3','PRSS12','COL5A2')
plot(c(1:10),main="L2_3_Cux2	Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('RORB','COBLL1','SCHLAP1','MME','PLCH1')
plot(c(1:10),main="L4_Rorb	Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('FEZF2','TLE4','LRRK1','ADRA1A','RORB')
plot(c(1:10),main="L5_6_Fezf2	Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('THEMIS','NR4A2','PXDN','OPRK1')
plot(c(1:10),main="L5_6_Themis	Excitatory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('GAD1','GAD2','SOX6','PVALB','SST','VIP','LHX6','NDNF','CALB2','SULF1')
plot(c(1:10),main="Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('ID2','LAMP5','CALB2','PAX6','RELN','SV2C','COL5A2','SEMA3C','FAM19A1','SPOCK1','MYO5B')
plot(c(1:10),main="Id2	Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('PVALB','NOS1','SULF1','LHX6','KCNS3','CRH','PLEKHH2','LGR5')
plot(c(1:10),main="Pvalb	Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('SST','NOS1','SEMA6A','FAM89A','LHX6','GRIK1')
plot(c(1:10),main="Sst	Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('VIP','TAC3','CALB2','LAMA3','FAM19A1','NPR3')
plot(c(1:10),main="Vip	Inhibitory")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('AQP4','GJB6','FGFR3')
plot(c(1:10),main="Astrocytes")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('CTSS','C1QB')
plot(c(1:10),main="Microglia")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('CSPG4','PDGFRA','VCAN')
plot(c(1:10),main="Oligodendrocyte.Precursors")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('MOG','MAG')
plot(c(1:10),main="Oligodendrocytes")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)

marker=c('DCN','PTGDS','ATP1A2','ITIH5','FLT1')
plot(c(1:10),main="Vascular")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95",cols = c("lightgrey", "orange","red","darkred"))
DotPlot(j, features =marker)



dev.off()


marker=c('RABL2','FOXJ1','PIFO','RSPH1','TEKT4','CALML4','CRYGN','HDC','FAM183B')
plot(c(1:10),family="YaHei Consolas",main="室管膜细胞 Ependymal cells")
VlnPlot(j, features =marker,ncol=1,pt.size = 0)
FeaturePlot(j, features =marker,min.cutoff="q10",max.cutoff="q95")


for (i in 1:length(ABC@meta.data$manusubcelltype2)) {
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
    ABC@meta.data$manusubcelltype2[i] = as.character(ABC@meta.data$manucelltype[i])
  }
}

