library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)

######


dir.create('GSEA')
setwd('./GSEA')
ABC = ABC
divide_by = 'manusubcelltype'
Idents(ABC) = ABC@meta.data$Group
disease.type = 'Focus'
control.type = 'Control'
OrgDb1="org.Hs.eg.db"   #鼠：org.Mm.eg.db 人：org.Hs.eg.db

eval(parse(text = paste('tempt_diff = SplitObject(ABC,split.by = "',divide_by,'")',sep = '')))
for (i in 1:length(tempt_diff)){
  dir.create(names(tempt_diff[i]))
  setwd(paste('./',names(tempt_diff[i]),sep = ''))
  
  Idents(tempt_diff[[i]]) = tempt_diff[[i]]@meta.data$Group
  tryCatch({
    markers <- FindAllMarkers(tempt_diff[[i]], only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  },warning = function(w){
    flag<<-FALSE
    print("warning")
  },error = function(e){
    flag<<-FALSE
    print("error")
  })
  write.csv(markers , "diff_gene.csv")
  #markers$gene = rownames(markers)
  Focus=rownames( markers [intersect(intersect(which(markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=0.5)),which( as.character(markers [,"cluster"])==disease.type)),])
  Control=rownames( markers [intersect(intersect(which(markers [,"p_val"]<0.05),which( markers [,"avg_log2FC"]>=0.5)),which( as.character(markers [,"cluster"])==control.type)),])

  
  dir.create(disease.type)
  setwd(paste('./',disease.type,sep=''))
  x = markers[Focus,c('gene','avg_log2FC')]
  tmp=gsub("\\..*","",x$gene)
  eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)
  genelist <- eg$ENTREZID
  tampt = 1
  eg$logFC=NA
  eg <- eg[!duplicated(eg$SYMBOL), ]
  for (w in 1:length(x$gene)) {
    if (x[w,'gene']==eg[tampt,'SYMBOL']){
      eg[tampt,'logFC']=as.numeric(x[w,'avg_log2FC'])
      tampt = tampt+1
    }
    if(tampt>length(eg$SYMBOL)){break}
  }
  eg <- eg[order(eg$logFC, decreasing = TRUE), ]
  geneList = eg[,'logFC']
  names(geneList) = as.character(eg[,'ENTREZID'])
  #GSEA分析——GO
  Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  #GSEA分析——KEGG
  KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  #GSEA分析——Reactome
  Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  if (!is.na(Go_Reactomeresult@result[1,1])) {
    png("all1.png",width = 2000,height = 2000,res = 216)
    print(gseaplot2(Go_Reactomeresult,1:min(10,length(Go_Reactomeresult@result$ID)),pvalue_table = TRUE))
    dev.off()
    pdf("all1.pdf",width = 10,height = 10)
    print(gseaplot2(Go_Reactomeresult,1:min(10,length(Go_Reactomeresult@result$ID)),pvalue_table = TRUE))
    dev.off()
    
    for (i in 1:min(10,length(Go_Reactomeresult@result$ID))){
      png(paste(as.character(i),".png",sep=''),width = 2000,height = 2000,res = 216)
      print(gseaplot2(Go_Reactomeresult,i,pvalue_table = TRUE))
      dev.off()
      pdf(paste(as.character(i),".pdf",sep=''),width = 10,height = 10)
      print(gseaplot2(Go_Reactomeresult,i,pvalue_table = TRUE))
      dev.off()
    }
  }
  
  setwd('../')
  
  dir.create(control.type)
  setwd(paste('./',control.type,sep=''))
  x = markers[Control,c('gene','avg_log2FC')]
  tmp=gsub("\\..*","",x$gene)
  eg=bitr(tmp,fromType="SYMBOL" ,toType=c("ENTREZID","ENSEMBL"),OrgDb=OrgDb1)
  genelist <- eg$ENTREZID
  tampt = 1
  eg$logFC=NA
  eg <- eg[!duplicated(eg$SYMBOL), ]
  for (w in 1:length(x$gene)) {
    if (x[w,'gene']==eg[tampt,'SYMBOL']){
      eg[tampt,'logFC']=as.numeric(x[w,'avg_log2FC'])
      tampt = tampt+1
    }
    if(tampt>length(eg$SYMBOL)){break}
  }
  eg <- eg[order(eg$logFC, decreasing = TRUE), ]
  geneList = eg[,'logFC']
  names(geneList) = as.character(eg[,'ENTREZID'])
  #GSEA分析——GO
  Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  #GSEA分析——KEGG
  KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  #GSEA分析——Reactome
  Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  if (!is.na(Go_Reactomeresult@result[1,1])) {
    png("all1.png",width = 2000,height = 2000,res = 216)
    print(gseaplot2(Go_Reactomeresult,1:min(10,length(Go_Reactomeresult@result$ID)),pvalue_table = TRUE))
    dev.off()
    pdf("all1.pdf",width = 10,height = 10)
    print(gseaplot2(Go_Reactomeresult,1:min(10,length(Go_Reactomeresult@result$ID)),pvalue_table = TRUE))
    dev.off()
    
    for (i in 1:min(10,length(Go_Reactomeresult@result$ID))){
      png(paste(as.character(i),".png",sep=''),width = 2000,height = 2000,res = 216)
      print(gseaplot2(Go_Reactomeresult,i,pvalue_table = TRUE))
      dev.off()
      pdf(paste(as.character(i),".pdf",sep=''),width = 10,height = 10)
      print(gseaplot2(Go_Reactomeresult,i,pvalue_table = TRUE))
      dev.off()
    }
  }
  setwd('../')
  
  setwd('../')
}

setwd('../')
