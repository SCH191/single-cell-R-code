############################################
res = H1
head(res)

with(res,plot(avg_log2FC,-log10(p_val),pch=20,main="Volcano plot",xlim=c(-2.5,2)))


with(subset(res, p_val_adj<.05 ), 
     points(avg_log2FC, 
            -log10(p-val), 
            pch=20, col="red"))

with(subset(res, abs(avg_log2FC)>1), 
     points(avg_log2FC, -log10(p-val), 
            pch=20, 
            col="orange"))

with(subset(res, p_val_adj<.05 & abs(avg_log2FC)>1), 
     points(avg_log2FC, 
            -log10(p-val), 
            pch=20, 
            col="green"))
# Label points with the textxy function from the calibrate plot
#install.packages("calibrate")
library(calibrate)
with(subset(res, p_val_adj<.05 & abs(avg_log2FC)>1), 
     textxy(avg_log2FC, -log10(p-val), 
            labs=Gene, cex=.8))


############################################################


library(EnhancedVolcano)

data = read.csv(file='test.csv',header = T,row.names = 1)
data = res

EnhancedVolcano(data, lab = data$gene, x = 'avg_log2FC' , y = 'p_val' , xlim = c(-20,20) , ylim = c(0,6) , pCutoff = 1E-50 , FCcutoff = 2)
#单组差异经典火山图#################################

library(ggpubr)
library(ggthemes)

ABC = ABC
pvalue_cut = 0.05
up_cut = 0.5
down_cut = -0.5

#
markers = read.csv('markers.csv',header = T)
#or
DefaultAssay(ABC) = 'RNA'
ABC <- NormalizeData(ABC)
ABC <- FindVariableFeatures(ABC, selection.method = "vst", nfeatures = 10000)
all.genes <- rownames(ABC)
ABC <- ScaleData(ABC, features = all.genes)
markers = FindAllMarkers(ABC, only.pos = F)

Idents(ABC) = ABC@meta.data$manusubcelltype
markers = FindMarkers(ABC,ident.1 = 'nident.1',ident.2 = 'nident.2',only.pos = F)
####
df = subset(markers , subset = cluster=='HS')

####
df$logFC = df$avg_log2FC
df$logP <- -log10(ifelse(df$p_val_adj == 0, 1e-300, df$p_val_adj))
head(df)

ggscatter(df , x = 'logFC',y='logP')

df$Group = 'not-significant'
df$Group[which((df$p_val_adj<pvalue_cut)&(df$logFC>up_cut))] = 'up-regulated'
df$Group[which((df$p_val_adj<pvalue_cut)&(df$logFC<down_cut))] = 'down-regulated'
table(df$Group)

df %>% top_n(n = 10, wt = avg_log2FC) -> top10
df %>% top_n(n = -10, wt = avg_log2FC) -> bellow10


ggscatter(df , x = 'logFC',y='logP',
          color = 'Group',palette = c('#327eba','#BBBBBB','#e06663'),
          size = 1,
          label = df$lebal,font.label = 8,
          repel = T,
          xlab = 'log2FoldChange',ylab = '-log10(Adjust P-value)')+theme_base()+
  geom_hline(yintercept = -log10(pvalue_cut),linetype='dashed')+
  geom_vline(xintercept = c(down_cut,up_cut),linetype='dashed')+
  geom_text_repel(data=bellow10,color = '#327eba',aes(x=logFC,y=logP,label=gene),force = 1.2,arrow = arrow(length = unit(0.008, "npc"),type = "open", ends = "last"))+
  geom_text_repel(data=top10,color = '#e06663',aes(x=logFC,y=logP,label=gene),force = 1.2,arrow = arrow(length = unit(0.008, "npc"),type = "open", ends = "last"))





#多组差异分析火山图####################################################

library(ggplot2)
library(tidyverse)
library(ggrepel)

ABC = ABC

#数据读入：
df <- read.csv("GroupVolcano test2.csv",header = T)
mapping <- setNames(c(0:(length(levels(as.factor(df$cluster)))-1)),levels(as.factor(df$cluster)))
df$num = mapping[df$cluster]
head(df)#查看前六行

#
df = FindAllMarkers(ABC, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
mapping <- setNames(c(0:(length(levels(as.factor(df$cluster)))-1)),levels(as.factor(df$cluster)))
df$num = mapping[df$cluster]
head(df)#查看前六行

#
DefaultAssay(ABC) = 'RNA'
Idents(ABC) = ABC@meta.data$Group
list = SplitObject(ABC,split.by = 'manusubcelltype')
Idents(ABC) = ABC@meta.data$manucelltype
df = data.frame()
for (i in 1:length(list)) {
  tampt = list[[i]]
  marker = FindAllMarkers(tampt, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  marker = subset(marker,subset = cluster=='HS')
  marker$cluster = as.character(names(list[i]))
  marker$num = as.character(i-1)
  df = rbind(df,marker)
}
head(df)#查看前六行


####    ###    ###   ###   ###
head(df)#查看前六行
df$label <- ifelse(df$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
head(df)

top10sig = data.frame()
for (i in 1:length(levels(as.factor(df$cluster)))) {
  sig <- filter(df,cluster==as.character((levels(as.factor(df$cluster))[i]))) %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
  max <- filter(df,cluster==as.character((levels(as.factor(df$cluster))[i]))) %>% distinct(gene,.keep_all = T) %>% top_n(1,(avg_log2FC))
  min <- filter(df,cluster==as.character((levels(as.factor(df$cluster))[i]))) %>% distinct(gene,.keep_all = T) %>% top_n(-1,(avg_log2FC))
  sig$max_logFC = max$avg_log2FC
  sig$min_logFC = min$avg_log2FC
  top10sig = rbind(top10sig,sig)
}
top10sig <- top10sig[order(top10sig$num), ]
df <- df[order(df$num), ]

head(top10sig)

df$size <- case_when(!(df$gene %in% top10sig$gene)~ 1,
                     df$gene %in% top10sig$gene ~ 2)

dt <- filter(df,size==1)
head(dt)

#绘制每个Cluster Top10以外基因的散点火山图：
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)
p

#叠加每个Cluster Top10基因散点(将散点适当放大强调）：
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p

#根据图p中avg_log2FC区间确定背景柱长度：
for (i in 0:(length(levels(as.factor(top10sig$cluster)))-1)) {
  max = top10sig[which(top10sig[,'num']==as.integer(i)),'max_logFC'][1]
  min = top10sig[which(top10sig[,'num']==as.integer(i)),'min_logFC'][1]
  if(i==0){
    max_list = c(max)
    min_list = c(min)
   } else{
     max_list = c(max_list,max)
     min_list = c(min_list,min)
   }
}
dfbar<-data.frame(x=c(0:(length(levels(as.factor(top10sig$cluster)))-1)),
                  y= max_list)
dfbar1<-data.frame(x=c(0:(length(levels(as.factor(top10sig$cluster)))-1)),
                   y= min_list)
#绘制背景柱：
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1

dt$num = as.integer(dt$num)
top10sig$num = as.integer(top10sig$num)


p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = num, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = num, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p2

#添加X轴的cluster色块标签：
dfcol<-data.frame(x=c(0:(length(levels(as.factor(top10sig$num)))-1)),
                  y=0,
                  label=levels(as.factor(top10sig$cluster)))
library(scales)
mycol = hue_pal()(length(levels(as.factor(top10sig$num))))
#mycol <- c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F","#7E61487F")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.4,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3

#给每个Cluster差异表达前Top10基因加上标签：
library(ggrepel)
p4 <- p3+
  geom_text_repel(
    data=top10sig,
    aes(x=num,y=avg_log2FC,label=gene),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )
p4

#散点颜色调整：
p5 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("red","black"))
p5

#修改X/Y轴标题和添加cluster数字：
p6 <- p5+
  labs(x="Cluster",y="average logFC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =6,
            color ="white")
p6

#自定义主题美化：
p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
p7








