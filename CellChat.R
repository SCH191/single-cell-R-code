library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(patchwork)
library(NMF)
library(tidyverse)

ABC = FCD
ABC

###################################

data.input  <- ABC@assays$RNA@data

identity = data.frame(group =ABC$manucelltype , row.names = names(ABC$manucelltype)) 
# create a dataframe consisting of the cell labels

unique(identity$group) 
# check the cell labels

#创建CellChat对象
cellchat <- createCellChat(object = data.input)
cellchat

cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

###########################################

CellChatDB <- CellChatDB.human 
#(out3 <- capture.output(str(CellChatDB)))
#out4 <- paste(out3, collapse="\n")
#mm(gsub("\\$","# ",gsub("\\.\\. ","#",out4)),type ="text")

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object

unique(CellChatDB$interaction$annotation)

##############################################预处理

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  

##############################################相互作用推断

cellchat <- computeCommunProb(cellchat)  
#注意这个函数如果你可以用就用，这个是作者的。

#mycomputeCommunProb <-edit(computeCommunProb)  # computeCommunProb内部似乎有一些bug，同一套数据在window10上没事，到了Linux上有报错。发现是computeExpr_antagonist这个函数有问题，(matrix(1, nrow = 1, ncol = length((group))))，中应为(matrix(1, nrow = 1, ncol = length(unique(group))))？ 不然矩阵返回的不对。de了它。
#environment(mycomputeCommunProb) <- environment(computeCommunProb)
#cellchat <- mycomputeCommunProb(cellchat)  # 这儿是我de过的。


cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#看看这结果
cellchat@netP$pathways
head(cellchat@LR$LRsig)

######################################################画图
setwd("G:/Cellchat")
netVisual_circle(cellchat@net$count, vertex.weight=groupSize,weight.scale=T,label.edge = F,title.name = "Number of interactions")


cellchat = aggregateNet(cellchat)
groupSize = as.numeric(table(cellchat@idents))
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight=groupSize,weight.scale=T,label.edge = F,title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight=groupSize,weight.scale=T,label.edge = F,title.name = "Interaction weights/strength")



####检查每种细胞发出的信号
#互作数量
mat = cellchat@net$count
par(mfrow=c(3,3),xpd=TRUE)
for(i in 1:nrow(mat)){
  mat2=matrix(0,nrow=nrow(mat),ncol=ncol(mat),dimnames=dimnames(mat))
  mat2[i,]=mat[i,]
  netVisual_circle(mat2,vertex.weight=groupSize,weight.scale=T,arrow.width=0.2,arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

#互作强度/概率
mat = cellchat@net$weight
par(mfrow=c(3,3),xpd=TRUE)
for(i in 1:nrow(mat)){
  mat2=matrix(0,nrow=nrow(mat),ncol=ncol(mat),dimnames=dimnames(mat))
  mat2[i,]=mat[i,]
  netVisual_circle(mat2,vertex.weight=groupSize,weight.scale=T,arrow.width=0.2,arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}



#############################################################单个配体信号通路###

levels(cellchat@idents)
cellchat@netP$pathways

pathways.show <- cellchat@netP$pathways[2]
vertex.receiver = c(1,2,3,4,5,6)
#层次图
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
#网络图
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
#和弦图
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#热图
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = c("Blue","RED") )







#多个配体-受体介导的细胞互作关系可视化
P=netVisual_bubble(cellchat, remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble.png",P,width=10,height=20)




##########################################################
cellchat = netAnalysis_computeCentrality(cellchat,slot.name = "netP")

netAnalysis_signalingRole_network(cellchat,signaling = "PTN",width = 15,height = 6,font.size = 10)

ht1 = netAnalysis_signalingRole_heatmap(cellchat,pattern = "outgoing",font.size = 8)
ht2 = netAnalysis_signalingRole_heatmap(cellchat,pattern = "incoming",font.size = 8)
ht1+ht2

selectK(cellchat,pattern = "outgoing")
nPatterns=3
cellchat = identifyCommunicationPatterns(cellchat,pattern = "outgoing",k=nPatterns)
netAnalysis_river(cellchat,pattern = "outgoing")
netAnalysis_dot(cellchat,pattern = "outgoing")

selectK(cellchat,pattern = "incoming")
nPatterns=3
cellchat = identifyCommunicationPatterns(cellchat,pattern = "outgoing",k=nPatterns)
netAnalysis_river(cellchat,pattern = "outgoing")
netAnalysis_dot(cellchat,pattern = "outgoing",dot.size = c(1,10))






#单组细胞###############################################################
seurat_object = ABC
ident_group_by = 'manusubcelltype'
object_species = 'human'
#
dir.create('cellchat0')
setwd('./cellchat0')


#ABC =sample1
#
eval(parse(text = paste("ABC@meta.data$",ident_group_by,"=droplevels(ABC@meta.data$",ident_group_by,")",sep = '')))
data.input  <- ABC@assays$RNA@data
eval(parse(text = paste("identity = data.frame(group =ABC$",ident_group_by," , row.names = names(ABC$",ident_group_by,")) ",sep = '')))
#identity = data.frame(group =ABC$manusubcelltype , row.names = names(ABC$manusubcelltype)) 
unique(identity$group) 
cellchat <- createCellChat(object = data.input)
cellchat
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
eval(parse(text = paste("CellChatDB <- CellChatDB.",object_species,sep = '')))
#CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
eval(parse(text = paste("cellchat <- projectData(cellchat, PPI.",object_species,")",sep = '')))
#cellchat <- projectData(cellchat, PPI.human)  
cellchat <- computeCommunProb(cellchat)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat = netAnalysis_computeCentrality(cellchat,slot.name = "netP")

#
h1 = netVisual_heatmap(cellchat)
h2 = netVisual_heatmap(cellchat,measure = 'weight')
png("1.png",width = 4000,height = 2000,res = 216)
h1+h2
dev.off()

for (j in 1:length(cellchat@netP$pathways)){
  
pathways.show <- cellchat@netP$pathways[j]


#try(ggsave(paste(pathways.show,"_2.png",sep = ''),netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.2),width=10,height=10))
png(paste(pathways.show,"_2.png",sep = ''),width = 2000,height = 2000,res = 216)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

try(ggsave(paste(pathways.show,"_4.png",sep = ''),netAnalysis_contribution(cellchat, signaling = pathways.show),width=10,height=10))
#png(paste(pathways.show,"_4.png",sep = ''),width = 2000,height = 2000,res = 216)
#netAnalysis_contribution(cellchat, signaling = pathways.show)
#dev.off()

#try(ggsave(paste(pathways.show,"_1.png",sep = ''),netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 12, height = 2.5, font.size = 10),width=10,height=10))
png(paste(pathways.show,"_1.png",sep = ''),width = 2000,height = 2000,res = 216)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 12, height = 2.5, font.size = 10)
dev.off()

#try(ggsave(paste(pathways.show,"_3.png",sep = ''),netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", signaling.name = pathways.show),width=10,height=10))
png(paste(pathways.show,"_3.png",sep = ''),width = 2000,height = 2000,res = 216)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", signaling.name = pathways.show)
dev.off()
}

setwd('../')

#两个Cell Chat对象对比################################################################
#eval(parse(text = paste("",sep = '')))

seurat_object = ABC
group_by = 'Group'
group1 = 'Control'
group2 = 'HS'
ident_group_by = 'manusubcelltype'
object_species = 'human'   #human   mouse
#
dir.create('cellchat1')
setwd('./cellchat1')
eval(parse(text = paste("brain = SplitObject(seurat_object,split.by = '",group_by,"')",sep = '')))
rm(seurat_object)
eval(parse(text = paste("sample1 = brain[['",group1,"']]",sep = '')))
eval(parse(text = paste("sample2 = brain[['",group2,"']]",sep = '')))
#sample1 = brain[['Youth']]
#sample2 = brain[['Infancy']]
rm(brain)
gc()
tampt =sample1
#
eval(parse(text = paste("tampt@meta.data$",ident_group_by,"=droplevels(tampt@meta.data$",ident_group_by,")",sep = '')))
data.input  <- tampt@assays$RNA@data
eval(parse(text = paste("identity = data.frame(group =tampt$",ident_group_by," , row.names = names(tampt$",ident_group_by,")) ",sep = '')))
#identity = data.frame(group =tampt$manusubcelltype , row.names = names(tampt$manusubcelltype)) 
unique(identity$group) 
cellchat <- createCellChat(object = data.input)
cellchat
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
eval(parse(text = paste("CellChatDB <- CellChatDB.",object_species,sep = '')))
#CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
eval(parse(text = paste("cellchat <- projectData(cellchat, PPI.",object_species,")",sep = '')))
#cellchat <- projectData(cellchat, PPI.human)  
cellchat <- computeCommunProb(cellchat)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat = netAnalysis_computeCentrality(cellchat,slot.name = "netP")
#
health.cellchat = cellchat



tampt =sample2
#
eval(parse(text = paste("tampt@meta.data$",ident_group_by,"=droplevels(tampt@meta.data$",ident_group_by,")",sep = '')))
data.input  <- tampt@assays$RNA@data
eval(parse(text = paste("identity = data.frame(group =tampt$",ident_group_by," , row.names = names(tampt$",ident_group_by,")) ",sep = '')))
#identity = data.frame(group =tampt$manusubcelltype , row.names = names(tampt$manusubcelltype)) 
unique(identity$group) 
cellchat <- createCellChat(object = data.input)
cellchat
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
eval(parse(text = paste("CellChatDB <- CellChatDB.",object_species,sep = '')))
#CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
eval(parse(text = paste("cellchat <- projectData(cellchat, PPI.",object_species,")",sep = '')))
#cellchat <- projectData(cellchat, PPI.human)  
cellchat <- computeCommunProb(cellchat)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat = netAnalysis_computeCentrality(cellchat,slot.name = "netP")
#
FCD.cellchat = cellchat

###########合并cellchat对象
#group.new = levels(as.factor(merge(brain[[1]],brain[[2]])@meta.data$manusubcelltype))
#group.new = levels(as.factor(merge@meta.data$manusubcelltype))

#eval(parse(text = paste("group.new = levels(as.factor(merge(brain[[1]],brain[[2]])@meta.data$",ident_group_by,"))",sep = '')))
#health.cellchat = liftCellChat(health.cellchat, group.new = group.new)
#FCD.cellchat = liftCellChat(FCD.cellchat, group.new = group.new)


#cellchat.list = c(Youth = health.cellchat,Infancy = FCD.cellchat)
eval(parse(text = paste("cellchat.list = c(",group1," = health.cellchat,",group2," = FCD.cellchat)",sep = '')))
saveRDS(cellchat.list,'cellchatlist.rds')
cellchat = mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = T)


#1.总体通讯数量强度对比
gg1 = compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = 'count')
gg2 = compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = 'weight')
p = gg1+gg2

try(ggsave("1.png",gg1+gg2,width=20,height=10))

#2.数量与强度差异网络图
png("2.png",width = 4000,height = 2000,res = 216)
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = 'weight')
dev.off()
#3.数量与强度差异热图
h1 = netVisual_heatmap(cellchat)
h2 = netVisual_heatmap(cellchat,measure = 'weight')
png("3.png",width = 4000,height = 2000,res = 216)

netVisual_heatmap(cellchat)+netVisual_heatmap(cellchat,measure = 'weight')
dev.off()
#4.细胞数量互作对比网络图
#weight.max = getMaxWeight(cellchat.list,attribute = c('idents','count'))


#5.保守和特异性信号通路的识别与可视化
gg1 = rankNet(cellchat,mode = 'comparison' , stacked = T,do.stat = T)
gg2 = rankNet(cellchat,mode = 'comparison' , stacked = F,do.stat = T)
p = gg1+gg2
try(ggsave("4.png",gg1+gg2,width=20,height=10))

eval(parse(text = paste("cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c('",group1,"', '",group2,"'))",sep = '')))
#cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Youth", "Infancy")) # set factor level

eval(parse(text = paste("pathways = cellchat@netP$",group1,"$pathways",sep = '')))
#pathways = cellchat@netP$Youth$pathways
for (j in 1:length(pathways)){
  pathways.show <- pathways[j] 
  
  try(ggsave(paste(pathways.show,"_1.png"),plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T),width=14,height=10))
  
  png(paste(pathways.show,"_2.png"),width = 4000,height = 2000,res = 216)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
  
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[1]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[1]))))
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[2]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[2]))))
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[1]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[1]))+netAnalysis_contribution(cellchat.list[[2]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[2]))))
  
  png(paste(pathways.show,"_3.png"),width = 4000,height = 2000,res = 216)
  try(weight.max <- getMaxWeight(cellchat.list, slot.name = c("netP"), attribute = pathways.show)) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
}
eval(parse(text = paste("pathways = cellchat@netP$",group2,"$pathways",sep = '')))
#pathways = cellchat@netP$Infancy$pathways
for (j in 1:length(pathways)){
  pathways.show <- pathways[j] 
  
  try(ggsave(paste(pathways.show,"_1.png"),plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T),width=14,height=10))
  
  png(paste(pathways.show,"_2.png"),width = 4000,height = 2000,res = 216)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
  
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[1]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[1]))))
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[2]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[2]))))
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[1]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[1]))+netAnalysis_contribution(cellchat.list[[2]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[2]))))
  
  png(paste(pathways.show,"_3.png"),width = 4000,height = 2000,res = 216)
  try(weight.max <- getMaxWeight(cellchat.list, slot.name = c("netP"), attribute = pathways.show)) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
}



eval(parse(text = paste("group.new = levels(as.factor(seurat_object@meta.data$",ident_group_by,"))",sep = '')))
health.cellchat = liftCellChat(health.cellchat, group.new = group.new)
FCD.cellchat = liftCellChat(FCD.cellchat, group.new = group.new)


#cellchat.list = c(Youth = health.cellchat,Infancy = FCD.cellchat)
eval(parse(text = paste("cellchat.list = c(",group1," = health.cellchat,",group2," = FCD.cellchat)",sep = '')))

cellchat = mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = T)


#1.总体通讯数量强度对比
gg1 = compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = 'count')
gg2 = compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = 'weight')
p = gg1+gg2

try(ggsave("1.png",gg1+gg2,width=20,height=10))

#2.数量与强度差异网络图
png("2.png",width = 4000,height = 2000,res = 216)
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = 'weight')
dev.off()
#3.数量与强度差异热图
h1 = netVisual_heatmap(cellchat)
h2 = netVisual_heatmap(cellchat,measure = 'weight')
png("3.png",width = 4000,height = 2000,res = 216)

netVisual_heatmap(cellchat)+netVisual_heatmap(cellchat,measure = 'weight')
dev.off()
#4.细胞数量互作对比网络图
#weight.max = getMaxWeight(cellchat.list,attribute = c('idents','count'))


#5.保守和特异性信号通路的识别与可视化
gg1 = rankNet(cellchat,mode = 'comparison' , stacked = T,do.stat = T)
gg2 = rankNet(cellchat,mode = 'comparison' , stacked = F,do.stat = T)
p = gg1+gg2
try(ggsave("4.png",gg1+gg2,width=20,height=10))


setwd('../')



#over#           #over#       #over#

#补充########################################
# cellchat.NL <- readRDS("/Users/suoqinjin/Library/CloudStorage/OneDrive-Personal/works/CellChat/tutorial/cellchat_humanSkin_NL.rds")
# cellchat.LS <- readRDS("/Users/suoqinjin/Library/CloudStorage/OneDrive-Personal/works/CellChat/tutorial/cellchat_humanSkin_LS.rds")
# object.list <- list(NL = cellchat.NL, LS = cellchat.LS)
# cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat = cellchat
object.list = cellchat.list
cellchat


## 比较交互总数和交互强度
# 为了回答细胞间通讯是否增强的问题，CellChat 比较了不同生物条件下推断的细胞间通讯网络的相互作用总数和相互作用强度。
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


##比较不同细胞群之间的相互作用数量和相互作用强度
#（A）圆形图显示两个数据集中不同细胞群之间的相互作用数量或相互作用强度的差异
#两个数据集之间细胞间通信网络中相互作用的差异数量或相互作用强度可以使用圆形图可视化，其中红色的
#（或者蓝色的) 彩色边缘代表增加（或者减少）与第一个数据集相比，第二个数据集中的信号。
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# (B) 热图显示两个数据集中不同细胞群之间相互作用或相互作用强度的差异数量
# CellChat 还可以使用热图更详细地显示交互数量或交互强度的差异。顶部彩色条形图表示热图中显示的每列绝对值的总和（传入信号）
# 右侧的彩色条形图表示每行绝对值（传出信令）的总和。因此，条形高度表示两个条件之间的相互作用数量或相互作用强度的变化程度。
# 在颜色栏中，红色的（或者蓝色的） 代表 增加（或者 减少) 与第一个数据集相比，第二个数据集中的信号。
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

# （C）圆形图显示多个数据集中不同细胞群之间的相互作用数量或相互作用强度
# 上面的差异网络分析只针对成对的数据集，如果有更多的数据集可供比较，CellChat 可以直接展示每个数据集中任意两个细胞群体之间的相互作用数量或相互作用强度。
# 为了更好地控制跨不同数据集的推断网络的节点大小和边权重，CellChat 计算每个细胞组的最大细胞数以及所有数据集的最大交互数（或交互权重）。
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# (D) 圆形图显示粗细胞类型之间相互作用的差异数量或相互作用强度
# 为了简化复杂的网络并深入了解细胞类型级别的细胞间通信，CellChat 根据定义的细胞组聚合细胞间通信。这里我们将细胞群分为三种细胞类型，然后重新合并 CellChat 对象的列表。
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#然后，我们可以显示每个数据集中任意两种细胞类型之间的相互作用数量或相互作用强度。
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
#类似地，CellChat 还可以使用圆形图显示任意两种细胞类型之间的相互作用数量或相互作用强度差异。红色（或蓝色）边缘表示第二个数据集与第一个数据集相比信号增加（或减少）。
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)


## 比较 2D 空间中的主要源和目标
# 通过比较 2D 空间中的传出和传入交互强度，可以轻松识别不同数据集之间发送或接收信号发生显着变化的细胞群。
# 通过以下选项 A 识别不同数据集之间发送或接收信号发生显着变化的细胞群，并通过以下选项 B 识别特定细胞群的信号变化。

#（A）识别在发送或接收信号方面发生显著变化的细胞群
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

# (B) 识别特定细胞群的信号传导变化
# 此外，我们可以识别NL和LS之间Inflam.DC和cDC1的特定信号传导变化。signaling.exclude指排除那些信号通路
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Micro_1", signaling.exclude = c(""))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Micro_2", signaling.exclude = c(""))
patchwork::wrap_plots(plots = list(gg1,gg2))



### 第二部分：识别具有不同网络架构和交互强度的改变的信令
# 然后，CellChat 可以根据多个生物条件下的功能或结构相似性来识别具有较大（或较小）网络差异的信号传导以及信号传导组。

## 根据功能/结构相似性识别差异较大（或较小）的信号网络以及信号组
# CellChat 根据不同条件下的功能和拓扑相似性对推断的通信网络进行联合流形学习和分类。注意：这种分析适用于两个以上的数据集。
# 通过量化不同条件下信号通路的细胞通信网络之间的相似性，该分析突出了潜在改变的信号通路。 CellChat 采用了网络生物学中的网络重新布线的概念，并假设不同通信网络之间的差异可能会影响不同条件下的生物过程。 UMAP 用于可视化信号关系并以直观的方式解释我们的信号输出，而不涉及条件分类。
# 功能相似性：高度的功能相似性表明主要的发送者和接收者相似，并且可以解释为两个信号通路或两个配体-受体对表现出相似和/或冗余的作用。注意：功能相似性分析不适用于具有不同细胞类型组成的多个数据集。
# 结构相似性：结构相似性用来比较它们的信令网络结构，不考虑发送者和接收者的相似性。注意：结构相似性分析适用于具有相同细胞类型组成或细胞类型组成截然不同的多个数据集。
# 在这里，我们可以基于功能相似性运行流形和分类学习分析，因为两个数据集具有相同的细胞类型组成。

# 根据功能相似性识别信号基团
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

#根据结构相似性识别信号基团
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)


# 计算并可视化学习的联合流形中的路径距离
# CellChat可以根据共享二维空间中的欧氏距离来识别差异较大（或较小）的信号网络。距离越大意味着两个数据集之间的通信网络在功能或结构相似性方面差异越大。
#应该注意的是，我们只计算两个数据集之间重叠信号通路的距离。这里不考虑仅在一个数据集中识别的那些信号通路。如果数据集超过三个，可以通过修改 comparison函数中的参数进行两两比较rankSimilarity。
rankSimilarity(cellchat, type = "functional")


## 识别具有不同相互作用强度的改变的信号传导
# 通过比较每个信号通路的信息流/相互作用强度，CellChat 通过在一种条件下改变其信息流来识别以下信号通路：(i) 关闭、(ii) 减少、(iii) 开启或 (iv) 增加与另一种情况相比。
# 通过遵循选项 A，根据整体信息流，并通过遵循选项 B，根据传出（或传入）信号模式，识别改变的信号传导途径或配体-受体对。

# (A) 比较每个信号通路或配体-受体对的整体信息流
# CellChat 可以通过简单地比较每个信号通路的信息流来识别保守的和上下文特定的信号通路，该信号通路由推断网络中所有细胞组对之间的通信概率之和定义（即网络中的总权重） ）。
# 此条形图可以以堆叠模式绘制，也可以不以堆叠模式绘制。根据推断网络内 NL 和 LS 皮肤之间整体信息流的差异对显著信号通路进行排序。
# 设置时do.stat = TRUE，进行配对 Wilcoxon 检验以确定两种情况之间的信号信息流是否存在显著差异。顶部红色信号通路在 NL 皮肤中富集，而这些绿色信号通路在 LS 皮肤中富集。
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

# (B) 比较与每个细胞群相关的传出（或传入）信号传导模式
# 上述分析总结了来自传出和传入信令的信息。 CellChat 还可以比较两个数据集之间的传出（或传入）信号传导模式，从而识别表现出不同信号传导模式的信号传导途径/配体受体。我们可以组合来自不同数据集的所有已识别的信号通路，从而并排比较它们，包括传出信号、传入信号以及通过将传出和传入信号聚合在一起的整体信号。
# CellChat 使用热图显示信号（信号通路或配体-受体对）对细胞组的贡献，包括传出信号或传入信号。在此热图中，colobar 表示跨细胞组的信号通路的相对信号强度（请注意，值是按行缩放的）。顶部彩色条形图通过汇总热图中显示的所有信号通路来显示细胞组的总信号强度。右侧灰色条形图通过汇总热图中显示的所有细胞组来显示信号通路的总信号强度。
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
#outgoing
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 18)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 18)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 18, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 18, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#all
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 18, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 18, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



### 第三部分：识别上调和下调的信号配体-受体对
## 通过比较通信概率来识别功能失调的信号
# CellChat 可以比较某些细胞组与其他细胞组的 LR 对介导的通信概率。这可以通过comparison在函数中 设置来完成netVisual_bubble。
netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(1:2),  comparison = c(1, 2), angle.x = 45)

# 此外，与另一数据集相比，CellChat 可以识别一个数据集中上调（增加）和下调（减少）的信号配体-受体对。
# 这可以通过在函数中 指定max.datasetand来完成 。增加的信令意味着与第一数据集相比，这些信令在第二数据集中具有更高的通信概率（强度）。
# 气泡图中显示的配体-受体对可以通过访问。min.datasetnetVisual_bubblegg1$data
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
gg1 + gg2



# 可视化已识别的上调和下调信号配体-受体对
# CellChat 可以使用气泡图（选项 A）、弦图（选项 B）或词云（选项 C）将已识别的上调和下调信号配体-受体对进行可视化。

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = group2
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")
# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05)#, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)#, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = group1,ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = group2,ligand.logFC = -0.05, receptor.logFC = NULL)

# （A）气泡图
# 然后，我们使用气泡图或弦图可视化上调和下调的信号配体-受体对。
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:2), targets.use =c(1:2), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:2), targets.use =c(1:2), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# (B) 和弦图
# 使用弦图可视化上调和下调的信号配体-受体对
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(1:1), targets.use = c(1:1), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(1:1), targets.use = c(1:1), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# （C）词云图
# 使用 wordcloud 可视化一种条件下与另一种条件下富集的配体、信号传导或配体-受体对的比较
# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)
# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)















#批处理#################################
setwd("G:/work/BC&DN/new/cellchat")
dir.create(as.character(i))
setwd(paste('./',as.character(i),sep=''))

#两个Cell Chat对象对比

sample1 = brain[['Control']]#Control
sample2 = brain[['Focus']]#Focus

sample1@meta.data$manusubcelltype = droplevels(sample1@meta.data$manusubcelltype)
sample2@meta.data$manusubcelltype = droplevels(sample2@meta.data$manusubcelltype)

ABC =sample1
#
data.input  <- ABC@assays$RNA@data
identity = data.frame(group =ABC$manusubcelltype , row.names = names(ABC$manusubcelltype)) 
unique(identity$group) 
cellchat <- createCellChat(object = data.input)
cellchat
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  
cellchat <- computeCommunProb(cellchat)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat = netAnalysis_computeCentrality(cellchat,slot.name = "netP")
#
health.cellchat = cellchat



ABC =sample2
#
data.input  <- ABC@assays$RNA@data
identity = data.frame(group =ABC$manusubcelltype , row.names = names(ABC$manusubcelltype)) 
unique(identity$group) 
cellchat <- createCellChat(object = data.input)
cellchat
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  
cellchat <- computeCommunProb(cellchat)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat = netAnalysis_computeCentrality(cellchat,slot.name = "netP")
#
FCD.cellchat = cellchat


#合并cellchat对象
#group.new = levels(as.factor(merge(brain[[1]],brain[[2]])@meta.data$manusubcelltype))
group.new = levels(as.factor(merge@meta.data$manusubcelltype))
health.cellchat = liftCellChat(health.cellchat, group.new = group.new)
FCD.cellchat = liftCellChat(FCD.cellchat, group.new = group.new)

cellchat.list = c(Countrol = health.cellchat,Focus = FCD.cellchat)

cellchat = mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = T)


#1.总体通讯数量强度对比
gg1 = compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = 'count')
gg2 = compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = 'weight')
p = gg1+gg2


png("1.png",width = 4000,height = 2000,res = 216)
p
dev.off()

#2.数量与强度差异网络图
png("2.png",width = 4000,height = 2000,res = 216)


par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = 'weight')
dev.off()
#3.数量与强度差异热图
h1 = netVisual_heatmap(cellchat)
h2 = netVisual_heatmap(cellchat,measure = 'weight')
png("3.png",width = 4000,height = 2000,res = 216)


h1+h2
dev.off()
#4.细胞数量互作对比网络图
#weight.max = getMaxWeight(cellchat.list,attribute = c('idents','count'))


#5.保守和特异性信号通路的识别与可视化
gg1 = rankNet(cellchat,mode = 'comparison' , stacked = T,do.stat = T)
gg2 = rankNet(cellchat,mode = 'comparison' , stacked = F,do.stat = T)
p = gg1+gg2
png("4.png",width = 4000,height = 2000,res = 216)


p
dev.off()
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Countrol", "Focus")) # set factor level


for (j in 1:length(cellchat@netP$Countrol$pathways)){
  pathways.show <- cellchat@netP$Countrol$pathways[j] 
  
  try(ggsave(paste(pathways.show,"_1.png"),plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T),width=14,height=10))
  
  png(paste(pathways.show,"_2.png"),width = 4000,height = 2000,res = 216)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
  
  png(paste(pathways.show,"_3.png"),width = 4000,height = 2000,res = 216)
  try(weight.max <- getMaxWeight(cellchat.list, slot.name = c("netP"), attribute = pathways.show)) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
}
for (j in 1:length(cellchat@netP$Focus$pathways)){
  pathways.show <- cellchat@netP$Focus$pathways[j] 
  
  try(ggsave(paste(pathways.show,"_1.png"),plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T),width=14,height=10))
  
  png(paste(pathways.show,"_2.png"),width = 4000,height = 2000,res = 216)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
  
  png(paste(pathways.show,"_3.png"),width = 4000,height = 2000,res = 216)
  try(weight.max <- getMaxWeight(cellchat.list, slot.name = c("netP"), attribute = pathways.show)) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
}

#每个患者组内实验组对照组细胞通讯比较##################################################################
ABC = merge
ABC@meta.data$manusubcelltype = droplevels(as.factor(ABC@meta.data$CellType))

dir.create('cellchat')
setwd('./cellchat')
tmp = SplitObject(ABC,split.by = 'patient')
for (i in 1:length(tmp)){
dir.create(as.character(names(tmp[i])))
setwd(paste('./',as.character(names(tmp[i])),sep=''))


brain = SplitObject(tmp[[i]],split.by = 'Group')
#两个Cell Chat对象对比

sample1 = brain[['Control']]#Control
sample2 = brain[['Focus']]#Focus

sample1@meta.data$manusubcelltype = droplevels(sample1@meta.data$manusubcelltype)
sample2@meta.data$manusubcelltype = droplevels(sample2@meta.data$manusubcelltype)

ABC =sample1
#
data.input  <- ABC@assays$RNA@data
identity = data.frame(group =ABC$manusubcelltype , row.names = names(ABC$manusubcelltype)) 
unique(identity$group) 
cellchat <- createCellChat(object = data.input)
cellchat
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  
cellchat <- computeCommunProb(cellchat)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat = netAnalysis_computeCentrality(cellchat,slot.name = "netP")
#
health.cellchat = cellchat



ABC =sample2
#
data.input  <- ABC@assays$RNA@data
identity = data.frame(group =ABC$manusubcelltype , row.names = names(ABC$manusubcelltype)) 
unique(identity$group) 
cellchat <- createCellChat(object = data.input)
cellchat
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  
cellchat <- computeCommunProb(cellchat)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat = netAnalysis_computeCentrality(cellchat,slot.name = "netP")
#
FCD.cellchat = cellchat


#合并cellchat对象
#group.new = levels(as.factor(merge(brain[[1]],brain[[2]])@meta.data$manusubcelltype))
group.new = levels(as.factor(merge@meta.data$manusubcelltype))
health.cellchat = liftCellChat(health.cellchat, group.new = group.new)
FCD.cellchat = liftCellChat(FCD.cellchat, group.new = group.new)

cellchat.list = c(Countrol = health.cellchat,Focus = FCD.cellchat)

cellchat = mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = T)

try(dev.off())
#1.总体通讯数量强度对比
gg1 = compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = 'count')
gg2 = compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = 'weight')
p = gg1+gg2

try(ggsave("1.png",gg1+gg2,width=20,height=10))

#2.数量与强度差异网络图
png("2.png",width = 4000,height = 2000,res = 216)
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = 'weight')
dev.off()
#3.数量与强度差异热图
h1 = netVisual_heatmap(cellchat)
h2 = netVisual_heatmap(cellchat,measure = 'weight')
png("3.png",width = 4000,height = 2000,res = 216)

netVisual_heatmap(cellchat)+netVisual_heatmap(cellchat,measure = 'weight')
dev.off()
#4.细胞数量互作对比网络图
#weight.max = getMaxWeight(cellchat.list,attribute = c('idents','count'))


#5.保守和特异性信号通路的识别与可视化
gg1 = rankNet(cellchat,mode = 'comparison' , stacked = T,do.stat = T)
gg2 = rankNet(cellchat,mode = 'comparison' , stacked = F,do.stat = T)
p = gg1+gg2
try(ggsave("4.png",gg1+gg2,width=20,height=10))


cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Countrol", "Focus")) # set factor level


for (j in 1:length(cellchat@netP$Countrol$pathways)){
  pathways.show <- cellchat@netP$Countrol$pathways[j] 
  
  try(ggsave(paste(pathways.show,"_1.png"),plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T),width=14,height=10))
  
  png(paste(pathways.show,"_2.png"),width = 4000,height = 2000,res = 216)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
  
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[1]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[1]))))
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[2]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[2]))))
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[1]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[1]))+netAnalysis_contribution(cellchat.list[[2]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[2]))))
  
  png(paste(pathways.show,"_3.png"),width = 4000,height = 2000,res = 216)
  try(weight.max <- getMaxWeight(cellchat.list, slot.name = c("netP"), attribute = pathways.show)) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
}

for (j in 1:length(cellchat@netP$Focus$pathways)){
  pathways.show <- cellchat@netP$Focus$pathways[j] 
  
  try(ggsave(paste(pathways.show,"_1.png"),plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T),width=14,height=10))
  
  png(paste(pathways.show,"_2.png"),width = 4000,height = 2000,res = 216)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
  
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[1]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[1]))))
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[2]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[2]))))
  try(ggsave(paste(pathways.show,"_4.png"),width = 40,height = 20,netAnalysis_contribution(cellchat.list[[1]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[1]))+netAnalysis_contribution(cellchat.list[[2]], signaling = pathways.show, signaling.name = paste(pathways.show, names(cellchat.list)[2]))))
  
  png(paste(pathways.show,"_3.png"),width = 4000,height = 2000,res = 216)
  try(weight.max <- getMaxWeight(cellchat.list, slot.name = c("netP"), attribute = pathways.show)) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(cellchat.list)) {
    try(netVisual_aggregate(cellchat.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat.list)[i])))
  }
  dev.off()
}

setwd('../')
}
setwd('../')


