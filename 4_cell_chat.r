library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(scRNAseq)
library(CellChat)
library(monocle)
library(SeuratDisk)
setwd('C:/Users/Zjy52/Desktop/single_cell/new-data/GSE117988')
##ser thread number 
future::plan("multiprocess", workers = 4)
##Cell DB
CellChatDB <-CellChatDB.human
##load data
load('../result/renamed_umap.RData')

##input data
data_input=rename_UMAP@assays$RNA@data
##meta
meta=read.csv('../result/UMAP_dataframe_re.csv')
names=meta[,1]
meta=meta[,-1]
rownames(meta)=names
colnames(meta)=c('res0.5','cluster','label')
##create cell chat 
cellchat=createCellChat(object = data_input,meta=meta,group.by = 'label')
cellchat <- addMeta(cellchat, meta = meta)
# set "labels" as default cell identity
cellchat <- setIdent(cellchat, ident.use = "label") 
# number of cells in each cell group
groupSize <- as.numeric(table(cellchat@idents)) 
##data prepare
CellChatDB.use <- subsetDB(CellChatDB,search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


##caculate interaction
cellchat <- computeCommunProb(cellchat)
##data fliter
cellchat <- filterCommunication(cellchat, min.cells = 10)
inter_data=cellchat@var.features$features.info

##dfnet data
df.net <- subsetCommunication(cellchat)


##caculate cell interaction level
cellchat <- computeCommunProbPathway(cellchat)


##cell chat plot

cellchat <- aggregateNet(cellchat)
##save cellcaht data
#save(cellchat,file='../result/cellchat.Rdata')
par(mfrow = c(1,2), xpd=TRUE)
##circle plot
##by count
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T,
                 label.edge= F,
                 title.name = "Number of interactions")
#by weight
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")




##Hierarchy plot
##select a pathway
vertex.receiver = seq(1,2,4,6)
pathways.show <- c("MIF")
netVisual_aggregate(cellchat, vertex.receiver=vertex.receiver,signaling = pathways.show)


par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

##savedata
write.csv(inter_data,'../result/cell_feature.csv',row.names = F)
write.csv(df.net,'../result/cell_signal.csv',row.names = F)

