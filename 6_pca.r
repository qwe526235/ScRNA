library(Seurat)
library(ggplot2)
library(scales)
setwd('C:/Users/Zjy52/Desktop/single_cell/data')
load('../seurat_qc.RData')
scedata <- ScaleData(scedata, vars.to.regress = c("percent.mt", "nFeature_RNA"), verbose = FALSE)
scedata <- RunPCA(scedata, npcs = 50, verbose = FALSE)
scedata <- FindNeighbors(scedata, reduction = "pca", dims = 1:50)
scedata <- FindClusters(scedata, 
                        resolution = seq(from = 0.1, 
                                         to = 1.0, 
                                         by = 0.1))
scedata <- RunUMAP(scedata, reduction = "pca", dims = 1:50)
p1<-DimPlot(scedata,label=T)
pdf('../6_umap_plot/6_umap_plot.pdf',width = 7,height = 7)
p1
dev.off()
save(scedata,file='../sce_pca.RData')

#Idents(scedata) <- "integrated_snn_res.1"
#scedata$seurat_clusters <- scedata@active.ident
#DimPlot(scedata,label = T,split.by = "orig.ident",ncol = 3)


