library(Seurat)
library(ggplot2)
library(scales)
setwd('C:/Users/Zjy52/Desktop/single_cell/data')
load('../sce_pca.RData')
DefaultAssay(scedata) <- "RNA"
all.markers  <- FindAllMarkers(scedata, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.75)
significant.markers  <- all.markers [all.markers $p_val_adj < 0.05, ]
#write.csv(significant.markers, file = "../significant.markers.csv")#保存
cluster=data.frame(table(significant.markers$cluster))
cluster=cluster[order(cluster$Freq),]
data=significant.markers[significant.markers$cluster %in% cluster$Var1[1],]
for(i in 2:length(cluster$Var1)){
  temp=significant.markers[significant.markers$cluster %in% cluster$Var1[i],]
  temp=temp[order(temp[,2]),]
  n=temp[1:min(cluster$Freq),]
  df=significant.markers[significant.markers$gene %in% n$gene,]
  data=rbind(data,df)
}
marker=data$gene
alldata <- ScaleData(scedata, 
                     features = marker, 
                     assay = "RNA")
p1<-DoHeatmap(alldata, 
          features = marker,
          group.by = "seurat_clusters",
          assay = "RNA",
          )
p1
