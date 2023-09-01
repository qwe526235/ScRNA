library(Seurat)
library(ggplot2)
library(scales)
setwd('C:/Users/Zjy52/Desktop/single_cell/')
data=read.csv('significant.markers.csv')
cluster=data.frame(table(data$cluster))
cluster=cluster[order(cluster$Freq),]
dk=data[data$cluster %in% cluster$Var1[1],]
for(i in 2:length(cluster$Var1)){
  temp=data[data$cluster %in% cluster$Var1[i],]
  temp=temp[order(temp[,2]),]
  n=temp[1:min(cluster$Freq),]
  df=data[data$gene %in% n$gene,]
  dk=rbind(dk,df)
}
for(i in 1:length(cluster$Var1)){
  k=cluster$Var1[i]
  dm=dk[dk$cluster==k,]
  list=paste(dm$gene[1],',')
  for(j in 2:(length(dm$gene)-1)){
    list=paste(list,dm$gene[2],',')
  }
  list=paste(list,dm$gene[length(dm$gene)])
  write.csv(list,paste('cell identity/marker_gene',k,'.csv'))
}