library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(scRNAseq)
library(CellChat)
library(monocle)
library(SeuratDisk)
setwd('C:/Users/Zjy52/Desktop/single_cell/new-data/GSE117988')

##tumor
t_file=list.files('C:/Users/Zjy52/Desktop/single_cell/new-data/GSE117988/pbmc')

file=paste0(getwd(),'/pbmc/',t_file,'/Matrix')
#read data
tulist = lapply(file,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                    project = strsplit(folder,'/')[[1]][length(strsplit(folder,'/')[[1]])-1],
                     min.cells = 3, min.features = 200)
})



##data merge
#dt<- merge(tulist[[1]], 
 #          y = c(tulist[[2]]), 
 #         add.cell.ids = c('SRR7722937','SRR7722938'), 
 #          project = "TUMOR")
slist=tulist[2:(length(t_file))]
#merge_by_group
dt=merge(tulist[[1]],y=c(slist),
         add.cell.ids=t_file,
         project='PBMC')
cell_gene=data.frame(dt@assays$RNA@data)
cell_test=cell_gene[1:20,1:20]
##total_merge
#dtt=merge(tulist[[1]],y=c(slist),
         #project='PBMC')


## add mt cells

dt[['percent.mt']]=PercentageFeatureSet(dt, pattern = "^MT-")

## check the quality dataframe
qc=data.frame(dt@meta.data,5)
qc=qc[,1:4]
qc_plot=VlnPlot(dt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
qc_plot
#outliers
find_outiers=function(x,iqr){
  low=quantile(x,probs=c(0.25))-IQR(x)*iqr
  high=quantile(x,probs=c(0.75))+IQR(x)*iqr
  temp=c(low,high)
}
nf=find_outiers(qc$nFeature_RNA,1.5)
nc=find_outiers(qc$nCount_RNA,1.5)
nm=find_outiers(qc$percent.mt,1.5)

##fliter by outliers or by arg
fliter_dt=subset(dt, subset = nFeature_RNA > 200 & nFeature_RNA < nf[2] & percent.mt < nc[2])
f_plot=VlnPlot(fliter_dt,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
f_plot

##normalize
n_data<- NormalizeData(fliter_dt, 
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)

## FindVariableFeatures
pbmc <- FindVariableFeatures(n_data,
                             selection.method = "vst",
                             nfeatures = 2000)



##top gene
gene_data=data.frame(pbmc@assays$RNA@meta.features)
top5 <- head(VariableFeatures(pbmc), 5)
plot1 <- VariableFeaturePlot(pbmc)
plot1
plot2 <- LabelPoints(plot = plot1, points = top5, repel = TRUE)
plot2




##data_scale
all.genes <- rownames(pbmc)
#pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc)
scale_data=data.frame(pbmc[["RNA"]]@scale.data)
##pca
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#JackStraw_optional
#pbmc <- JackStraw(pbmc, num.replicate = 100)
#jack_data=data.frame(pbmc@reductions$pca@jackstraw$empirical.p.values)
##pbmc <- ScoreJackStraw(pbmc, dims = 1:10)
#ackStrawPlot(pbmc, dims = 1:10)
##pca_select
pca_plot=ElbowPlot(pbmc)
pca_plot
#caculate 85% contribution dimension
pca_ca=Stdev(object = pbmc,reduction = 'pca')
for(i in 1:length(pca_ca)){
  if(sum(pca_ca[1:i])/sum(pca_ca)>=0.85){
    pc_number=i
    break
  }
}

#cell cluster
pbmc <- FindNeighbors(pbmc,dims =1:pc_number)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)



##UNAP
UMAP<- RunUMAP(pbmc,dim=1:pc_number)
u_plot=DimPlot(UMAP, reduction = "umap", label = TRUE)
u_plot
# TSNE
tsne <- RunTSNE(pbmc, dims = 1:pc_number)
t_plot=DimPlot(tsne, reduction = "tsne", label = TRUE)
t_plot
##umap_data
um_data=data.frame(cell=row.names(UMAP@meta.data),
                   res=UMAP@meta.data$RNA_snn_res.0.5,
                   cluster=UMAP@meta.data$seurat_clusters)
##tsne_data
ts_data=data.frame(cell=row.names(tsne@meta.data),
                   res=tsne@meta.data$RNA_snn_res.0.5,
                   cluster=tsne@meta.data$seurat_clusters)

##marker_gene
#total_marker_gene
##umap_marker
##only positive gene
#um_marker= FindAllMarkers(UMAP,only.pos = TRUE, 
   #                           min.pct = 0.25, 
   #                           logfc.threshold = 0.25)

um_marker= FindAllMarkers(UMAP, 
                          min.pct = 0.25, 
                           logfc.threshold = 0.25)
u_marker_plot=um_marker[,c(2,6:7)]
##top-text
top_u_marker=um_marker %>%
  group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) %>%ungroup()
##top-point
u_point=um_marker %>%
  group_by(cluster) %>% top_n(n =100, wt = avg_log2FC) %>% ungroup()
u_ma=ggplot()+geom_violin(data=u_marker_plot,aes(x=cluster,y=avg_log2FC,fill=cluster))+
  #geom_jitter(data=u_point,aes(x=cluster,y=avg_log2FC,color=cluster),alpha=0.5)+
ggrepel::geom_text_repel(data=top_u_marker,aes(x=cluster,y=avg_log2FC,label=gene),
                         size=2.5)


u_ma
##tsene-optional-same as UMAP
data.ls <- data.frame(listDatasets())

##cell annotation
#SingleR reference
reference=HumanPrimaryCellAtlasData()
#ref=SummarizedExperiment(reference)
testdata <- GetAssayData(UMAP, slot="data")
clusters <- UMAP@meta.data$seurat_clusters

pred.grun <- SingleR(test=testdata, 
                     ref=reference,
                     labels=reference$label.main,clusters = clusters)
id=data.frame(cluster=rownames(pred.grun),cell=pred.grun$labels)
new_id=id$cell
names(new_id) <- levels(UMAP)
is.na(new_id)
rename_UMAP<- RenameIdents(UMAP, new_id)
DimPlot(rename_UMAP, reduction = "umap", label = TRUE, pt.size = 0.5)



##save result
##
##renamed UMAP DATA
save(rename_UMAP,file='../result/renamed_umap.RData')
## UMAP Data befroe rename
save(UMAP,file='../result/UMAP_ori.RData')

##no-rename UMAP dataframe
write.csv(um_data,'../result/UMAP_dataframe_nore.csv',row.names = F)
##renameUMAP dataframe
um_data_re=left_join(um_data,id,by=c('cluster'='cluster'))
write.csv(um_data_re,'../result/UMAP_dataframe_re.csv',row.names = F)

##umap ggplot data norename
umap_pos=Embeddings(UMAP,'umap')
umap_data_nore=cbind(um_data,umap_pos)
umap_data_re=left_join(umap_data_nore,id,by='cluster')
write.csv(umap_data_re,'../result/umap_ggplot_data.csv',row.names = F)


##top gene
write.csv(u_marker_plot,'../result/umap_gene.csv',row.names =F )
