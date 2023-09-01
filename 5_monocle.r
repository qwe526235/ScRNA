library(Seurat)
library(ggplot2)
library(dplyr)
library(SeuratDisk)
library(monocle3)
setwd('C:/Users/Zjy52/Desktop/single_cell/new-data/GSE117988')
#load data
pbmc=load('../result/renamed_umap.RData')
data <- GetAssayData(rename_UMAP, assay ='RNA', slot = 'counts')
cell_metadata <- rename_UMAP@meta.data
gene_annotation <-data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <-rownames(data)
cds <- new_cell_data_set(data,cell_metadata =cell_metadata, gene_metadata =gene_annotation)
##monocel reduction
cds <- preprocess_cds(cds, num_dim = 50)
##umap
umap<- reduce_dimension(cds,preprocess_method = "PCA")
umap=cluster_cells(umap)
umap=learn_graph(umap)
umap=order_cells(umap)
##tsne
#tsne <- reduce_dimension(cds, reduction_method="tSNE",preprocess_method = "PCA")

##seurat position
#cds.embed <- umap@int_colData$reducedDims$UMAP
#int.embed <- Embeddings(rename_UMAP, reduction = "umap")
#int.embed <- int.embed[rownames(cds.embed),]
#umap@int_colData$reducedDims$UMAP <- int.embed

##top 10 exprssion gene
Track_genes <- graph_test(umap,neighbor_graph="principal_graph", cores=4)
Track_genes <-Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-5)
Track_genes_sig <- Track_genes %>%
  top_n(n=10, morans_I) %>% 
  pull(gene_short_name) %>% 
  as.character()
plot_genes_in_pseudotime(umap[Track_genes_sig,],
                         color_cells_by="seurat_clusters",
                         min_expr=0.5, ncol= 2)

