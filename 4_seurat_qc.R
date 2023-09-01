library(Seurat)
setwd('C:/Users/Zjy52/Desktop/single_cell/data')

folders=list.files('C:/Users/Zjy52/Desktop/single_cell/data')
folders
scList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder,
                     min.cells = 3, min.features = 200)
})
PD<- merge(scList[[1]], 
            y = c(scList[[2]]), 
            add.cell.ids = c('PD1','PD2'), 
            project = "PD")
Ctrl<-merge(scList[[3]],
            y=c(scList[[4]]),
            add.cell.ids=c('Ctrl1','Ctrl2'))
##preqc_plot
PD[["percent.mt"]] <- PercentageFeatureSet(PD,pattern = "^MT-")
Ctrl[["percent.mt"]] <- PercentageFeatureSet(Ctrl,pattern = "^MT-")


preQC_PD <- VlnPlot(PD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    ncol = 3, 
                    group.by = "orig.ident", 
                    pt.size = 0)
pdf('../4_preqc_plot/PD_pre_qc.pdf',width = 10,height = 5)
preQC_PD
dev.off()

preQC_Ctrl <- VlnPlot(Ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    ncol = 3, 
                    group.by = "orig.ident", 
                    pt.size = 0)
pdf('../4_preqc_plot/Ctrl_pre_qc.pdf',width = 10,height = 5)
preQC_Ctrl
dev.off()

##QC
##1st
PD_QC<- subset(PD, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
Ct_QC<-subset(Ctrl,subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
sampleList <- list(PD_QC, Ct_QC)
ifnb.list <- lapply(X = sampleList, FUN = function(x) {
  x <- NormalizeData(x,normalization.method='LogNormalize')
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
scedata <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

scedata<- IntegrateData(anchorset = scedata)


save(scedata, file = "../seurat_qc.RData")

