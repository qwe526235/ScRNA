library(Seurat)
library(ggplot2)
library(scales)
setwd('C:/Users/Zjy52/Desktop/single_cell/data')
load('../seurat_qc.RData')
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s_genes<-CaseMatch(s.genes,match=rownames(scedata))
g2m_genes<-CaseMatch(g2m.genes,match=rownames(scedata))
scedata <- CellCycleScoring(scedata, 
                            s.features = s_genes, 
                            g2m.features = g2m_genes,
                            set.ident = TRUE)
p1<-VlnPlot(scedata,features = c("S.Score","G2M.Score"),group.by = "orig.ident")
pdf('../5_cell_period_plot/5_cell_period_vln.pdf',width = 7,height = 5)
p1
dev.off()

p2<-ggplot(scedata@meta.data,aes(G2M.Score,S.Score,color = Phase))+
  geom_point()
pdf('../5_cell_period_plot/5_cell_period_point.pdf',width = 7,height = 5)
p2
dev.off()

p3<- ggplot(scedata@meta.data,aes(x=old.ident,fill= Phase))+
  geom_bar(position ="fill")+
  xlab(NULL)+ylab("PERCENT")+
  scale_y_continuous(labels = percent,breaks=seq(0,1,by=0.2))
pdf('../5_cell_period_plot/5_cell_period_bar.pdf',width = 7,height = 5)
p3
dev.off()

p4 <- ggplot(scedata@meta.data,aes(x=orig.ident,fill= Phase))+
  geom_bar(position ="stack")+xlab(NULL)+ylab("Cell")
pdf('../5_cell_period_plot/5_cell_period_bar_2.pdf',width = 7,height = 5)
p4
dev.off()
