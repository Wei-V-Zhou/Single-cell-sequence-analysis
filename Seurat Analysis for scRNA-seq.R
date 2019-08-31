# Load libraries
library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(reticulate)

# Load Data and Create a Seurat Object
data<-read.table("m4T1.raw.featurecounts.txt")
counts<-data[,5:ncol(data)]
pbmc<-CreateSeuratObject(
  counts=as.matrix(counts),
  project="10X_PBMC",min.cells=3,min.features=200,assay="RNA")

# Cell QC and Filtering
Mito.genes<-grep(pattern="^mt-",x=rownames(pbmc@assays[["RNA"]]),value=TRUE)
percent.mito<-Matrix::colSums(pbmc@assays[["RNA"]][Mito.genes, ])/Matrix::colSums(pbmc@assays[["RNA"]])
pbmc<-AddMetaData(object=pbmc,metadata=percent.mito,col.name="percent.mito")
pbmc$percent.mito<-percent.mito
VlnPlot(object=pbmc,features=c("nFeature_RNA","nCount_RNA","percent.mito"),ncol=3)
par(mfrow=c(1,2))
FeatureScatter(object=pbmc,feature1="nCount_RNA",feature2="percent.mito")
FeatureScatter(object=pbmc,feature1="nCount_RNA",feature2="nFeature_RNA")
pbmc<-subset(x=pbmc,subset=nFeature_RNA>7000 & nFeature_RNA<12500 & percent.mito>-Inf & percent.mito<0.1)

# Data Normalization
pbmc<-NormalizeData(
  object=pbmc,
  normalization.method="LogNormalize",
  scale.factor=10000)

# Detect Variable Genes
pbmc<-FindVariableFeatures(
  object=pbmc,
  mean.function=ExpMean,
  dispersion.function=LogVMR,
  x.low.cutoff=0.0125,
  x.high.cutoff=3,
  y.cutoff=0.5,
  nfeatures=2000)
head(x=HVFInfo(object=pbmc))

# Data Scaling and Removing
pbmc<-ScaleData(object=pbmc,vars.to.regress=c("nCounts_RNA","percent.mito"))

# Perform linear dimensional reduction
pbmc<-RunPCA(object=pbmc,npcs=30,verbose=FALSE)
DimPlot(object=pbmc,reduction="pca")
FeaturePlot(object=pbmc,features="Cd19")
FeatureScatter(object=pbmc,feature1="Cd14",feature2="PC_1")
FeatureScatter(object=pbmc,feature1="Cd14",feature2="Cd19")
CellScatter(
  object=pbmc,
  cell1="X4T1.1.10",
  cell2="X4T1.1.11")
VariableFeaturePlot(object=pbmc)
VlnPlot(object=pbmc,features=c("Ccl5","Cd14","Cd19"))
RidgePlot(object=pbmc,feature=c("Ccl5","Cd14","Cd19"))
DimHeatmap(object=pbmc,reduction="pca",cells=200,balanced=TRUE)

# Determine significant principal component
pbmc<-JackStraw(object=pbmc,reduction="pca",dims=20,num.replicate=100,prop.freq=0.1,verbose=FALSE)
pbmc<-ScoreJackStraw(object=pbmc,dims=1:20,reuction="pca")
JackStrawPlot(object=pbmc,dims=1:20,reduction="pca")
ElbowPlot(object=pbmc)

# Cell Clustering
pbmc<-FindNeighbors(pbmc,reduction="pca",dims=1:20)
pbmc<-FindClusters(pbmc,resolution=0.5,algorithm=1)

# Perform non-linear dimensional reduction
pbmc<-RunTSNE(object=pbmc,dim.use=1:10,do.fast=TRUE)
DimPlot(object=pbmc,reduction="tsne")

# Run UMAP
# install python package "umap-learn"
pbmc<-RunUMAP(pbmc,reduction="pca",dims=1:20)
DimPlot(pbmc,reduction="umap",pt.size=1)
DimPlot(pbmc,reduction="umap",split.by="seurat_clusters")

# Find Differentially Expressed genes
# Find Cluster1 markers
cluster1.markers<-FindMarkers(object=pbmc,ident.1=1,min.pct=0.25)
print(x=head(x=cluster1.markers))
# Find cluster5 markers distinguishing from clusters 0 and 3
cluster5.markers<-FindMarkers(object=pbmc,ident.1=5,ident.2=c(0,3),min.pct=0.25)
print(x=head(x=cluster5.markers,n=5))
# Find markers for every cluster compared to all other cells
pbmc.markers<-FindAllMarkers(object=pbmc,only.pos=TRUE,min.pct=0.25,thresh.use=0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2,avg_logFC)
cluster2.markers<-FindMarkers(object=pbmc,ident.1=2,thresh.use=0.25,test.use="roc",only.pos=TRUE)
VlnPlot(object=pbmc,features=c("Ighm","Igkc"))
FeaturePlot(object=pbmc,feature=c("Ighm","Igkc"),cols=c("grey","blue"),reduction="umap")
top10<-pbmc.markers %>% group_by(cluster) %>% top_n(10,avg_logFC)
DoHeatmap(object=pbmc,features=top10$gene,label=TRUE)

# Assign the celltype identity to the cluster
current.cluster.ids<-c(0,1,2,3,4,5,6,7,8)
new.cluster.ids<-c("CD4 T Cells","CD8 T Cells","CD14 T Cells","CD4 B Cells","CD8 B Cells","CD14 B Cells","Tumour Cells","Macrophage","Vessel")
pbmc@active.ident<-plyr::mapvalues(x=pbmc@active.ident,from=current.cluster.ids,to=new.cluster.ids)
DimPlot(object=pbmc,reduction="tsne",do.label=TRUE,pt.size=0.5)

# Further subdivide celltypes
pbmc<-StashIdent(object=pbmc,save.name="ClusterNames_0.6")
pbmc<-FindClusters(object=pbmc,reduction.type="pca",dims.use=1:10,resolution=0.8,print.output=FALSE)
plot1<-DimPlot(object=pbmc,reduction="tsne",do.return=TRUE,no.legend=TRUE,do.label=TRUE)
plot2<-DimPlot(object=pbmc,reduction="tsne",do.return=TRUE,group.by="ClusterNames_0.6",no.legend=TRUE,do.label=TRUE)
plot_grid(plot1,plot2)
# Find discriminating markers
FeaturePlot(object=pbmc,features=c("Ccr7","Cd19"),cols=c("green","blue"))
# Save Data
pbmc<-SetIdent(object=pbmc,value="ClusterNames_0.6")
saveRDS(pbmc,file="m4T1.10X_pbmc_final.rds")