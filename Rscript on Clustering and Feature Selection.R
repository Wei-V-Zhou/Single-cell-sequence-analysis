# Load Libraries
library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(reticulate)

# Load Data and Create a Seurat Object
load("m4T1.D4.sceset.qc.new.Rdata")
counts<-counts(sce.all)
pbmc<-CreateSeuratObject(
  counts=as.matrix(counts),
  project="10X_PBMC",min.cells=3,min.features=200,assay="RNA")

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
pbmc<-ScaleData(object=pbmc)

# Perform linear dimensional reduction
pbmc<-RunPCA(object=pbmc,npcs=30,verbose=FALSE)
DimPlot(object=pbmc,reduction="pca")
FeaturePlot(object=pbmc,features="Cd16")
FeatureScatter(object=pbmc,feature1="Cd14",feature2="PC_1")
FeatureScatter(object=pbmc,feature1="Cd14",feature2="Cd19")
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
