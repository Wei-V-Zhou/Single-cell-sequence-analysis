#=======================================================================================#
# The usages of all packages:                                                           #
#      dplyr: provide a flexible grammar of data manipulation                           #
#     export: export the R plot into PPT, Word, jpg, png, tif, pdf, etc. table to excel #
#     Seurat: create an object for single-cell expression data to reduce the dimensions #
#    cowplot: supply a simple code for easy theme change while plotting                 #
#  tidyverse: assemble the packages of data manipulation and visualization              #
# reticulate: supply the R interface to Python for UMAP clustering                      #
#=======================================================================================#

##################
# Load libraries #
##################
library(dplyr)
library(export)
library(Seurat)
library(cowplot)
library(tidyverse)
library(reticulate)

########################################
# Load Data and Create a Seurat Object #
########################################
# Read the txt data using data frame
data<-read.table("m4T1.raw.featurecounts.txt")
# Select the gene counts from five to last columns in data 
counts<-data[,5:ncol(data)]
dim(counts)
# [1] 38538   951
# Create a single cell-seq object with raw data and make a basic minimum gene-cutoff which could
# filter the genes expressed in  at least 3 cells and choose the cells with at least 200 genes
sce<-CreateSeuratObject(
  counts=as.matrix(counts),
  min.cells=3,
  min.features=200,
  assay="scRNA",
  project="Bone metastasis")
dim(sce)
# [1] 30542   949

#########################
# Cell QC and Filtering #
#########################
# Find the mitochondrial genes using "^mt-" pattern in assays with counts matrix
Mito.genes<-grep(pattern="^mt-",x=rownames(sce@assays[["scRNA"]]),value=TRUE)
Mito.genes
# [1] "mt-Tf"   "mt-Rnr1" "mt-Tv"   "mt-Rnr2" "mt-Tl1"  "mt-Nd1"  "mt-Ti"   "mt-Tq"   "mt-Tm"  
# [10] "mt-Nd2"  "mt-Tw"   "mt-Ta"   "mt-Tn"   "mt-Tc"   "mt-Ty"   "mt-Co1"  "mt-Ts1"  "mt-Co2" 
# [19] "mt-Tk"   "mt-Atp8" "mt-Atp6" "mt-Co3"  "mt-Nd3"  "mt-Tr"   "mt-Nd4"  "mt-Th"   "mt-Ts2" 
# [28] "mt-Tl2"  "mt-Nd5"  "mt-Nd6"  "mt-Te"   "mt-Cytb" "mt-Tt"   "mt-Tp"   "mt-Nd4l"
# Calculate the percentage of mitochondrial genes in every cell and store it in percentage.mito using AddMetaData
percent.mito<-Matrix::colSums(sce@assays[["scRNA"]][Mito.genes, ])/Matrix::colSums(sce@assays[["scRNA"]])
sce<-AddMetaData(object=sce,metadata=percent.mito,col.name="percent.mito")
# In case the above function doesn't work, the following step was taken
sce$percent.mito<-percent.mito
head(sce$percent.mito)
# X4T1.1.10  X4T1.1.11  X4T1.1.12  X4T1.1.15  X4T1.1.17  X4T1.1.18 
# 0.02579563 0.01460066 0.02434285 0.01319374 0.02295164 0.03178161
# Show the Violin Plot of the nFeature_RNA, nCount_RNA, and the percent.mito in 3 columns
VlnPlot(object=sce,features=c("nFeature_scRNA","nCount_scRNA","percent.mito"),ncol=3)
par(mfrow=c(1,2))
# GenePlot is used to visualize gene-gene relationships or anything calculated by the object for further filtering
FeatureScatter(object=sce,feature1="nCount_scRNA",feature2="percent.mito",cols="Blue")
FeatureScatter(object=sce,feature1="nCount_scRNA",feature2="nFeature_scRNA")
# Remove the cells expressed the genes (>12500, <5000 or percent.mito>0.1)
sce<-subset(x=sce,subset=nFeature_scRNA>5000 & nFeature_scRNA<12500 
            & percent.mito>-Inf & percent.mito<0.1)
dim(sce)
# [1] 30542   911

######################
# Data Normalization #
######################
# Normalize the data, log-transform the result and mutiply this by a scale factor
sce<-NormalizeData(
  object=sce,
  normalization.method="LogNormalize",
  scale.factor=10000)

#########################
# Detect Variable Genes #
#########################
# Calculate highly variable genes and focuse on these for downstream analysis: Calculate the average expression 
# and dispersion for each gene, place these genes into bins and calculate a z-score for dispersion to control for 
# the relationship between variablity and average expression
sce<-FindVariableFeatures(
  object=sce,
  mean.function=ExpMean,
  dispersion.function=LogVMR,
  x.low.cutoff=0.0125,
  x.high.cutoff=3,
  y.cutoff=0.5)
head(x=HVFInfo(object=sce))

#############################
# Data Scaling and Removing #
#############################
# Scale and center the data to regress out the noise-like siganal, nCounts_scRNA and percent.mito
sce<-ScaleData(object=sce,vars.to.regress=c("nCounts_scRNA","percent.mito"))

########################################
# Perform linear dimensional reduction #
########################################
# Perform the PCA reduction for scaled data
sce<-RunPCA(object=sce,npcs=30,verbose=FALSE)
# Examine and visualize the PCA result
DimPlot(object=sce,reduction="pca")
# Dimensional reduction with the color by a quantitative feature
FeaturePlot(object=sce,features="Cd19")
# Scatter plot across single cells without geneplot
FeatureScatter(object=sce,feature1="Cd14",feature2="PC_1")
FeatureScatter(object=sce,feature1="Cd14",feature2="Cd19")
# Scatter plot across individual features without cellplot
CellScatter(
  object=sce,
  cell1="X4T1.1.10",
  cell2="X4T1.1.11")
# Divide the counts into non-variable and variable count
VariableFeaturePlot(object=sce)
# Violin nad Ridge plots
VlnPlot(object=sce,features=c("Ccl5","Cd14","Cd19"))
# Export the graph to PPT
# graph2ppt(file="Violin plot.pptx",width=7,height=5)
RidgePlot(object=sce,feature=c("Ccl5","Cd14","Cd19"))
# Heatmaps
DimHeatmap(object=sce,reduction="pca",cells=200,balanced=TRUE)

#############################################
# Determine significant principal component #
#############################################
# Determine the statistically significant principal components
# Identify those significant PCs that have a strong enrichment of low p-values
sce<-JackStraw(object=sce,reduction="pca",dims=27,num.replicate=100,prop.freq=0.1,verbose=TRUE)
# Visualize the distribution of p-values in each significant PC 
sce<-ScoreJackStraw(object=sce,dims=1:27,reuction="pca")
JackStrawPlot(object=sce,dims=1:27,reduction="pca")
# Determine which PCs to use and draw the cutoff according to the standard deviations
ElbowPlot(object=sce)

###################
# Cell Clustering #
###################
# Calculate the k-nearest neighbors and construct the SNN graph
sce<-FindNeighbors(sce,reduction="pca",dims=1:20)
sce<-FindClusters(sce,resolution=0.5,algorithm=1)

############################################
# Perform non-linear dimensional reduction #
############################################
# Visualize and explore the datasets using tSNE
sce<-RunTSNE(object=sce,dim.use=1:27)
DimPlot(object=sce,reduction="tsne",do.label=TRUE)
# Get quick information from a scatterplot by hovering over points towards some features
plot<-DimPlot(object=sce,reduction="tsne",do.label=TRUE)
HoverLocator(plot=plot,information=FetchData(object=sce,vars="Gnai3"))

############
# Run UMAP #
############
# install python package "umap-learn"
sce<-RunUMAP(sce,reduction="pca",dims=1:21)
DimPlot(sce,reduction="umap",pt.size=1)
DimPlot(sce,reduction="umap",split.by="seurat_clusters")

#######################################
# Find Differentially Expressed genes #
#######################################
# Find Cluster1 markers
cluster1.markers<-FindMarkers(object=sce,ident.1=1,min.pct=0.25)
print(x=head(x=cluster1.markers))
# Find cluster5 markers distinguishing from clusters 0 and 3
cluster5.markers<-FindMarkers(object=sce,ident.1=5,ident.2=c(0,3),min.pct=0.25)
print(x=head(x=cluster5.markers,n=6))
# Find markers for every cluster compared to all other cells
sce.markers<-FindAllMarkers(object=sce,only.pos=TRUE,min.pct=0.25,thresh.use=0.25)
sce.markers %>% group_by(cluster) %>% top_n(3,avg_logFC)
cluster2.markers<-FindMarkers(object=sce,ident.1=2,thresh.use=0.25,test.use="roc",only.pos=TRUE)
VlnPlot(object=sce,features=c("Ighm","Igkc"))
FeaturePlot(object=sce,feature=c("Ighm","Igkc"),cols=c("grey","blue"),reduction="umap")
top10<-sce.markers %>% group_by(cluster) %>% top_n(10,avg_logFC)
DoHeatmap(object=sce,features=top10$gene,label=TRUE)

###############################################
# Assign the celltype identity to the cluster #
###############################################
current.cluster.ids<-c(0,1,2,3,4,5,6,7,8)
new.cluster.ids<-c("CD4 T Cells","CD8 T Cells","CD14 T Cells","CD4 B Cells","CD8 B Cells","CD14 B Cells","Tumour Cells","Macrophage","Vessel")
pbmc@active.ident<-plyr::mapvalues(x=pbmc@active.ident,from=current.cluster.ids,to=new.cluster.ids)
DimPlot(object=pbmc,reduction="tsne",do.label=TRUE,pt.size=0.5)

###############################
# Further subdivide celltypes #
###############################
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

#============================#
#       Musician: Resonance  #
#           Date: 2019/08/02 #
# Revised author: Resonance  #
#           Time: 2019/08/08 #
#============================#