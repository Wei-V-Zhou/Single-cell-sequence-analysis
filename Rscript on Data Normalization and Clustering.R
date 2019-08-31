# Data Normalization
# Load libraries
library(SingleCellExperiment)
library(scRNA.seq.funcs)
library(scater)
library(scran)
options(stringsAsFactors=FALSE)
set.seed(1234567)

# Load data from the cell QC & filtering results
sce.QC<-readRDS("m4T1.sce.QC.rds")
endog_genes<-!rowData(sce.QC)$is_feature_control

# Three normalization methods for scRNA-seq
# CPM normalization
logcounts(sce.QC)<-log2(calculateCPM(sce.QC,use_size_factors=FALSE)+1)
tmp<-runPCA(sce.QC[endog_genes, ],exprs_values="logcounts")
plotPCA(tmp,colour_by="CellType",size_by="total_features_by_counts",shape_by="use")
# # scran normalization
# qclust<-quickCluster(sce.QC,min.size=30)
# sce.QC<-computeSumFactors(sce.QC,sizes=15,clusters=qclust)
# sce.QC<-normalize(sce.QC)
# plotPCA(sce.QC[endog_genes, ],colour_by="CellType",size_by="total_features_by_counts",shape_by="use")
# # Downsampling normalization
# logcounts(sce.QC)<-log2(Down_Sample_Matrix(counts(sce.QC))+1)
# plotPCA(sce.QC[endog_genes, ],colour_by="CellType",size_by="total_features_by_counts",shape_by="use")
assay(sce.QC,"logcounts")<-logcounts(sce.QC)
assay(sce.QC,"logcounts_raw")<-NULL
reducedDim(sce.QC)<-NULL
saveRDS(sce.QC,file="m4T1.sce.norm.rds")

# Clustering Analysis
# Load libraries
library(SingleCellExperiment)
library(pcaMethods)
library(SC3)
library(scater)
library(pheatmap)
library(mclust)
set.seed(1234567)

# Load data from the cell normalization results
sce.norm<-readRDS("m4T1.sce.norm.rds")
table(colData(sce.norm)$CellType)
plotPCA(sce.norm,colour_by="CellType")

# Three clustering methods for scRNA-seq
# SC clustering
sce.norm<-sc3(sce.norm,ks=10,biology=TRUE,n_cores=1)
sc3_plot_markers(sce.norm,k=10,show_pdata="CellType")
sc3_plot_expression(sce.norm,k=10,show_pdata="CellType")
sc3_plot_consensus(sce.norm,k=10,show_pdata="CellType")
# tSNE+k-means clustering
sce.norm<-runTSNE(sce.norm,rand_seed=1)
colData(sce.norm)$tSNE_kmeans<-as.character(kmeans(sce.norm@reducedDims$TSNE,centers=7)$clust)
plotTSNE(sce.norm,colour_by="tSNE_kmeans",shape_by="CellType")
# SINCERA/hierarchical clustering
sce.norm<-sc3(sce.norm,ks=10,biology=TRUE,n_cores=1)
input<-logcounts(sce.norm[rowData(sce.norm)$sc3_gene_filter, ])
dat<-apply(input,1,function(y) scRNA.seq.funcs::z.transform.helper(y))
dd<-as.dist((1-cor(t(dat),method="pearson"))/2)
hc<-hclust(dd,method="average")
num.singleton<-0
kk<-1
for(i in 2:dim(dat)[2]){
  clusters<-cutree(hc,k=i)
  clustersizes<-as.data.frame(table(clusters))
  singleton.clusters<-which(clustersizes$Freq<2)
  if (length(singleton.clusters)<=num.singleton){kk<-i}
  else{break;}
}
cat(kk)
pheatmap(t(dat),cluster_cols=hc,cutree_cols=kk,kmeans_k=100,show_rownames=FALSE)

