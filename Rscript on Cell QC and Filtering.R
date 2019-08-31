# Load libraries
library(SingleCellExperiment)
library(scater)
library(scran)
options(stringsAsFactors=FALSE)

# Load data from the expression matrix and annotation
data<-read.table("m4T1.raw.featurecounts.txt")
anno<-read.csv("m4T1.cell.annotation.csv")
feature<-data[,1:4]
counts<-data[,5:ncol(data)]

# Create the SCE object
sce<-SingleCellExperiment(assays=list(counts=as.matrix(counts)),
                          colData=anno)
rowData(sce)$feature_symbol<-feature[,3]

# Remove genes that are not expressed in any cell
keep_feature<-rowSums(counts(sce)>0)>0
sce<-sce[keep_feature,]
dim(sce)

# Remove mcherry & GFP gene
remove_cell_mCherry<-counts(sce)["mCherry",]==0 & counts(sce)["luciGFP",]==0
table(remove_cell_mCherry)
colData(sce)$use<-!remove_cell_mCherry
table(colData(sce)$use)
remove_gene_mcherry<-rownames(sce) %in% c("mCherry","luciGFP")
rowData(sce)$use<-!remove_gene_mcherry
table(rowData(sce)$use)
sce<-sce[rowData(sce)$use,colData(sce)$use]
dim(sce)

# Define control features: ERCC & mitochondrial genes
isSpike(sce,'ERCC')<-grepl("^ERCC",rownames(sce))
isSpike(sce,"mt")<-rownames(sce) %in% c("mt-Atp6","mt-Atp8","mt-Co1","mt-Co2","mt-Co3",
    "mt-Cytb","mt-Nd1","mt-Nd2","mt-Nd3","mt-Nd4","mt-Nd4I","mt-Nd5","mt-Nd6","mt-Rnr1",
    "mt-Rnr2","mt-Ta","mt-Tc","mt-Te","mt-Tf","mt-Th","mt-Ti","mt-Tk","mt-TI1","mt-TI2",
    "mt-Tm","mt-Tn","mt-Tp","mt-Tq","mt-Tr","mt-Ts1","mt-Ts2","mt-Tt","mt-Tv","mt-Tw","mt-Ty")

# Calculate the quality metrics:
sce<-calculateQCMetrics(sce,
                        feature_controls=list(
                          ERCC=isSpike(sce,"ERCC"),
                          MT=isSpike(sce,"mt"))
                        )

# Cell QC
hist(sce$total_counts,breaks=100)
abline(v=25000,col="red")
hist(sce$total_features_by_counts,breaks = 100)
abline(v=7000,col="red")
plotColData(sce,x="total_features_by_counts",y="pct_counts_ERCC")
plotColData(sce,x="total_features_by_counts",y="pct_counts_MT")

# Cell filtering
# filter_by_expr_features>5000
filter_by_expr_features<-(sce$total_features_by_counts>5000)
table(filter_by_expr_features)
# filter_total_counts>30000
filter_by_total_counts<-(sce$total_counts>5e5 & sce$total_counts<1e7)
table(filter_by_total_counts)
# filter_by_ERCC
filter_by_ERCC<-(sce$pct_counts_ERCC<10)
table(filter_by_ERCC)
# filter_by_MT
library(dplyr)
library(Seurat)
pbmc<-CreateSeuratObject(counts=counts,project="pbmc3k",min.cells=3,min.features=200)
pbmc[["percent.mt"]]<-PercentageFeatureSet(object=pbmc,pattern="^mt-")
VlnPlot(object=pbmc,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
filter_by_MT<-(sce$pct_counts_MT<10)
table(filter_by_MT)
# filter_by_map
pct.map<-read.table("lane1-13.Alinge.static.txt",sep="\t",header=T)
head(pct.map)
colnames(pct.map)[1]<-c("CellName")
which(remove_cell_mCherry)
pct.map<-subset(pct.map,CellName!="4T1-1-35" & CellName!="4T1-3-73")
filter_by_map<-pct.map[,3]/pct.map[,2]*100>40
table(filter_by_map)

# Manual filtering
sce$use<-(filter_by_expr_features 
          & filter_by_total_counts 
          & filter_by_ERCC
          & filter_by_MT 
          & filter_by_map)
table(sce$use)
colnames(sce)[!sce$use]
filter_genes<-apply(
  counts(sce[,colData(sce)$use]),1,
  function(x) length(x[x>1])>=2
)
rowData(sce)$use<-filter_genes
table(rowData(sce)$use)
# Automatic filtering
sce<-runPCA(sce,use_coldata=TRUE,detect_outliers=TRUE)
reducedDimNames(sce)
table(sce$outlier)
plotReducedDim(sce,use_dimred="PCA_coldata",size_by="total_features_by_counts",
               shape_by="use",colour_by="outlier")
# Compare the two filterings
library(limma)
man=colnames(sce)[!sce$use]
auto=colnames(sce)[sce$outlier]
venn.diag<-vennCounts(cbind(colnames(sce) %in% man,colnames(sce) %in% auto))
vennDiagram(venn.diag,names=c("Manual","Automatic"),counts.col=c("red"),circle.col=c("blue","green"))

# Gene analysis
# Gene expression
plotHighestExprs(sce,exprs_values="counts")
# Gene filtering
keep_feature<-nexprs(sce[,colData(sce)$use],byrow=TRUE,detection_limit=1)>=2
rowData(sce)$use<-keep_feature
table(keep_feature)

# Save the data
# Save the corresponding data after Cell QC and filtering
sce.QC<-sce[rowData(sce)$use,colData(sce)$use]
# Create an additional slot with log-transformed counts and remove PCA results
assay(sce.QC,"logcounts_raw")<-log2(counts(sce.QC)+1)
reducedDim(sce)<-NULL
# Output the data
saveRDS(sce.QC,file="m4T1.sce.QC.rds")
