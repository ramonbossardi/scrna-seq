library(Seurat)
library(dplyr)
library(Matrix)
library(readxl)
library(cowplot)
library(patchwork) #combine plots
library(hdf5r)
library(SingleR)
library(celldex)
library(RColorBrewer)
library(ggplot2)
library(SingleCellExperiment)

setwd("~/Analises R/single cell/Singer")

save.image(file = "whole_aorta.Rdata")
load("whole_aorta.Rdata")

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117963

#Download files .h5 from GSE - Copy the link ftp
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117963/suppl/GSE117963_10X_whole_aorta_filtered_gene_bc_matrices_h5.h5",
              destfile = "GSE117963_10X_whole_aorta")

#Transform in different matrix
matrix_whole_aorta <- Read10X_h5("~/Analises R/single cell/Singer/GSE117963_10X_whole_aorta",use.names = T)
srat_whole_aorta <- CreateSeuratObject(matrix_whole_aorta,project = "whole_aorta")

#remove matrices to save memory:
rm(monaco.fine, monaco.ref, monaco.main, sce, hpca.fine, hpca.main, hpca.ref, dice.fine, dice.main, dice.ref, CIPR_all_results, CIPR_top_results)

#number of UMI reads detected per cell (nCount_RNA), and the number of expressed (detected) genes per same cell (nFeature_RNA).
VlnPlot(srat_whole_aorta, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Conventional way is to scale it to 10,000 (as if all cells have 10k UMIs overall), and log2-transform the obtained values.
FeatureScatter(srat_whole_aorta,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

srat_whole_aorta<- NormalizeData(srat_whole_aorta)
#discovers the most variable features (genes) - these are usually most interesting for downstream analysis.
srat_whole_aorta <- FindVariableFeatures(srat_whole_aorta, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes:
top10 <- head(VariableFeatures(srat_whole_aorta), 10)
top10 

#Plot variable features with and without labels:
  
plot1 <- VariableFeaturePlot(srat_whole_aorta)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

#ScaleData converts normalized gene expression to Z-score (values centered at 0 and with variance of 1).
all.genes <- rownames(srat_whole_aorta)
srat_whole_aorta <- ScaleData(srat_whole_aorta, features = all.genes)

#Do PCA, which is a common way of linear dimensionality reduction. By default we use 2000 most variable genes.
srat_whole_aorta <- RunPCA(srat_whole_aorta, features = VariableFeatures(object = srat_whole_aorta))

# visualizing both cells and features that define the PCA
print(srat_whole_aorta[["pca"]], dims = 1:10, nfeatures = 20)

VizDimLoadings(srat_whole_aorta, dims = 1:9, reduction = "pca") & 
  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))

ElbowPlot(srat_whole_aorta)

#It’s often good to find how many PCs can be used without much information loss.
DimPlot(srat_whole_aorta, reduction = "pca")

#do clustering. Higher resolution leads to more clusters (default is 0.8). 
#It would be very important to find the correct cluster resolution in the future, since cell 
#type markers depends on cluster definition.
srat_whole_aorta <- FindNeighbors(srat_whole_aorta, dims = 1:10)

srat_whole_aorta <- FindClusters(srat_whole_aorta, resolution = 0.8)

#visualization purposes, we also need to generate UMAP reduced dimensionality representation
srat_whole_aorta <- RunUMAP(srat_whole_aorta, dims = 1:10, verbose = F)

#cluster sizes
table(srat_whole_aorta@meta.data$seurat_clusters)

#DimPlot uses UMAP by default, with Seurat clusters as identity:
DimPlot(srat_whole_aorta,label.size = 4,repel = T,label = T)

#In order to control for clustering resolution and other possible artifacts, we will take a close 
#look at two minor cell populations: 1) dendritic cells (DCs), 2) platelets, aka thrombocytes. Let’s 
#visualise two markers for each of this cell type: LILRA4 and TPM2 for DCs, and PPBP and GP1BB for platelets.

FeaturePlot(srat_whole_aorta, features = c("Camk2d", "Camk2g", "Myh11", "Acta2"))
FeaturePlot(srat_whole_aorta, features = c("Tagln", "Itga8", "Col1a2", "Ly6a", "Gli1"))

#visualize other confounders: # needo to have two dataset
#FeaturePlot(srat_whole_aorta, features = "Doublet_score") & theme(plot.title = element_text(size=10))

#cell cycle scores, as described here. This has to be done after normalization and scaling. 
#Seurat has a built-in list, cc.genes.updated.2019 (newer), that defines genes involved in cell cycle.

cc.genes.updated.2019

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

#https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html
srat_whole_aorta <- CellCycleScoring(srat_whole_aorta, s.features = s.genes, g2m.features = g2m.genes)
table(srat_whole_aorta[[]]$Phase)

FeaturePlot(srat_whole_aorta,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))

#8.4 Differential expression and marker selection

DimPlot(srat_whole_aorta, label = T)

#First, let’s set the active assay back to “RNA,” and re-do the normalization and scaling 
#(since we removed a notable fraction of cells that failed QC):

DefaultAssay(srat_whole_aorta) <- "RNA"
srat_whole_aorta <- NormalizeData(srat_whole_aorta)
srat_whole_aorta <- FindVariableFeatures(srat_whole_aorta, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srat_whole_aorta)
srat_whole_aorta <- ScaleData(srat_whole_aorta, features = all.genes)

#The following function allows to find markers for every cluster by comparing it to all 
#remaining cells, while reporting only the positive ones.
all.markers <- FindAllMarkers(srat_whole_aorta, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)

dim(all.markers)
table(all.markers$cluster)

#Look 3 markers in each cluster
top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers

#Cell type annotation using SingleR
#automaic annotation with SingleR is workflow-agnostic (can be used with Seurat, SCE, etc)
monaco.ref <- celldex::MonacoImmuneData()
# hpca.ref <- celldex::HumanPrimaryCellAtlasData()
# dice.ref <- celldex::DatabaseImmuneCellExpressionData()

#Let’s convert our Seurat object to single cell experiment (SCE) for convenience. 
#After this, using SingleR becomes very easy:

sce <- as.SingleCellExperiment(DietSeurat(srat_whole_aorta))
sce

monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
# hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
# hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
# dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
# dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)

#summary of general cell type annotations.
table(monaco.main$pruned.labels)

#table(hpca.main$pruned.labels)
#table(dice.main$pruned.labels)

#The finer cell types annotations are you after, the harder they are to get reliably. 
#This is where comparing many databases, as well as using individual markers from literature, would all be very valuable.

table(monaco.fine$pruned.labels)

# table(hpca.fine$pruned.labels)
# table(dice.fine$pruned.labels)

#add the annotations to the Seurat object metadata 

srat_whole_aorta@meta.data$monaco.main <- monaco.main$pruned.labels
srat_whole_aorta@meta.data$monaco.fine <- monaco.fine$pruned.labels

# srat_whole_aorta@meta.data$hpca.main   <- hpca.main$pruned.labels
# srat_whole_aorta@meta.data$dice.main   <- dice.main$pruned.labels
# srat_whole_aorta@meta.data$hpca.fine   <- hpca.fine$pruned.labels
# srat_whole_aorta@meta.data$dice.fine   <- dice.fine$pruned.labels

srat_whole_aorta <- SetIdent(srat_whole_aorta, value = "monaco.fine")
DimPlot(srat_whole_aorta, label = T , repel = T, label.size = 3) + NoLegend()


## starts<-rep(100,40)
## fx<-function(nstart) kmeans(Boston, 4, nstart=nstart)
## numCores<-detectCores()
## numCores
## system.time(results<-lapply(starts, fx))
## system.time(results<-mclapply(starts, fx, mc.cores = numCores))

# Integration dataset

#Download files .h5 from GSE - Copy the link ftp

download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117963/suppl/GSE117963_10X_whole_aorta_filtered_gene_bc_matrices_h5.h5",
              destfile = "GSE117963_10X_whole_aorta")

setwd("~/Analises R/single cell/Singer/positive_whole")

save.image(file = "whole_aorta_medial_layer.Rdata")
load("whole_aorta_medial_layer.Rdata")

#Transform in different matrix

download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117963/suppl/GSE117963_10X_whole_aorta_filtered_gene_bc_matrices_h5.h5",
              destfile = "GSE117963_10X_whole_aorta")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117963/suppl/GSE117963_10X_lineage_positive_filtered_gene_bc_matrices_h5.h5",
              destfile = "GSE117963_10X_lineage_positive")

whole_aorta.data <- Read10X_h5("~/Analises R/single cell/Singer/positive_whole/GSE117963_10X_whole_aorta",use.names = T)
lineage_positive.data <- Read10X_h5("~/Analises R/single cell/Singer/positive_whole/GSE117963_10X_lineage_positive",use.names = T)

rm(whole_aorta.data)
rm(lineage_positive.data)


# Initialize the Seurat object with the raw (non-normalized data).
s1 <- CreateSeuratObject(counts = whole_aorta.data, project = "whole_aorta", min.cells = 3, min.features = 200)
s2 <- CreateSeuratObject(counts = lineage_positive.data, project = "Medial layer", min.cells = 3, min.features = 200)
s1$stim <- 's1'
s2$stim <- 's2'
aorta <-merge(s1,y=s2, add.cell.ids = c("s1","s2"), project="aorta")
aorta

head(aorta)
levels(aorta)










