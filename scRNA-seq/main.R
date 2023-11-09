## Specify the filename, format and labels.
## --Arguments--
##     filename: filename of the scRNA-seq dataset
##     format: data format of the scRNA-seq dataset (support two data formats, csv and 10X)
##     labels: whether the file contains the true label of cells ('1' is Yes, '0' is No)
filename = 'pbmc3k'
format = '10X'
labels = '1'

## Generate Seurat Object from the raw data
source('BuildSeuratObject.R')
SeuratData <- BuildSeuratObject(filename, format, labels)

## Preprocess scRNA-seq using standard Seurat pipeline
## --Arguments--
##    UMAP_Dim: The dimension of UMAP space to embed into. (Default: 2, Recommended: 2~5)
source('SeuratPreprocess.R')
UMAP_Dim = 5
SeuratData <- SeuratPreprocess(SeuratData, filename, UMAP_Dim)

## Cluster the cells using CDC algorithm
## --Arguments--
##     k: k of KNN (Default: 30, Recommended: 30~50) 
##     ratio: percentile ratio of internal points (Default: 0.9, Recommended: 0.75~0.95, 0.55~0.65 for pbmc3k)
source('CDC.R')
k = 40
ratio = 0.62
Idents(SeuratData) <- CDC(SeuratData@reductions[["umap"]]@cell.embeddings, k, ratio)

## Plot the clustering result
DimPlot(SeuratData, pt.size=1) + NoLegend()

## Evaluate the clustering accuracy using ARI metric
ARI <- mclust::adjustedRandIndex(Idents(SeuratData), SeuratData@meta.data[["Cluster"]])



