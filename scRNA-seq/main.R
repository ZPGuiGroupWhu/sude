## --Arguments--
##     'filename': filename of the scRNA-seq dataset
##     'labels': whether the file contains the true annotations of cells ('1' is Yes, '0' is No)

filename = 'NdpKo'
labels = '1'

## Generate Seurat Object from the raw data
source('BuildSeuratObject.R')
SeuratData <- BuildSeuratObject(filename, format, labels)

## Preprocess scRNA-seq using standard Seurat pipeline
source('SeuratPreprocess.R')
SeuratData <- SeuratPreprocess(SeuratData, filename)

## Plot the clustering result
Idents(SeuratData)<-SeuratData@meta.data[["Cluster"]]
DimPlot(SeuratData, pt.size=1) + NoLegend()
