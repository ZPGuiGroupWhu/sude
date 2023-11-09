library(readr)
library(Seurat)

## This code aims to read the scRNA-seq data and convert it to a Seurat object. 
## The input data is required to be organized as the following standard csv or 10X formats.

BuildSeuratObject <- function (filename, datatype, labels){
  data <- Read10X(data.dir = filename)
  seuratSCE <- CreateSeuratObject(data, project=filename)
  if(labels=='1'){
    sample_info_df <- read.csv(paste(filename,'/Labels.csv',sep = ""))
    seuratSCE$Cluster <- sample_info_df$x
  }
  return(seuratSCE)
}

