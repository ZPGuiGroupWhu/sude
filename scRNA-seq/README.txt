---Attention---
Before running the code, please decompress the two .zip files.


---Depends---
## Please click Tools->Global Options->Packages, change CRAN repository to a near mirror. Then, execute the following code:
## Install packages from CRAN.
install.packages(c("argparse", "assertthat", "ClusterR", "dotCall64", "fields", "geometry", "ggplot2", "gtools", "mclust", "prodlim", "RcppHungarian", "readr", "reshape2", "Seurat", "spam", "uwot"))
## Determine whether the package "BiocManager" exists, if not, install this package.
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
## Install packages from Bioconductor.
BiocManager::install(c("BiocGenerics", "BiocSingular", "GenomeInfoDb", "GenomicRanges", "IRanges", "MatrixGenerics", "S4Vectors", "scran", "scuttle", "SingleCellExperiment", "SummarizedExperiment"), force = TRUE, update = TRUE, ask = FALSE)