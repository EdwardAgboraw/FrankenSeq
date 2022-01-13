
#Package Installation:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")


#General Framework
install.packages("shiny")
install.packages("Seurat")
install.packages("shinythemes")
install.packages("shinycssloaders")
install.packages("magrittr")

BiocManager::install("SingleCellExperiment")

#For Feature Selection
install.packages("DUBStepR")

BiocManager::install("M3Drop")

#For Dimension Reduction
install.packages("VGAM")

remotes::install_github('satijalab/seurat-wrappers')

BiocManager::install("scry")

#For Cluster Validation
install.packages("NbClust")
install.packages("ggplot2")
install.packages("factoextra")

#For Cluster Analysis
BiocManager::install("HGC")
BiocManager::install("monocle")
BiocManager::install("SC3")

install.packages("scDHA")

BiocManager::install("bluster")
install.packages("dynamicTreeCut")
