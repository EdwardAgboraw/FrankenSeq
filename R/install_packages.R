
#Basic Framework Installation Package

#' basic framework
#'
#' @return
#' @export
#'
#' @examples
basic_framework = function() {

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")


    if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")

    #Background:

    install.packages("dplyr")
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

    #For Cluster Analysis
    BiocManager::install("HGC")
    BiocManager::install("monocle")
    BiocManager::install("SC3")

    BiocManager::install("bluster")
    install.packages("dynamicTreeCut")

}


#' cluster_validation
#'
#' @return
#' @export
#'
#' @examples
cluster_validation = function() {

  #For Cluster Validation
  install.packages("NbClust")
  install.packages("ggplot2")
  install.packages("factoextra")

}


#' deep_learning
#'
#' @return
#' @export
#'
#' @examples
deep_learning = function() {

  install.packages("scDHA")

}




