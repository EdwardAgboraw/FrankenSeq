#' Seurat (SNN Graph Based) Clustering Function
#'
#' @param sO the underlying Seurat object
#' @param noDim the number of dimensions to be use in the clustering step
#' @param res the Seurat resolution parameter
#' @param reductionMethod the type of graph (UMAP or T-SNE) used to represent the results
#' @param method the technique used to compress the data (PCA or GLM-PCA)
#'
#'@import Seurat
#'
#'@importFrom shiny reactiveValues
#'
#'
#'
#' @export
#'
#'
seuratClustering = function(sO, noDim, res, reductionMethod, method) { #functional

    results = reactiveValues()

    sO = FindNeighbors(sO, dims = 1:noDim, reduction = method)#, reduction = "DR")

    sO = FindClusters(sO, resolution = res, reduction = method)

    if (reductionMethod == "umap") {
        sO = RunUMAP(sO, dims = 1:noDim, reduction = method)
        clusterPlot = DimPlot(sO, reduction = "umap")

        results$plot = clusterPlot

    } else if (reductionMethod == "tsne") {
        sO = RunTSNE(object = sO, dims =1:noDim, reduction = method)
        clusterPlot = DimPlot(sO, reduction = "tsne")

        results$plot = clusterPlot

    } else {

        clusterPlot = DimPlot(sO, reduction = method)
        results$plot = clusterPlot

    }

    results$data = sO

    out = results

}


HGC_clustering = function(sO, k, noDim, reductionMethod, method) { #functional

    results = reactiveValues()

    sO = FindNeighbors(sO, dims = 1:noDim, reduction = method)

    sO = HGC::FindClusteringTree(sO, graph.type = "SNN")

    sO.tree <- sO@graphs$ClusteringTree

    hgc_clusters = dendextend::cutree(sO.tree, k = k)

    hgc_clusters = as.factor(hgc_clusters)

    sO@meta.data$hgc_clusters = hgc_clusters

    sO@active.ident <- as.factor(hgc_clusters)

    if (reductionMethod == "umap") {
        sO = RunUMAP(sO, dims = 1:noDim, reduction = method)

        hgc_plot = DimPlot(sO, reduction = reductionMethod, group.by = "hgc_clusters")

        results$plot = hgc_plot

    } else if (reductionMethod == "tsne") {

        sO = RunTSNE(object = sO, dims = 1:noDim, reduction = method)

        hgc_plot = DimPlot(sO, reduction = reductionMethod, group.by = "hgc_clusters")

        results$plot = hgc_plot

    } else if (reductionMethod == "pca") {

        hgc_plot = DimPlot(sO, reduction = reductionMethod, group.by = "hgc_clusters")

        results$plot = hgc_plot

    }

    results$data = sO


    out = results

}

#'@import SingleCellExperiment
#'@import dynamicTreeCut
h_clustering = function(sO, linkM, noDim, reductionMethod, method, numClusters, hc_option) {

    results = reactiveValues()

    if (reductionMethod == "umap") {
        sO = RunUMAP(sO, dims = 1:noDim, reduction = method)
    } else {
        sO = RunTSNE(object = sO, dims = 1:noDim, reduction = method)
    }


    sce = as.SingleCellExperiment(sO)

    SummarizedExperiment::rowData(sce)$feature_symbol = rownames(sO)

    counts(sce) = as.matrix(counts(sce))

    logcounts(sce) = as.matrix(logcounts(sce))

    mat = reducedDim(sce, toupper(method))

    mat = mat[,1:noDim]  #only pass the top components on to the clustering algorithm

    if(hc_option == "Use determined Cluster Number") {

        hp = bluster::HclustParam(method= linkM, cut.params = list(k = numClusters))    #cut.dynamic=TRUE)

    } else {

        hp = bluster::HclustParam(method= linkM, cut.dynamic=TRUE)

    }


    sce.HClust <- bluster::clusterRows(mat, hp)

    names(sce.HClust) <- colnames(sO)
    sO = AddMetaData(sO, sce.HClust, col.name = "h_clusters")

    sO@active.ident <- as.factor(sce.HClust)

    hc_Plot = DimPlot(sO, reduction = reductionMethod, group.by = "h_clusters")

    results$data = sO
    results$plot = hc_Plot

    out = results
}

#'
#'@import SingleCellExperiment
km_clustering = function(sO, k, noDim, reductionMethod, method) {

    results = reactiveValues()

    #require(bluster)

    if (reductionMethod == "umap") {
        sO = RunUMAP(sO, dims = 1:noDim, reduction = method)
    } else {
        sO = RunTSNE(object = sO, dims = 1:noDim, reduction = method)
    }

    set.seed(100)

    sce = as.SingleCellExperiment(sO)

    SummarizedExperiment::rowData(sce)$feature_symbol = rownames(sO)

    counts(sce) = as.matrix(counts(sce))

    logcounts(sce) = as.matrix(logcounts(sce))

    mat = reducedDim(sce, toupper(method))

    mat = mat[,1:noDim] #only pass the top components on to the clustering algorithm

    sce.kmeans <- bluster::clusterRows(mat, bluster::KmeansParam(k))

    sO@meta.data$km_clusters = sce.kmeans

    sO@active.ident <- as.factor(sce.kmeans)

    km_Plot = DimPlot(sO, reduction = reductionMethod, group.by = "km_clusters")


    results$data = sO
    results$plot = km_Plot

    out = results


}


monocleClustering = function(sO, noK, noDim, reductionMethod, method) { #functional

    results = reactiveValues()

    sO = RunTSNE(object = sO, dims = 1:noDim, reduction = method)

    cds = as.CellDataSet(sO, assay = "RNA", reduction = "tsne") # easily convert Seurat to CDS

    cds@expressionFamily = VGAM::uninormal() # The correct expression family for normalized/scaled/etc... data

    #extract variable features from Seurat and set Ordering Filter in Monocle
    var_genes = sO[["RNA"]]@var.features
    cds = monocle::setOrderingFilter(cds, var_genes)

    cds = monocle::clusterCells(cds, num_clusters = noK + 1)

    #Monocle-Seurat Conversion
    Mclusters <- Biobase::phenoData(cds)$Cluster
    names(Mclusters) <- colnames(sO)
    sO = AddMetaData(sO, Mclusters, col.name = "MClusters")
    sO@active.ident <- as.factor(Mclusters)

    if (reductionMethod == "umap") {
        sO = RunUMAP(sO, dims = 1:noDim, reduction = method)
    }

    #Mplot = plot_cell_clusters(cds)
    Mplot = DimPlot(sO, reduction = reductionMethod, group.by = "MClusters")

    results$data = sO
    results$plot = Mplot

    out = results

}


#'
#'@import SingleCellExperiment
consensus_clustering = function(sO, k, noDim, reductionMethod, method) {

    results = reactiveValues()

    #Seurat to SingleCellExperiment

    sce = as.SingleCellExperiment(sO)

   SummarizedExperiment::rowData(sce)$feature_symbol = rownames(sO)

    counts(sce) = as.matrix(counts(sce))

    logcounts(sce) = as.matrix(logcounts(sce))

    #Feature Selection Filter

    var_genes = sO[["RNA"]]@var.features

    sce = sce[var_genes,]

    #Run SC3 and extract results

    sce = SC3::sc3(sce, ks = k, gene_filter = FALSE)

    results_table = colData(sce)

    clusters = results_table[ , ncol(results_table)]

    sO@meta.data$sc3_clusters = clusters

    sO@active.ident = clusters

    Idents(sO) <- "sc3_clusters"

    #Create Output Graph

    if (reductionMethod == "umap") {
        sO = RunUMAP(sO, dims = 1:noDim, reduction = method)
    } else {
        sO = RunTSNE(object = sO, dims = 1:noDim, reduction = method)
    }

    cc_Plot = DimPlot(sO, reduction = reductionMethod, group.by = "sc3_clusters")


    results$data = sO
    results$plot = cc_Plot

    out = results


}

autoEncoderClusterring = function(sO, noDim, kv, core_num, reductionMethod, method) {

    results = reactiveValues()

    #require(scDHA)

    expr_mat = GetAssayData(object = sO, slot = "counts") #extract data

    expr_mat = as.matrix(expr_mat)

    expr_mat = t(expr_mat) #transpose the data

    expr_mat = log2(expr_mat + 1)

    expr_mat = methods::as(expr_mat, "dgCMatrix")

    expr_results = scDHA::scDHA(expr_mat, seed = 1, sparse = TRUE, k = kv, ncores = core_num) #clustering step

    clusters = expr_results$cluster

    sO@meta.data$HA_Clusters = clusters

    sO@active.ident <- as.factor(clusters)

    Idents(sO) <- "HA_Clusters"

    #Create Output Graph

    if (reductionMethod == "umap") {
        sO = RunUMAP(sO, dims = 1:noDim, reduction = method)
    } else {
        sO = RunTSNE(object = sO, dims = 1:noDim, reduction = method)
    }

    ha_Plot = DimPlot(sO, reduction = reductionMethod, group.by = "HA_Clusters")

    results$data = sO
    results$plot = ha_Plot

    out = results


}

#'
#' @import ggplot2
clusterValidation = function(sO, plot, maxK) {

    all.genes = rownames(sO)

    sO = ScaleData(sO, features = all.genes)

    ScaledData = sO@assays$RNA@scale.data

    vf = sO@assays$RNA@var.features

    SmallData = ScaledData[vf,]

    if(plot == "Elbow Plot") {

        ElbowPlot = factoextra::fviz_nbclust(SmallData, kmeans, method = "wss", k.max = maxK) + theme_minimal() + ggtitle("the Elbow Method")

        return(ElbowPlot)

    } else {

        SilhouettePlot = factoextra::fviz_nbclust(SmallData, kmeans, method = "silhouette", k.max = maxK) + theme_minimal() + ggtitle("The Silhouette Plot")

        return(SilhouettePlot)
    }

}


#'
#'@import SingleCellExperiment
sc3_estimate = function(sO) {

    sce = as.SingleCellExperiment(sO)

    sce = SC3::sc3_estimate_k(sce)

    clusNum = sce@metadata$sc3$k_estimation

    out = clusNum

}

generateTable = function(sO, tableType) {

    if (tableType == "Summary Report") {

        summaryTable <- table(Idents(sO))

        summaryTable <- as.data.frame(summaryTable)

        colnames(summaryTable) <- c("Cluster", "Frequency")

        #resultsTable(dataSummary)

        return(summaryTable)

    } else {

        fullTable <- as.data.frame(sO@active.ident)

        colnames(fullTable) <- "Cluster"

        fullTable <- cbind(rownames(fullTable), data.frame(fullTable, row.names = NULL))

        colnames(fullTable) <- c("Cell Label", "Cluster")

        #resultsTable(fullTable)

        return(fullTable)

    }

}

degAnalysis = function(sO, logFoldChange, minPercentage, outputOption) {

    sO.markers <- FindAllMarkers(sO, only.pos = TRUE, min.pct = minPercentage, logfc.threshold = logFoldChange)


    if (outputOption == "Heatmap") {

        sO.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) -> top10

        degHeatMap <- DoHeatmap(sO, features = top10$gene) + NoLegend()

        return(degHeatMap)

    } else {

        degTable <- sO.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

        return(degTable)
    }

}





