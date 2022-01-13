
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

    require(HGC)
    require(Seurat)

    sO = FindNeighbors(sO, dims = 1:noDim, reduction = method)

    sO = FindClusteringTree(sO, graph.type = "SNN")

    sO.tree <- sO@graphs$ClusteringTree

    hgc_clusters = cutree(sO.tree, k = k)

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

h_clustering = function(sO, linkM, noDim, reductionMethod, method, numClusters, hc_option) {

    results = reactiveValues()

    require(bluster)
    require(dynamicTreeCut)

    if (reductionMethod == "umap") {
        sO = RunUMAP(sO, dims = 1:noDim, reduction = method)
    } else {
        sO = RunTSNE(object = sO, dims = 1:noDim, reduction = method)
    }


    sce = as.SingleCellExperiment(sO)

    rowData(sce)$feature_symbol = rownames(sO)

    counts(sce) = as.matrix(counts(sce))

    logcounts(sce) = as.matrix(logcounts(sce))

    mat = reducedDim(sce, toupper(method))

    mat = mat[,1:noDim]  #only pass the top components on to the clustering algorithm

    if(hc_option == "Use determined Cluster Number") {

        hp = HclustParam(method= linkM, cut.params = list(k = numClusters))    #cut.dynamic=TRUE)

    } else {

        hp = HclustParam(method= linkM, cut.dynamic=TRUE)

    }


    sce.HClust <- clusterRows(mat, hp)

    names(sce.HClust) <- colnames(sO)
    sO = AddMetaData(sO, sce.HClust, col.name = "h_clusters")

    sO@active.ident <- as.factor(sce.HClust)

    hc_Plot = DimPlot(sO, reduction = reductionMethod, group.by = "h_clusters")

    results$data = sO
    results$plot = hc_Plot

    out = results
}

km_clustering = function(sO, k, noDim, reductionMethod, method) {

    results = reactiveValues()

    require(bluster)

    if (reductionMethod == "umap") {
        sO = RunUMAP(sO, dims = 1:noDim, reduction = method)
    } else {
        sO = RunTSNE(object = sO, dims = 1:noDim, reduction = method)
    }

    set.seed(100)

    sce = as.SingleCellExperiment(sO)

    rowData(sce)$feature_symbol = rownames(sO)

    counts(sce) = as.matrix(counts(sce))

    logcounts(sce) = as.matrix(logcounts(sce))

    mat = reducedDim(sce, toupper(method))

    mat = mat[,1:noDim] #only pass the top components on to the clustering algorithm

    sce.kmeans <- clusterRows(mat, KmeansParam(k))

    sO@meta.data$km_clusters = sce.kmeans

    sO@active.ident <- as.factor(sce.kmeans)

    km_Plot = DimPlot(sO, reduction = reductionMethod, group.by = "km_clusters")


    results$data = sO
    results$plot = km_Plot

    out = results


}


monocleClustering = function(sO, noK, noDim, reductionMethod, method) { #functional

    results = reactiveValues()

    require(VGAM)

    sO = RunTSNE(object = sO, dims = 1:noDim, reduction = method)

    cds = as.CellDataSet(sO, assay = "RNA", reduction = "tsne") # easily convert Seurat to CDS

    cds@expressionFamily = uninormal() # The correct expression family for normalized/scaled/etc... data

    #extract variable features from Seurat and set Ordering Filter in Monocle
    var_genes = sO[["RNA"]]@var.features
    cds = setOrderingFilter(cds, var_genes)

    cds = monocle::clusterCells(cds, num_clusters = noK + 1)

    #Monocle-Seurat Conversion
    Mclusters <- phenoData(cds)$Cluster
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


consensus_clustering = function(sO, k, noDim, reductionMethod, method) {

    results = reactiveValues()

    #Seurat to SingleCellExperiment

    sce = as.SingleCellExperiment(sO)

    rowData(sce)$feature_symbol = rownames(sO)

    counts(sce) = as.matrix(counts(sce))

    logcounts(sce) = as.matrix(logcounts(sce))

    #Feature Selection Filter

    var_genes = sO[["RNA"]]@var.features

    sce = sce[var_genes,]

    #Run SC3 and extract results

    sce = sc3(sce, ks = k, gene_filter = FALSE)

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

    require(scDHA)

    expr_mat = GetAssayData(object = sO, slot = "counts") #extract data

    expr_mat = t(expr_mat) #transpose the data

    expr_mat = log2(expr_mat + 1)

    expr_mat = as(expr_mat, "dgCMatrix")

    expr_results = scDHA(expr_mat, seed = 1, sparse = TRUE, k = kv - 1, ncores = core_num) #clustering step

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



clusterValidation = function(sO, plot, maxK) {

    require(NbClust)
    require(ggplot2)
    require(factoextra)

    all.genes = rownames(sO)

    sO = ScaleData(sO, features = all.genes)

    ScaledData = sO@assays$RNA@scale.data

    vf = sO@assays$RNA@var.features

    SmallData = ScaledData[vf,]

    if(plot == "Elbow Plot") {

        ElbowPlot = fviz_nbclust(SmallData, kmeans, method = "wss", k.max = maxK) + theme_minimal() + ggtitle("the Elbow Method")

        return(ElbowPlot)

    } else {

        SilhouettePlot = fviz_nbclust(SmallData, kmeans, method = "silhouette", k.max = maxK) + theme_minimal() + ggtitle("The Silhouette Plot")

        return(SilhouettePlot)
    }

}


sc3_estimate = function(sO) {

    require(Seurat)
    require(SingleCellExperiment)

    sce = as.SingleCellExperiment(sO)

    sce = sc3_estimate_k(sce)

    clusNum = sce@metadata$sc3$k_estimation

    out = clusNum

}
