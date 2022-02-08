
featureSelection = function(sO, fs_method, number_of_features) {

    require(Seurat)
    require(M3Drop)
    require(bluster)
    require(DUBStepR)

    results = reactiveValues()

    if (fs_method == "HVG - Seurat (vst)") {

        sO = FindVariableFeatures(sO, selection.method = "vst", nfeatures = number_of_features)

        feature_genes = sO@assays$RNA@var.features

        fsplot = VariableFeaturePlot(sO)

        top10 = head(VariableFeatures(sO), 10)

        fs_plot = LabelPoints(plot = fsplot, points = top10, repel = TRUE)

        results$data = sO

        results$plot = fs_plot

        results$featurelist = feature_genes


    }  else if (fs_method ==   "HVG - Seurat (mvp)" ) {

        sO = FindVariableFeatures(sO, selection.method = "mvp")

        feature_genes = sO@assays$RNA@var.features

        fsplot = VariableFeaturePlot(sO)

        top10 = head(VariableFeatures(sO), 10)

        fs_plot = LabelPoints(plot = fsplot, points = top10, repel = TRUE)

        results$data = sO

        results$plot = fs_plot

        results$featurelist = feature_genes

    } else if (fs_method ==   "HVG - Seurat (Dispersion)" ) {

        sO = FindVariableFeatures(sO, selection.method = "disp", nfeatures = number_of_features)

        feature_genes = sO@assays$RNA@var.features

        fsplot = VariableFeaturePlot(sO)

        top10 = head(VariableFeatures(sO), 10)

        fs_plot = LabelPoints(plot = fsplot, points = top10, repel = TRUE)

        results$data = sO

        results$plot = fs_plot

        results$featurelist = feature_genes

    } else if (fs_method ==   "Drop-Out Based Feature Selection - M3Drop" ) {

        counts = sO@assays$RNA@counts

        count_mat = NBumiConvertData(counts, is.counts=TRUE)

        count_mat = NBumiConvertData(count_mat, is.counts=TRUE)

        DANB_fit = NBumiFitModel(count_mat)

        NBDropFS = NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", ntop = number_of_features, suppress.plot=TRUE)

        feature_genes = row.names(NBDropFS)

        sO@assays$RNA@var.features = feature_genes

        top10 = head(feature_genes, 10)

        fsplot = VariableFeaturePlot(sO)

        fs_plot = LabelPoints(plot = fsplot, points = top10, repel = TRUE)

        results$data = sO

        results$plot = fs_plot

        results$featurelist = feature_genes


    } else if (fs_method == "Gene-Gene Correlation Feature Selection (DubStepR)") {

        dubstepR.out = DUBStepR(input.data = sO@assays$RNA@data, min.cells = 0.05*ncol(sO), optimise.features = TRUE, k = 10, num.pcs = 20, error = 0)

        feature_genes = dubstepR.out$optimal.feature.genes

        top10 = head(feature_genes, 10)

        sO@assays$RNA@var.features = feature_genes

        fsplot = VariableFeaturePlot(sO)

        fs_plot = LabelPoints(plot = fsplot, points = top10, repel = TRUE)

        results$data = sO

        results$plot = fs_plot

        results$featurelist = feature_genes

    } else if (fs_method == "Filter by Gene Expression") {

        sce = as.SingleCellExperiment(sO)

        LogCounts = as.array(logcounts(sce))

        exprsn = rowMeans(LogCounts)

        keep = order(exprsn, decreasing = TRUE)[seq_len(number_of_features)]

        sce[keep, ]

        feature_genes = rownames(sce)

        sO@assays$RNA@var.features = feature_genes

        top10 = head(feature_genes, 10)

        fsplot = VariableFeaturePlot(sO)

        fs_plot = LabelPoints(plot = fsplot, points = top10, repel = TRUE)

        #outputs
        results$data = sO

        results$plot = fs_plot

        results$featurelist = feature_genes

    } else if (fs_method == "Feature Selection by Deviance (Scry)") {

        require(scry)

        m = GetAssayData(sO, slot = "counts", assay = "RNA")

        #select features

        deviantFeatures = scry::devianceFeatureSelection(m)

        dev_ranked_genes = rownames(sO)[order(deviantFeatures, decreasing = TRUE)]

        topdev = head(dev_ranked_genes, number_of_features)

        top10 = head(topdev, 10)

        sO@assays$RNA@var.features = topdev

        #create plot
        fsplot = VariableFeaturePlot(sO)

        fs_plot = LabelPoints(plot = fsplot, points = top10, repel = TRUE)

        #outputs

        results$data = sO

        results$plot = fs_plot

        results$featurelist = topdev


    }

    out = results

}

#The output of this function is always a seurat object that contains the processed data, and a list of selected genes stored in var.features.

