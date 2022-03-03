
featureSelection = function(sO, fsMethod, featureNumber) {

    results = reactiveValues()

    if (fsMethod == "HVG - Seurat (vst)") {

        sO = FindVariableFeatures(sO, selection.method = "vst", nfeatures = featureNumber)

        feature_genes = sO@assays$RNA@var.features

        vfPlot = VariableFeaturePlot(sO)

        top10 = head(VariableFeatures(sO), 10)

        fsPlot = LabelPoints(plot = vfPlot, points = top10, repel = TRUE)

        #output

        results$data = sO

        results$plot = fsPlot

        results$featurelist = feature_genes


    }  else if (fsMethod ==   "HVG - Seurat (mvp)" ) {

        sO = FindVariableFeatures(sO, selection.method = "mvp")

        feature_genes = sO@assays$RNA@var.features

        vfPlot = VariableFeaturePlot(sO)

        top10 = head(VariableFeatures(sO), 10)

        fsPlot = LabelPoints(plot = vfPlot, points = top10, repel = TRUE)

        #output

        results$data = sO

        results$plot = fsPlot

        results$featurelist = feature_genes

    } else if (fsMethod ==   "HVG - Seurat (Dispersion)" ) {

        sO = FindVariableFeatures(sO, selection.method = "disp", nfeatures = featureNumber)

        feature_genes = sO@assays$RNA@var.features

        vfPlot = VariableFeaturePlot(sO)

        top10 = head(VariableFeatures(sO), 10)

        fsPlot = LabelPoints(plot = vfPlot, points = top10, repel = TRUE)

        #outputs

        results$data = sO

        results$plot = fsPlot

        results$featurelist = feature_genes

    } else if (fsMethod ==   "Drop-Out Based Feature Selection - M3Drop" ) {

        counts = sO@assays$RNA@counts

        count_mat = M3Drop::NBumiConvertData(counts, is.counts=TRUE)

        count_mat = M3Drop::NBumiConvertData(count_mat, is.counts=TRUE)

        DANB_fit = M3Drop::NBumiFitModel(count_mat)

        NBDropFS = M3Drop::NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", ntop = featureNumber, suppress.plot=TRUE)

        feature_genes = row.names(NBDropFS)

        sO@assays$RNA@var.features = feature_genes

        top10 = head(feature_genes, 10)

        vfPlot = VariableFeaturePlot(sO)

        fsPlot = LabelPoints(plot = vfPlot, points = top10, repel = TRUE)

        #outputs

        results$data = sO

        results$plot = fsPlot

        results$featurelist = feature_genes


    } else if (fsMethod == "Gene-Gene Correlation Feature Selection (DubStepR)") {

        dubstepR.out = DUBStepR::DUBStepR(input.data = sO@assays$RNA@data, min.cells = 0.05*ncol(sO), optimise.features = TRUE, k = 10, num.pcs = 20, error = 0)

        feature_genes = dubstepR.out$optimal.feature.genes

        top10 = head(feature_genes, 10)

        sO@assays$RNA@var.features = feature_genes

        vfPlot = VariableFeaturePlot(sO)

        fsPlot = LabelPoints(plot = vfPlot, points = top10, repel = TRUE)

        #outputs

        results$data = sO

        results$plot = fsPlot

        results$featurelist = feature_genes

    } else if (fsMethod == "Filter by Gene Expression") {

        sce = as.SingleCellExperiment(sO)

        LogCounts = as.array(logcounts(sce))

        exprsn = rowMeans(LogCounts)

        keep = order(exprsn, decreasing = TRUE)[seq_len(featureNumber)]

        sce[keep, ]

        feature_genes = rownames(sce)

        sO@assays$RNA@var.features = feature_genes

        top10 = head(feature_genes, 10)

        vfPlot = VariableFeaturePlot(sO)

        fsPlot = LabelPoints(plot = vfPlot, points = top10, repel = TRUE)

        #outputs
        results$data = sO

        results$plot = fsPlot

        results$featurelist = feature_genes

    } else if (fsMethod == "Feature Selection by Deviance (Scry)") {

        #require(scry)

        m = GetAssayData(sO, slot = "counts", assay = "RNA")

        #select features

        deviantFeatures = scry::devianceFeatureSelection(m)

        dev_ranked_genes = rownames(sO)[order(deviantFeatures, decreasing = TRUE)]

        topdev = head(dev_ranked_genes, featureNumber)

        top10 = head(topdev, 10)

        sO@assays$RNA@var.features = topdev

        #create plot
        vfPlot = VariableFeaturePlot(sO)

        fsPlot = LabelPoints(plot = vfPlot, points = top10, repel = TRUE)

        #outputs

        results$data = sO

        results$plot = fsPlot

        results$featurelist = topdev


    }

    out = results

}

#The output of this function is always a seurat object that contains the processed data, and a list of selected genes stored in var.features.

