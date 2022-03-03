
PCA_DimR = function(sO, outputType) {

    dimReduc = reactiveValues()

    all.genes = rownames(sO)

    sO = ScaleData(sO, features = all.genes)

    sO = RunPCA(sO, features = VariableFeatures(object = sO))

    r_dims = sO[["pca"]]

    dimReduc$data = sO

    if (outputType == "PCA") {

        dimReduc$plot = DimPlot(sO, reduction = "pca")

    } else if (outputType == "Elbow Plot") {

        dimReduc$plot = ElbowPlot(sO)

    }

    dimReduc$method = "pca"

    dimReduc$rdims = r_dims

    out = dimReduc

}

GLM_PCA_DimR = function(sO, L) {

    dimReduc = reactiveValues()

    sO = SeuratWrappers::RunGLMPCA(sO, L = L) #no need to do any data scaling before-hand.

    r_dims = sO[["glmpca"]]

    dimReduc$data = sO

    dimReduc$plot = DimPlot(sO, reduction = "glmpca", pt.size = 0.5)

    dimReduc$rdims = r_dims

    dimReduc$method = "glmpca"

    out = dimReduc

}

#'
#'@import SingleCellExperiment
GLM_PCA_Residuals_DimR = function(sO, outputType) {

    dimReduc = reactiveValues()

    all.genes = rownames(sO)

    #Generate Residual Data

    sce = as.SingleCellExperiment(sO)
    SingleCellExperiment::rowData(sce)$feature_symbol = rownames(sO)

    sce = scry::nullResiduals(sce, assay="counts", type="deviance")

    DCsce = sce[sO@assays$RNA@var.features,]

    nr.data = SummarizedExperiment::assay(DCsce,"binomial_deviance_residuals") #

    #Run Seurat PCA

    metadata = sO@meta.data

    selectedfeatures = sO@assays$RNA@var.features

    sO_2 = CreateSeuratObject(counts = nr.data, meta.data = metadata)

    sO_2@assays$RNA@var.features = selectedfeatures

    sO_2 = ScaleData(sO_2, features = all.genes)

    sO_2 = RunPCA(sO_2, features = VariableFeatures(object = sO_2)) #finish

    #Transfer Data to Main Seurat Object

    sO[["rpca"]] = sO_2[["pca"]]

    #output graphs

    if (outputType ==  "PCA") {

        dimReduc$plot = DimPlot(sO, reduction = "rpca")

    } else if (outputType == "Elbow Plot") {

        dimReduc$plot = ElbowPlot(sO, reduction = "rpca")

    }

    dimReduc$data = sO

    dimReduc$rdims = sO[["rpca"]]

    dimReduc$method = "rpca"

    out = dimReduc


}

#the output of these methods should include a seurat object that contains the dimension reduction results.
