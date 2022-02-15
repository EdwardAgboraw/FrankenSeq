
LoadData = function(filetype, filepath) {

    DataProcessing = reactiveValues()

    if(filetype == "RDS (Existing Seurat Object)") {

        rawdata = readRDS(file = filepath)

        rawdata[["percent.mt"]] = PercentageFeatureSet(rawdata, pattern = "^MT-")

    } else if (filetype == "CSV (Gene x Cell Table)") {

        data = read.csv(file = filepath)

        rawdata = CreateSeuratObject(counts = data) #assay = "RNA", names.field = NULL)

        rawdata[["percent.mt"]] = PercentageFeatureSet(rawdata, pattern = "^MT-")

    } else if (filetype == "SingleCellExperiment Object (RDS)") {

        rawdata = readRDS(file = filepath)

        rawdata = as.Seurat(rawdata, data = NULL)

        rawdata[["percent.mt"]] = PercentageFeatureSet(rawdata, pattern = "^MT-")
    }

    out = rawdata

}

FilterData = function(rawdata, minNum, maxNum, percentMT, nFeatures) {

    FilteredData = subset(rawdata, subset = nFeature_RNA > minNum & nFeature_RNA < maxNum & percent.mt < percentMT)
    FilteredData = NormalizeData(FilteredData, normalization.method = "LogNormalize", scale.factor = 10000)

    FilteredData = FindVariableFeatures(FilteredData, selection.method = "vst", nfeatures = nFeatures) #necessary for compatiblity with downstream algorithms

    out = FilteredData
}
