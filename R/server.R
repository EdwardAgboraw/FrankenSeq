
#The Server function - this is where the actual R-Script goes. This is what tells the app what to do.
library(Seurat)
library(monocle)
library(SingleCellExperiment)
library(SC3)
library(DUBStepR)
library(M3Drop)

#source('clusterMethods.R')
#source('featureSelection.R')
#source('dimensionReduction.R')
#source('LoadData.R')

server = function(input,output) {

    options(shiny.maxRequestSize=300*1024^2) #increase the max file limit

    #QUALITY CONTROL

    RawData = reactiveVal()
    filtered_Data = reactiveVal()

    output$dgraph = renderPlot({

        if(is.null(input$rdata)) {

            return(NULL)

        }

        #input raw data
        rdata = input$rdata
        filePath = input$rdata$datapath

        fileType = input$FileType

        rawdata = LoadData(fileType, filePath)

        RawData(rawdata)

        #Data Processing
        minNum = input$Min_Features
        maxNum = input$Max_Features
        percentMT = input$Mitochondria
        nFeatures = input$nFeatures

        filterdata = FilterData(RawData(), minNum, maxNum, percentMT, nFeatures)

        filtered_Data(filterdata)


        if(input$DPOptions == "Quality Control") {

            QCplots = VlnPlot(RawData(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

            return(QCplots)

        } else if (input$DPOptions == "Data Filtration") {

            NewQCplots = VlnPlot(filtered_Data(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

            return(NewQCplots)

        }

    })

    #FEATURE SELECTION

    FeatureGenes = reactiveVal()# DATASET 2

    FeatureList = reactiveVal()

    output$FSPlot = renderPlot({

        if(is.null(filtered_Data())) {

            return(NULL)

        }

        fs = input$FS
        num_of_features = input$nFeatures

        fs_results = featureSelection(filtered_Data(), fs, num_of_features)

        FeatureGenes(fs_results$data)

        FeatureList(fs_results$featurelist)

        fs_plot = fs_results$plot

        return(fs_plot)


    })

    #DIMENSION REDUCTION

    DimRData = reactiveVal() # DATASET 3

    DimMethod = reactiveVal()

    output$drPlot = renderPlot({

        if(is.null(FeatureGenes())) {

            return(NULL)

        }

        dimR = input$dr

        O_type = input$DRG_options

        L = input$L_value

        if(dimR == "PCA (Seurat)") {

            dr_data = PCA_DimR(FeatureGenes(), O_type)

            DimRData(dr_data$data)

            DimMethod(dr_data$method)

            return(dr_data$plot)

        } else if(dimR == "GLM PCA") {

            dr_data = GLM_PCA_DimR(FeatureGenes(), L)

            DimRData(dr_data$data)

            DimMethod(dr_data$method)

            return(dr_data$plot)

        } else if(dimR == "GLM PCA (Residuals)") {

            dr_data = GLM_PCA_Residuals_DimR(FeatureGenes(), O_type)

            DimRData(dr_data$data)

            DimMethod(dr_data$method)

            return(dr_data$plot)

        }

    })

    #Cluster Validation

    output$cvPlot = renderPlot({

        if(input$CVOptions == "Estimate K with SC3" || is.null(DimRData())) {

            return(NULL)

        }

        library(NbClust)
        library(ggplot2)
        library(factoextra)

        plottype = input$CVOptions

        maxK = input$maxK

        cvData = DimRData()

        clusterValPlot = clusterValidation(cvData, plottype, maxK)

        return(clusterValPlot)
    })

    #SC3 Estimate K

    output$sc3_estimated_k_value <- renderText({

        if(input$CVOptions != "Estimate K with SC3" || is.null(DimRData())) {

            return(NULL)

        }

        clusValData = DimRData()

        k = sc3_estimate(clusValData)

        string = paste("Estimated number of clusters = ", k)

        return(string)


    })


    #CLUSTER ANALYSIS

    Clustered_Data = reactiveVal() # DATASET 4

    output$SClusters = renderPlot({

        if(is.null(DimRData())) {

            return(NULL)

        }

        sO = DimRData()#this object should have both the current feature gene and reduced dimension sets already loaded.

        method1 = DimMethod()

        method = input$CM

        reductionMethod = input$cv

        noDim = input$dimensions

        res = input$resolution

        k_value = input$k_value

        core_number = input$ncore

        linkM = input$LM

        hcOption = input$HCoptions

        ha_fs_option = input$DLoptions

        ha_cl_option= input$DLoptions2

        if (method == "K-Nearest Neighbor (Seurat)") {

            knn = seuratClustering(sO, noDim, res, reductionMethod, method1)

            Clustered_Data(knn$data)

            return(knn$plot)

        } else if (method == "Graph Based Hierarchical Clustering (HGC)") {

            hgc = HGC_clustering(sO, k_value, noDim, reductionMethod, method1)

            Clustered_Data(hgc$data)

            return(hgc$plot)

        } else if (method == "Hierarchical Clustering") {

            hc = h_clustering(sO, linkM, noDim, reductionMethod, method1, k_value, hcOption)

            Clustered_Data(hc$data)

            return(hc$plot)

        } else if (method == "K-means Clustering") {

            km = km_clustering(sO, k_value, noDim, reductionMethod, method1)

            Clustered_Data(km$data)

            return(km$plot)

        }  else if (method == "Density Peak Clustering (Monocle)") {

            cd = monocleClustering(sO, k_value, noDim, reductionMethod, method1)

            Clustered_Data(cd$data)

            return(cd$plot)

        } else if (method == "Consensus Clustering (SC3)") {

            cc = consensus_clustering(sO, k_value, noDim, reductionMethod, method1)

            Clustered_Data(cc$data)

            return(cc$plot)

        } else if (method == "Hierarchical AutoEncoder") {

            ha = autoEncoderClusterring(sO, noDim, k_value, core_number, reductionMethod, method1) #ha_fs_option, ha_cl_option)

            Clustered_Data(ha$data)

            return(ha$plot)
        }

    })

    #OUTPUT DATA

    FinalDataTable = reactiveVal()

    output$table = renderTable({

        if(is.null(Clustered_Data())) {

            return(NULL)

        }

        tableChoice = input$TableOptions

        sO = Clustered_Data()

        if (tableChoice == "Summary Report") {

            Dsummary = table(Idents(sO))

            Dsummary = as.data.frame(Dsummary)

            colnames(Dsummary) = c("Cluster", "Frequency")

            FinalDataTable(Dsummary)

            return(Dsummary)

        } else {

            results = as.data.frame(sO@active.ident)

            colnames(results) = "Cluster"

            results = cbind(rownames(results), data.frame(results, row.names = NULL))

            colnames(results) = c("Cell Label", "Cluster")

            #results = results[order(results$Cluster), ]

            FinalDataTable(results)

            return(results)

        }

    })

    output$Pipeline <- renderUI({

        str1 <- paste("<B>Feature Selection Algorithm:<B> ", input$FS)
        str2 <- paste("<B>Dimension Reduction Method:<B> ", input$dr)
        str3 <- paste("<B>Unsupervised Clustering Method:<B> ", input$CM)

        HTML(paste(str1, str2, str3, sep = '<br/>'))

    })

    output$dlb <- downloadHandler(

        filename = "FinalData.csv",

        contentType = "csv",

        content = function(file) {

            write.csv(FinalDataTable(), file, row.names = FALSE)
        }
    )

    output$FSdlb <- downloadHandler(

        filename = "Selected_Features.txt",

        contentType = "csv",

        content = function(file) {

            write.csv(FeatureList(), file, row.names = FALSE)
        }
    )

    output$SOdlb <- downloadHandler(

        filename = "FinalSeuratObject.txt",

        contentType = "RDS",

        content = function(file) {

            saveRDS(Clustered_Data(), file)
        }
    )

    #Differentially Expressed Gene Analysis

    output$degPlot = renderPlot({

        if(is.null(Clustered_Data())) {

            return(NULL)

        }

        sO = Clustered_Data()

        lfc = input$LFC

        minpct = input$minPC

        sO.markers <- FindAllMarkers(sO, only.pos = TRUE, min.pct = minpct, logfc.threshold = lfc)

        sO.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) -> top10

        degHeatMap = DoHeatmap(sO, features = top10$gene) + NoLegend()

        return(degHeatMap)

    })

    output$degTable = renderTable({

        sO = Clustered_Data()

        lfc = input$LFC

        minpct = input$minPC

        sO.markers <- FindAllMarkers(sO, only.pos = TRUE, min.pct = minpct, logfc.threshold = lfc)

        DEG_Table = sO.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

        return(DEG_Table)

    })

    #Deep Learning

    DeepLearning_Data = reactiveVal()

    output$DL_Plot = renderPlot({

        if(is.null(DimRData())) {

            return(NULL)

        }

        #Variables

        sO = DimRData()

        k_value = input$DL_Knumber

        noDim = input$DL_dimensions

        core_number = input$DL_ncore

        method1 = DimMethod()

        reductionMethod = input$DL_cv

        #DL Function Call

        ha = autoEncoderClusterring(sO, noDim, k_value, core_number, reductionMethod, method1) #ha_fs_option, ha_cl_option)

        DeepLearning_Data(ha$data)

        return(ha$plot)

    })


    output$dl_seurat <- downloadHandler(

        filename = "DL_SeuratObject.txt",

        contentType = "RDS",

        content = function(file) {

            saveRDS(DeepLearning_Data(), file)
        }
    )



    DL_DataTable = reactiveVal()

    output$DL_Table = renderTable({

        if(is.null(DeepLearning_Data())) {

            return(NULL)

        }

        sO = DeepLearning_Data()

        tableChoice = input$DL_TableOptions

        if (tableChoice == "Summary Report") {

            Dsummary = table(Idents(sO))

            Dsummary = as.data.frame(Dsummary)

            colnames(Dsummary) = c("Cluster", "Frequency")

            DL_DataTable(Dsummary)

            return(Dsummary)

       } else {

            results = as.data.frame(sO@active.ident)

            colnames(results) = "Cluster"

            results = cbind(rownames(results), data.frame(results, row.names = NULL))

            colnames(results) = c("Cell Label", "Cluster")

            results = results[order(results$Cluster), ]

            DL_DataTable(results)

            return(results)
        }

    })


    output$dl_table_download <- downloadHandler(

        filename = "DL_Data.csv",

        contentType = "csv",

        content = function(file) {

            write.csv(DL_DataTable(), file, row.names = FALSE)
        }
    )


    output$DL_Heatmap = renderPlot({

        if(is.null(DeepLearning_Data())) {

            return(NULL)

        }

        #Variables

        sO = DeepLearning_Data()

        dl_lfc = input$DL_LFC

        dl_minpct = input$DL_minPC

        dl_sO.markers <- FindAllMarkers(sO, only.pos = TRUE, min.pct = dl_minpct, logfc.threshold = dl_lfc)

        dl_sO.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) -> top10

        dl_degHeatMap = DoHeatmap(sO, features = top10$gene) + NoLegend()

        return(dl_degHeatMap)

    })

    DL_Biomarker_Table = reactiveVal()

    output$DL_Biomarkers = renderTable({

        if(is.null(DeepLearning_Data())) {

            return(NULL)

        }

        #Variables

        sO = DeepLearning_Data()

        dl_lfc = input$DL_LFC

        dl_minpct = input$DL_minPC

        dl_sO.markers <- FindAllMarkers(sO, only.pos = TRUE, min.pct = dl_minpct, logfc.threshold = dl_lfc)

        dl_sO.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) -> top10

        dl_DEG_Table = dl_sO.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

        DL_Biomarker_Table(dl_DEG_Table)

        return(dl_DEG_Table)
    })


    output$dl_bmTable <- downloadHandler(

        filename = "DL_BiomarkerTable.csv",

        contentType = "csv",

        content = function(file) {

            write.csv(DL_Biomarker_Table(), file, row.names = FALSE)
        }
    )

}
