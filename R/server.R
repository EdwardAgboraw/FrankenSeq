#'The Server Function
#'
#' @param input standard shiny parameter
#'
#' @param output standard shiny parameter
#'
#'@import magrittr
#'@import dplyr
#'@import shiny
#'@return None
server = function(input,output) {

    options(shiny.maxRequestSize=300*1024^2) #increase the max file limit

    #QUALITY CONTROL

    rawData <- reactiveVal()
    filteredData <- reactiveVal()

    output$dgraph = renderPlot({

        if(input$Help_QC == "Help") {

            return(NULL)

        }

        #input raw data
        rdata <- input$rdata
        filePath <- input$rdata$datapath

        fileType <- input$FileType

        if (fileType == 'PBMC3k Test Data') {

            inputData = LoadTestData()

        } else {

            inputData <- LoadData(fileType, filePath)

        }


        rawData(inputData)

        #Data Processing
        minNum = input$Min_Features
        maxNum = input$Max_Features
        percentMT = input$Mitochondria
        nFeatures = input$nFeatures

        cleanData <- FilterData(rawData(), minNum, maxNum, percentMT, nFeatures)

        filteredData(cleanData)


        if(input$DPOptions == "Quality Control") {

            qcPlots <- VlnPlot(rawData(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

            return(qcPlots)

        } else if (input$DPOptions == "Data Filtration") {

            newQCplots <- VlnPlot(filteredData(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

            return(newQCplots)

        }

    })

    #FEATURE SELECTION

    featureGenes <- reactiveVal()# DATASET 2

    featureList <- reactiveVal()

    output$FSPlot = renderPlot({

        if(is.null(filteredData()) || input$Help_FS == "Help") {

            return(NULL)

        }

        fs <- input$FS
        featureNumber <- input$nFeatures

        fs_results <- featureSelection(filteredData(), fs, featureNumber)

        featureGenes(fs_results$data)

        featureList(fs_results$featurelist)

        fs_plot <- fs_results$plot

        return(fs_plot)


    })

    #DIMENSION REDUCTION

    dimRData <- reactiveVal() # DATASET 3

    dimRMethod <- reactiveVal()

    output$drPlot = renderPlot({

        if(is.null(featureGenes()) || input$Help_DR == "Help") {

            return(NULL)

        }

        dimR <- input$dr

        graphType_DR <- input$DRG_options

        L <- input$L_value

        if(dimR == "PCA (Seurat)") {

            dr_data <- PCA_DimR(featureGenes(), graphType_DR)

            dimRData(dr_data$data)

            dimRMethod(dr_data$method)

            return(dr_data$plot)

        } else if(dimR == "GLM PCA") {

            dr_data <- GLM_PCA_DimR(featureGenes(), L)

            dimRData(dr_data$data)

            dimRMethod(dr_data$method)

            return(dr_data$plot)

        } else if(dimR == "GLM PCA (Residuals)") {

            dr_data <- GLM_PCA_Residuals_DimR(featureGenes(), graphType_DR)

            dimRData(dr_data$data)

            dimRMethod(dr_data$method)

            return(dr_data$plot)

        }

    })

    #Cluster Validation

    output$cvPlot = renderPlot({

        if(input$CVOptions == "Estimate K with SC3" || is.null(dimRData()) || input$Help_CV == "Help") {

            return(NULL)

        }

        graphType_CV <- input$CVOptions

        maxK <- input$maxK

        cvData <- dimRData()

        clusterValPlot <- clusterValidation(cvData, graphType_CV, maxK)

        return(clusterValPlot)
    })

    #SC3 Estimate K

    output$sc3_estimated_k_value <- renderText({

        if(input$CVOptions != "Estimate K with SC3" || is.null(dimRData()) || input$Help_CV == "Help") {

            return(NULL)

        }

        clusValData <- dimRData()

        k <- sc3_estimate(clusValData)

        string <- paste("Estimated number of clusters = ", k)

        return(string)


    })


    #CLUSTER ANALYSIS

    clusterData <- reactiveVal() # DATASET 4

    output$SClusters = renderPlot({

        if(is.null(dimRData()) || input$Help_CA == "Help") {

            return(NULL)

        }

        sO <- dimRData() #this object should have both the current feature gene and reduced dimension sets already loaded.

        reductionMethod <- dimRMethod() #formerly method1

        clusteringMethod <- input$CM #formerly method

        graphType_CA <- input$cv #formerly reductionMethod

        noDim <- input$dimensions

        res <- input$resolution

        kValue <- input$k_value

        coreNumber <- input$ncore

        linkM <- input$LM

        hcOption <- input$HCoptions

        ha_fs_option <- input$DLoptions

        ha_cl_option <- input$DLoptions2

        if (clusteringMethod == "K-Nearest Neighbor (Seurat)") {

            knn <- seuratClustering(sO, noDim, res, graphType_CA, reductionMethod)

            clusterData(knn$data)

            return(knn$plot)

        } else if (clusteringMethod == "Graph Based Hierarchical Clustering (HGC)") {

            hgc <- HGC_clustering(sO, kValue, noDim, graphType_CA, reductionMethod)

            clusterData(hgc$data)

            return(hgc$plot)

        } else if (clusteringMethod == "Hierarchical Clustering") {

            hc <- h_clustering(sO, linkM, noDim, graphType_CA, reductionMethod, kValue, hcOption)

            clusterData(hc$data)

            return(hc$plot)

        } else if (clusteringMethod == "K-means Clustering") {

            km <- km_clustering(sO, kValue, noDim, graphType_CA, reductionMethod)

            clusterData(km$data)

            return(km$plot)

        }  else if (clusteringMethod == "Density Peak Clustering (Monocle)") {

            cd <- monocleClustering(sO, kValue, noDim, graphType_CA, reductionMethod)

            clusterData(cd$data)

            return(cd$plot)

        } else if (clusteringMethod == "Consensus Clustering (SC3)") {

            cc <- consensus_clustering(sO, kValue, noDim, graphType_CA, reductionMethod)

            clusterData(cc$data)

            return(cc$plot)

    }
        })

    #OUTPUT DATA

    resultsTable = reactiveVal()

    output$table = renderTable({

        if(is.null(clusterData()) || input$Help_R == "Help") {

            return(NULL)

        }

        tableChoice = input$TableOptions

        sO = clusterData()

        resultsTable = generateTable(sO, tableChoice)

        return(resultsTable)

    })

    output$Pipeline <- renderUI({

        str1 <- paste("<B>Feature Selection Algorithm:<B> ", input$FS)
        str2 <- paste("<B>Dimension Reduction Method:<B> ", input$dr)
        str3 <- paste("<B>Unsupervised Clustering Method:<B> ", input$CM)

        HTML(paste(str1, str2, str3, sep = '<br/>'))

    })

    output$dlb <- downloadHandler(

        filename <- "FinalData.csv",

        contentType <- "csv",

        content <- function(file) {

            utils::write.csv(resultsTable(), file, row.names = FALSE)
        }
    )

    output$FSdlb <- downloadHandler(

        filename <- "Selected_Features.txt",

        contentType <- "csv",

        content <- function(file) {

            utils::write.csv(featureList(), file, row.names = FALSE)
        }
    )

    output$SOdlb <- downloadHandler(

        filename <- "FinalSeuratObject.txt",

        contentType<- "RDS",

        content <- function(file) {

            saveRDS(clusterData(), file)
        }
    )

    #Differentially Expressed Gene Analysis

    output$degPlot = renderPlot({

        if(is.null(clusterData()) || input$Help_DA == "Help") {

            return(NULL)

        }

        sO <- clusterData()

        lfc <- input$LFC

        minPCT <- input$minPC

        outputType <- input$DEG_Options

        degResults = degAnalysis(sO, lfc, minPCT, outputType)

        return(degResults)

    })

    degMarkerTable = reactiveVal()

    output$degTable = renderTable({

        if(is.null(clusterData()) || input$Help_DA == "Help") {

            return(NULL)

        }

        sO <- clusterData()

        lfc <- input$LFC

        minPCT <- input$minPC

        outputType <- input$DEG_Options

        degResults = degAnalysis(sO, lfc, minPCT, outputType)

        degMarkerTable(degResults)

        return(degResults)

    })

    output$deg_biomarkers <- downloadHandler(

        filename <- "ClusterBiomarkers.csv",

        contentType <- "csv",

        content <- function(file) {

            utils::write.csv(degMarkerTable, file, row.names = FALSE)
        }
    )


    #Deep Learning

    deepLearningData <- reactiveVal()

    output$DL_Plot <- renderPlot({

        if(is.null(dimRData()) || input$Help_DL == "Help") {

            return(NULL)

        }

        #Variables

        sO <- dimRData()

        kValue <- input$DL_Knumber

        noDim <- input$DL_dimensions

        coreNumber <- input$DL_ncore

        reductionMethod <- dimRMethod() #formerly method1

        graphType_DL <- input$DL_cv #formerly reductionMethod

        #DL Function Call

        deepLearningResults <- autoEncoderClusterring(sO, noDim, kValue, coreNumber, graphType_DL, reductionMethod)

        deepLearningData(deepLearningResults$data)

        return(deepLearningResults$plot)

    })


    output$dl_seurat <- downloadHandler(

        filename <- "DL_SeuratObject.txt",

        contentType <- "RDS",

        content <- function(file) {

            saveRDS(deepLearningData(), file)
        }
    )



    resultsTable_DL = reactiveVal()

    output$DL_Table = renderTable({

        if(is.null(deepLearningData()) || input$Help_DL == "Help") {

            return(NULL)

        }

        sO = deepLearningData()

        tableChoice = input$DL_TableOptions

        if (tableChoice == "Summary Report") {

            Dsummary <- table(Idents(sO))

            Dsummary <- as.data.frame(Dsummary)

            colnames(Dsummary) <- c("Cluster", "Frequency")

            resultsTable_DL(Dsummary)

            return(Dsummary)

       } else {

            results <- as.data.frame(sO@active.ident)

            colnames(results) <- "Cluster"

            results <- cbind(rownames(results), data.frame(results, row.names = NULL))

            colnames(results) <- c("Cell Label", "Cluster")

            results <- results[order(results$Cluster), ]

            resultsTable_DL(results)

            return(results)
        }

    })


    output$dl_table_download <- downloadHandler(

        filename <- "DL_Data.csv",

        contentType <- "csv",

        content <- function(file) {

            utils::write.csv(resultsTable_DL(), file, row.names = FALSE)
        }
    )


    output$DL_Heatmap = renderPlot({

        if(is.null(deepLearningData()) || input$Help_DL == "Help") {

            return(NULL)

        }

        #Variables

        sO <- deepLearningData()

        lfc_DL <- input$DL_LFC

        minPCT_DL <- input$DL_minPC

        degHeatMap_DL = degAnalysis(sO, lfc_DL, minPCT_DL, "Heatmap")

        #sO.markers_DL <- FindAllMarkers(sO, only.pos = TRUE, min.pct = minPCT_DL, logfc.threshold = lfc_DL)

        #sO.markers_DL %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) -> top10

        #degHeatMap_DL <- DoHeatmap(sO, features = top10$gene) + NoLegend()

        return(degHeatMap_DL)

    })

    biomarkerTable_DL = reactiveVal()

    output$DL_Biomarkers = renderTable({

        if(is.null(deepLearningData()) || input$Help_DL == "Help") {

            return(NULL)

        }

        #Variables

        sO <- deepLearningData()

        lfc_DL <- input$DL_LFC

        minPCT_DL <- input$DL_minPC

        #sO.markers_DL <- FindAllMarkers(sO, only.pos = TRUE, min.pct = dl_minpct, logfc.threshold = dl_lfc)

        #sO.markers_DL %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) -> top10

        #degTable_DL <- dl_sO.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

        degTable_DL = degAnalysis(sO, lfc_DL, minPCT_DL, "Table")

        biomarkerTable_DL(degTable_DL)

        return(degTable_DL)
    })


    output$dl_bmTable <- downloadHandler(

        filename <- "DL_BiomarkerTable.csv",

        contentType <- "csv",

        content <- function(file) {

            utils::write.csv(biomarkerTable_DL(), file, row.names = FALSE)
        }
    )

    #The Help Text

    output$QC_Help_Text <- renderText(

        "This is the Data Processing Tab. Accepted inputs include Seurat/Single Cell Experiment Objects (as RDS files)
        and Gene x Cell matrices (as CSV files). 'Quality Control' generates violin plots visualizing the QC metrics of the raw input data.
        'Data Filtration' trims cells from the data according to user defined Feature Count and % Mitochondrial DNA limits."

    )

    output$FS_Help_Text <- renderText(

        "This is the Feature Selection Tab. "

    )

    output$DR_Help_Text <- renderText(

        "This is the Dimension Reduction Tab."

    )

    output$CV_Help_Text <- renderText(

        "This is the Cluster Validation Tab. This tab estimates the number of clusters in the data using either the Elbow or Silhouette Methods.
        Cluster number can also be directly estimated using the Tracy-Widom Theory of Random Matrices (via SC3). This tab is both time-intensive and optional -
        this part of the pipeline can safely be skipped."

    )

    output$CA_Help_Text <- renderText(

        "This is the Cluster Analysis Tab."

    )

    output$R_Help_Text <- renderText(

        "This is the R Help Text."

    )

    output$DA_Help_Text <- renderText(

        "This is the Differentially Expressed Gene (DEG) Analysis Tab."

    )

    output$DL_Help_Text <- renderText(

        "This is the Deep Learning Tab. This Tab represents a fully independent analysis pipeline."

    )

}
