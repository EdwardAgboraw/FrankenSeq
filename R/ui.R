
library(shinythemes)
library(shinycssloaders)
library(magrittr)

ui = fluidPage(theme = shinytheme("spacelab"),
               navbarPage(
                   "FrankenSeq",

                   tabPanel("Data Processing", # Data Processing Tab

                            sidebarPanel(

                                radioButtons(inputId = "Help_QC", label = "Select an Option", choices = c("Help", "Run Tab"), selected = "Help"),

                                selectInput(inputId = "FileType", label = "Select Input File Type", choices = c("RDS (Existing Seurat Object)", "CSV (Gene x Cell Table)", "SingleCellExperiment Object (RDS)", "Read Count Data (CSV)"), selected = "RDS (Existing Seurat Object)"),

                                fileInput(inputId = "rdata", label = "Select an input file", accept = ".rds", buttonLabel = "Browse...", placeholder = "No file selected..."),

                                radioButtons(inputId = "DPOptions", label = "Select an Analysis", choices = c("Quality Control", "Data Filtration")),

                                #Conditional Panel
                                conditionalPanel(
                                    condition = "input.DPOptions == 'Data Filtration'",
                                    numericInput(inputId = "Max_Features", label = "Maximum Number of Features", max = 5000, min = 1, value = 3000, step = 1),
                                    numericInput(inputId = "Min_Features", label = " Minimum Number of Features", max = 5000, min = 1, value = 200, step = 1),
                                    numericInput(inputId = "Mitochondria", label = " Max % Mitochondrial DNA", max = 100, min = 1, value = 5, step = 1),
                                ),

                            ),

                            mainPanel(

                                conditionalPanel( #QC Help Text

                                    condition = "input.Help_QC == 'Help'",

                                    textOutput("QC_Help_Text")
                                ),

                                plotOutput("dgraph") %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),

                            )
                   ),

                   tabPanel("Feature Selection", #Feature Selection Tab

                            sidebarPanel(

                                radioButtons(inputId = "Help_FS", label = "Select an Option", choices = c("Help", "Run Tab"), selected = "Help"),

                                selectInput(inputId = "FS", label = "Choose a Feature Selection method", choices = c("HVG - Seurat (vst)", "HVG - Seurat (mvp)", "HVG - Seurat (Dispersion)", "Drop-Out Based Feature Selection - M3Drop", "Gene-Gene Correlation Feature Selection (DubStepR)", "Feature Selection by Deviance (Scry)","Filter by Gene Expression"), selected = "HVG - Seurat (vst)"),

                                #conditional Feature Selection Panels
                                conditionalPanel(
                                    condition = "input.FS == 'HVG - Seurat (vst)' || input.FS == 'HVG - Seurat (Dispersion)' || input.FS == 'Feature Selection by Deviance (Scry)' || input.FS == 'Filter by Gene Expression' || input.FS == 'Drop-Out Based Feature Selection - M3Drop'",
                                    numericInput("nFeatures", "Number of Top Features", min = 1, value = 1556, step = 1)

                                ),

                            ),

                            mainPanel(

                                conditionalPanel( #DL Help Text

                                    condition = "input.Help_FS == 'Help'",

                                    textOutput("FS_Help_Text")
                                ),

                                plotOutput("FSPlot") %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),
                            )

                   ),

                   tabPanel("Dimension Reduction", #Dimension Reduction Tab

                            sidebarPanel(

                                radioButtons(inputId = "Help_DR", label = "Select an Option", choices = c("Help", "Run Tab"), selected = "Help"),

                                selectInput(inputId = "dr", label = "Choose a Dimension Reduction method", choices = c("PCA (Seurat)", "GLM PCA (Residuals)", "GLM PCA"), selected = "PCA (Seurat)"),

                                conditionalPanel(
                                    condition = "input.dr == 'GLM PCA'",

                                    numericInput("L_value", "Choose the number of dimensions to return: ", min = 2, value = 10, step = 1),

                                ),

                            ),

                            mainPanel(

                                conditionalPanel(
                                    condition = "input.dr == 'PCA (Seurat)' || input.dr == 'GLM PCA (Residuals)' ",

                                    selectInput("DRG_options", "Select a visualition option", choices = c("Elbow Plot", "PCA"), selected = "PCA"),

                                ),

                                conditionalPanel( #DR Help Text

                                    condition = "input.Help_DR == 'Help'",

                                    textOutput("DR_Help_Text")
                                ),

                                plotOutput("drPlot") %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),
                            )

                   ),


                   tabPanel("Cluster Validation", #Cluster Validation Tab

                            sidebarPanel(

                                radioButtons(inputId = "Help_CV", label = "Select an Option", choices = c("Help", "Run Tab"), selected = "Help"),

                                radioButtons(inputId = "CVOptions", label = "Select a Plot", choices = c("Estimate K with SC3","Elbow Plot", "Silhouette Plot")),

                                numericInput(inputId = "maxK", label = "Select maxinum number of clusters to evaluate;", value = 15, step = 1)

                            ),

                            mainPanel(

                                conditionalPanel( #CV Help Text

                                    condition = "input.Help_CV == 'Help'",

                                    textOutput("CV_Help_Text")
                                ),

                                textOutput("sc3_estimated_k_value") %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),

                                plotOutput("cvPlot")  %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),

                            ),

                   ),



                   tabPanel("Cluster Analysis", #Cluster Analysis Tab

                            sidebarPanel(

                                radioButtons(inputId = "Help_CA", label = "Select an Option", choices = c("Help", "Run Tab"), selected = "Help"),

                                selectInput(inputId = "CM", label = "Choose an Unsupervised Clustering Algorithm", choices = c("K-Nearest Neighbor (Seurat)", "Graph Based Hierarchical Clustering (HGC)", "Hierarchical Clustering", "K-means Clustering", "Density Peak Clustering (Monocle)", "Consensus Clustering (SC3)")),

                                numericInput(inputId = "dimensions", label = "Dimensionality: ", value = 10, max = 30, min = 2, step = 1),

                                #Conditional Panels

                                conditionalPanel(

                                    condition = "input.CM == 'K-Nearest Neighbor (Seurat)'",

                                    numericInput(inputId = "resolution", label = "Resolution Parameter", value = 0.5, max = 1.5, step = 0.05),
                                ),


                                conditionalPanel(

                                    condition = "input.CM == 'Graph Based Hierarchical Clustering (HGC)' || input.CM == 'K-means Clustering' || input.CM == 'Consensus Clustering (SC3)' || input.CM == 'Density Peak Clustering (Monocle)' || input.CM == 'Hierarchical Clustering'",

                                    numericInput(inputId = "k_value", label = "Desired Number of Clusters", value = 4, min = 1, step = 1),


                                ),

                                conditionalPanel(

                                    condition = "input.CM == 'Hierarchical Clustering'",

                                    selectInput(inputId = "LM", label = "Choose a Linkage Method", choices = c("ward.D", "ward.D2", "complete", "mcquitty", "average", "median", "centroid"), selected = "ward.D2"),

                                    radioButtons(inputId = "HCoptions", label = "Select a methodology", choices = c("Use determined Cluster Number","Use dynamic tree cut")),

                                ),


                            ),

                            mainPanel(

                                selectInput(inputId = "cv", label = "T-SNE or UMAP Cluster Visualization", choices = c("umap", "tsne"), selected = "umap"),

                                conditionalPanel( #CA Help Text

                                    condition = "input.Help_CA == 'Help'",

                                    textOutput("CA_Help_Text")
                                ),


                                plotOutput("SClusters") %>% withSpinner(color="#0dc5c1", hide.ui = FALSE)
                            ),
                   ),

                   tabPanel("Results", #Clustering Results Tab

                            sidebarPanel(

                                radioButtons(inputId = "Help_R", label = "Select an Option", choices = c("Help", "Run Tab"), selected = "Help"),

                                radioButtons(inputId = "TableOptions", label = "Select a Table", choices = c("Summary Report", "Full Data Table")),

                                helpText("Current Pipeline"),
                                htmlOutput("Pipeline"),

                                helpText("Click here to download the Table as a CSV"),
                                downloadButton("dlb", "Download Table"),

                                helpText("Click here to download a List of the Current Selected Features"),
                                downloadButton("FSdlb", "Download Features"),

                                helpText("Click here to download the Current Seurat Object as an RDS"),
                                downloadButton("SOdlb", "Download Seurat Object"),
                            ),

                            mainPanel(

                                conditionalPanel( #DL Help Text

                                    condition = "input.Help_R == 'Help'",

                                    textOutput("R_Help_Text")
                                ),

                                tableOutput("table"),
                            ),

                   ),

                   tabPanel("DEG Analysis", #DEG ANALYSIS Tab

                            sidebarPanel(

                                radioButtons(inputId = "Help_DA", label = "Select an Option", choices = c("Help", "Run Tab"), selected = "Help"),

                                radioButtons(inputId = "DEG_Options", label = "Select an Option", choices = c("Heatmap","Table")),

                                numericInput(inputId = "LFC", label = "Select a Log Fold Change Threshold;", value = 0.25, step = 0.05),

                                numericInput(inputId = "minPC", label = "Select a minimum percent expression value; ", value = 0.25, step = 0.05),


                                conditionalPanel(

                                    condition = "input.DEG_Options == 'Table'",

                                    helpText("Click here to download the Table as a CSV"),
                                    downloadButton("deg_biomarkers", "Download Table"),

                                ),



                            ),

                            mainPanel(

                                conditionalPanel( #DA Help Text

                                    condition = "input.Help_DA == 'Help'",

                                    textOutput("DA_Help_Text")
                                ),

                                conditionalPanel(

                                    condition = "input.DEG_Options == 'Heatmap'",

                                    plotOutput('degPlot', width = "100%") %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),

                                ),

                                conditionalPanel(

                                    condition = "input.DEG_Options == 'Table'",

                                    tableOutput('degTable') %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),

                                ),

                            ),

                   ),

                   tabPanel("Deep Learning", #Deep Learning Tab

                            sidebarPanel(

                                radioButtons(inputId = "Help_DL", label = "Select an Option", choices = c("Help", "Run Tab"), selected = "Help"),

                                radioButtons(inputId = "DL_Options", label = "Select an Option", choices = c("Cluster Analysis","Cluster Data Table", "Heatmap", "Cluster Biomarker Table"), selected = "Cluster Analysis"),


                                conditionalPanel(

                                    condition = "input.DL_Options == 'Cluster Analysis' || input$Help_DL == 'Help'",

                                    numericInput(inputId = "DL_Knumber", label = "Select the desired number of clusters;", value = 8, step = 1, min = 1),

                                    numericInput(inputId = "DL_dimensions", label = "Dimensionality: ", value = 10, max = 30, min = 2, step = 1),

                                    numericInput(inputId = "DL_ncore", label = "Choose the number of cores to use in the analysis", value = 4, min = 1, step = 1),

                                    selectInput(inputId = "DL_cv", label = "T-SNE, or UMAP Cluster Visualization", choices = c("umap", "tsne"), selected = "umap"),

                                    helpText("Click here to download the complete Seurat Object as an RDS file"),
                                    downloadButton("dl_seurat", "Download Object"),

                                ),

                                conditionalPanel(

                                    condition = "input.DL_Options == 'Cluster Data Table'",

                                    radioButtons(inputId = "DL_TableOptions", label = "Select a Table", choices = c("Summary Report", "Full Data Table")),

                                    helpText("Click here to download the current Table Option as a CSV"),
                                    downloadButton("dl_table_download", "Download Table"),

                                ),


                                conditionalPanel(

                                    condition = "input.DL_Options == 'Heatmap' || input.DL_Options == 'Cluster Biomarker Table'",

                                    numericInput(inputId = "DL_LFC", label = "Select a Log Fold Change Threshold;", value = 0.25, step = 0.05),

                                    numericInput(inputId = "DL_minPC", label = "Select a minimum percent expression value; ", value = 0.25, step = 0.05),


                                ),

                                conditionalPanel(

                                    condition = "input.DL_Options == 'Cluster Biomarker Table'",

                                    helpText("Click here to download the Cluster Biomarker Table as a PNG"),
                                    downloadButton("dl_bmTable", "Download Table"),


                                ),

                            ),

                            mainPanel(

                                conditionalPanel( #DL Help Text

                                    condition = "input.Help_DL == 'Help'",

                                    textOutput("DL_Help_Text")
                                ),

                                conditionalPanel( #DL Cluster Plot

                                    condition = "input.DL_Options == 'Cluster Analysis'",

                                    plotOutput('DL_Plot', width = "100%") %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),

                                ),

                                conditionalPanel( #DL Data Tables

                                    condition = "input.DL_Options == 'Cluster Data Table'",

                                    tableOutput('DL_Table') %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),

                                ),

                                conditionalPanel( #DL Heatmap

                                    condition = "input.DL_Options == 'Heatmap'",

                                    plotOutput('DL_Heatmap', width = "100%") %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),

                                ),

                                conditionalPanel( #DL Biomarker Table

                                    condition = "input.DL_Options == 'Cluster Biomarker Table'",

                                    tableOutput('DL_Biomarkers') %>% withSpinner(color="#0dc5c1", hide.ui = FALSE),

                                ),



                            ),

                   ),


               )

)
