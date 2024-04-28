
#' Launches the GUI
#'
#'
#' @export
#' @import shiny
#' @import shinythemes
#' @import shinycssloaders
launch_FrankenSeq = function() {

    #if necessary, install torch and its dependencies
    torch::install_torch(reinstall = FALSE)

    shinyApp(

        ui = ui,
        server = server

        )

    }


utils::globalVariables(c("nFeature_RNA","percent.mt","kmeans","cluster","avg_log2FC","dl_lfc","dl_sO.markers","dl_minpct'"))
