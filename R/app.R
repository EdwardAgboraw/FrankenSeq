
#' Launches the GUI
#'
#' @return
#' @export
#'
#'
launch_FrankenSeq = function() {

shinyApp(

    ui = ui,
    server = server

    )

}


utils::globalVariables(c("nFeature_RNA","percent.mt","kmeans","cluster","avg_log2FC","dl_lfc","dl_sO.markers","dl_minpct'"))
