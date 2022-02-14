
#' Launches the GUI
#'
#' @return
#' @export
#'
#' @examples
launch_FrankenSeq = function() {

    require(dplyr)
    require(shiny)
    require(Seurat)
    require(shinythemes)
    require(shinycssloaders)
    require(magrittr)

    require(DUBStepR)

    require(VGAM)

    require(dynamicTreeCut)

    require(NbClust)
    require(ggplot2)
    require(factoextra)

    require(scDHA)

    require(M3Drop)
    require(SingleCellExperiment)
    require(HGC)
    require(monocle)
    require(SC3)
    require(bluster)
    require(scry)
    require(SeuratWrappers)

shinyApp(

    ui = ui,
    server = server

    )

}


