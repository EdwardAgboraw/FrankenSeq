
library(shiny)


source('ui.R')
source('server.R')

launch_FrankenSeq = function() {

shinyApp(

    ui = ui_2,
    server = server

    )

}


