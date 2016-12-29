library(shiny)
source("/Users/Weitinglin/Documents/Repository/code_in_lab/00_microarry_function.R")
load("/Users/Weitinglin/Documents/Repository/code_in_lab/ttresult.Rdata")


# Preprocess --------------------------------------------------------------




# UI part -----------------------------------------------------------------


ui <- tagList( 
   
    navbarPage(title = "Project",
        tabPanel(
            title = "Introduction",
            column(12,
                   includeMarkdown("DESCRIPTION.md")
                   )
        ),
        tabPanel(
            title = "Workflow",
            img(src="workflow.png")
        ),
        tabPanel(
            title = "Gene",
            dataTableOutput(outputId = "queryResult")
        )
        )
    ) 
   


# SERVER part -------------------------------------------------------------


server <- function(input, output){
    output$queryResult <-  renderDataTable({
        searchHarmonizome(c("CD44","ALDH1A3","CD9","CDKN2A","DPP4","HSPB1","KIT","NANOG"))
    })
    
    
    }



shinyApp(ui = ui, server = server)
