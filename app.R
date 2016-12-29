library(shiny)
library(VennDiagram)


library(shiny)

ui <- fluidPage( 
    
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = 'file1',label = 'choose file to upload',
                      accept =c('text/csv',"text/comma-separated-values","text/plain")
            )),
        mainPanel(
            dataTableOutput(outputId = "demo1")
        )
    )
)

source("/Users/Weitinglin/Documents/R_scripts/Lab/temple_shiny.R")
server <- function(input, output){
    
    rawdata <- reactive({input$file1})
    
    output$demo1 <- renderDataTable({
        if (is.null(rawdata()))
            return(NULL)
        
        query  <- read.csv(rawdata()$datapath, header = FALSE, sep=" ")
        result <- KEGG(query$V1)
        result
    })
    
    
}



shinyApp(ui = ui, server = server)
