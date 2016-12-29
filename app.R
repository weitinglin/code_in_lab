library(shiny)
library(VennDiagram)


library(shiny)

ui <- tagList( 
   
    navbarPage(title = "Project",
        tabPanel(
            "Introduction",
            column(12,
                   includeMarkdown("DESCRIPTION.md")
                   )
        )
    ) 
   
)

server <- function(input, output){
    }



shinyApp(ui = ui, server = server)
