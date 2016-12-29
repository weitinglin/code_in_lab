library(shiny)


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
        )
    ) 
   
)

server <- function(input, output){
    }



shinyApp(ui = ui, server = server)
