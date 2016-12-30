library(shiny)
source("/Users/Weitinglin/Documents/Repository/code_in_lab/00_microarry_function.R")
load("/Users/Weitinglin/Documents/Repository/code_in_lab/total_ttest_result.Rdata")
load("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/data/intermediate/hgu133plus2.RData")
load("/Users/Weitinglin/Documents/Repository/code_in_lab/total_probe_dataframe.Rdata")


# Preprocess --------------------------------------------------------------

# UP-regulation 


t.greater.CLF_CLS <- t.greater.CLF_CLS %>% mutate(Case = "P6+fibroblast > P6 nofilter", Method = "Nofilter")
t.greater.Sphere_CLF <- t.greater.Sphere_CLF %>% mutate(Case = "P6-Sphere > P6+fibroblast nofilter", Method = "Nofilter")
t.greater.Sphere_CLS <- t.greater.Sphere_CLS %>% mutate(Case = "P6-Sphere > p6 nofilter", Method = "Nofilter")
total.t.greater <- bind_rows(t.greater.CLF_CLS, t.greater.Sphere_CLF, t.greater.Sphere_CLS)


# DOWN-regulation 
t.less.CLF_CLS <- t.less.CLF_CLS %>% mutate(Case = "P6+fibroblast < P6 nofilter", Method = "Nofilter")
t.less.Sphere_CLF <- t.less.Sphere_CLF %>% mutate(Case = "P6-Sphere < P6+fibroblast nofilter", Method = "Nofilter")
t.less.Sphere_CLS <- t.less.Sphere_CLS %>% mutate(Case = "P6-Sphere < p6 nofilter", Method = "Nofilter")
total.t.less <- bind_rows(t.less.CLF_CLS, t.less.Sphere_CLF, t.less.Sphere_CLS)

# group_without_filter 
total_without_filter <- bind_rows(total.t.greater, total.t.less)


# UP-regulation 
f.t.greater.CLF_CLS <- f.t.greater.CLF_CLS %>% mutate(Case = "P6+fibroblast > P6", Method = "Filter")
f.t.greater.Sphere_CLF <- f.t.greater.Sphere_CLF %>% mutate(Case = "P6-Sphere > P6+fibroblast", Method = "Filter")
f.t.greater.Sphere_CLS <- f.t.greater.Sphere_CLS %>% mutate(Case = "P6-Sphere > p6", Method = "Filter")
total.f.t.greater <- bind_rows(f.t.greater.CLF_CLS, f.t.greater.Sphere_CLF, f.t.greater.Sphere_CLS)


# DOWN-regulation 
f.t.less.CLF_CLS <- f.t.less.CLF_CLS %>% mutate(Case = "P6+fibroblast < P6", Method = "Filter")
f.t.less.Sphere_CLF <- f.t.less.Sphere_CLF %>% mutate(Case = "P6-Sphere < P6+fibroblast", Method = "Filter") 
f.t.less.Sphere_CLS <- f.t.less.Sphere_CLS %>% mutate(Case = "P6-Sphere < p6", Method = "Filter") 
total.f.t.less <- bind_rows(f.t.less.CLF_CLS, f.t.less.Sphere_CLF, f.t.less.Sphere_CLS)



# group_with_filter 
total_with_filter    <- bind_rows(total.f.t.greater, total.f.t.less)


total_ttest_result <- bind_rows(total_without_filter, total_with_filter)




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
          title = "Result",
          fluidRow(
              selectInput(inputId = "Result.filter",
                          label = "Filter:",
                          choices = c("Nofilter" = "Nofilter","Filter"="Filter"))
          ),
          fluidRow()
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
