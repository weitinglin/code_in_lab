library(shiny)
source("/Users/Weitinglin/Documents/Repository/code_in_lab/00_microarry_function.R")
load("/Users/Weitinglin/Documents/Repository/code_in_lab/total_ttest_result.Rdata")
load("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/data/intermediate/hgu133plus2.RData")
load("/Users/Weitinglin/Documents/Repository/code_in_lab/total_probe_dataframe.Rdata")


# Preprocess --------------------------------------------------------------

# # UP-regulation 
# 
# 
# t.greater.CLF_CLS <- t.greater.CLF_CLS %>% mutate(Case = "P6+fibroblast > P6", Method = "Nofilter")
# t.greater.Sphere_CLF <- t.greater.Sphere_CLF %>% mutate(Case = "P6-Sphere > P6+fibroblast", Method = "Nofilter")
# t.greater.Sphere_CLS <- t.greater.Sphere_CLS %>% mutate(Case = "P6-Sphere > p6", Method = "Nofilter")
# total.t.greater <- bind_rows(t.greater.CLF_CLS, t.greater.Sphere_CLF, t.greater.Sphere_CLS)
# 
# 
# # DOWN-regulation 
# t.less.CLF_CLS <- t.less.CLF_CLS %>% mutate(Case = "P6+fibroblast < P6", Method = "Nofilter")
# t.less.Sphere_CLF <- t.less.Sphere_CLF %>% mutate(Case = "P6-Sphere < P6+fibroblast", Method = "Nofilter")
# t.less.Sphere_CLS <- t.less.Sphere_CLS %>% mutate(Case = "P6-Sphere < p6", Method = "Nofilter")
# total.t.less <- bind_rows(t.less.CLF_CLS, t.less.Sphere_CLF, t.less.Sphere_CLS)
# 
# # group_without_filter 
# total_without_filter <- bind_rows(total.t.greater, total.t.less)
# 
# 
# # UP-regulation 
# f.t.greater.CLF_CLS <- f.t.greater.CLF_CLS %>% mutate(Case = "P6+fibroblast > P6", Method = "Filter")
# f.t.greater.Sphere_CLF <- f.t.greater.Sphere_CLF %>% mutate(Case = "P6-Sphere > P6+fibroblast", Method = "Filter")
# f.t.greater.Sphere_CLS <- f.t.greater.Sphere_CLS %>% mutate(Case = "P6-Sphere > p6", Method = "Filter")
# total.f.t.greater <- bind_rows(f.t.greater.CLF_CLS, f.t.greater.Sphere_CLF, f.t.greater.Sphere_CLS)
# 
# 
# # DOWN-regulation 
# f.t.less.CLF_CLS <- f.t.less.CLF_CLS %>% mutate(Case = "P6+fibroblast < P6", Method = "Filter")
# f.t.less.Sphere_CLF <- f.t.less.Sphere_CLF %>% mutate(Case = "P6-Sphere < P6+fibroblast", Method = "Filter") 
# f.t.less.Sphere_CLS <- f.t.less.Sphere_CLS %>% mutate(Case = "P6-Sphere < p6", Method = "Filter") 
# total.f.t.less <- bind_rows(f.t.less.CLF_CLS, f.t.less.Sphere_CLF, f.t.less.Sphere_CLS)
# 
# 
# 
# # group_with_filter 
# total_with_filter    <- bind_rows(total.f.t.greater, total.f.t.less)
# 
# 
# total_ttest_result <- bind_rows(total_without_filter, total_with_filter)




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
              column(1),
              column(11,selectInput(inputId = "Result.filter",
                          label = "Pre-filter before the t-test:",
                          choices = c("Nofilter" = "Nofilter","Filter by variance at 50%"="Filter")))
          ),
          fluidRow(
              column(1),
              column(5,sliderInput(inputId = "Result.A.upper",
                                    label = "Filter with expression level(A) post ttest: Lower than",
                                    min = 1, max = 12, step = 0.5, value = 12)),
              column(5,sliderInput(inputId = "Result.A.lower",
                                    label = "Filter with expression level(A) post ttest: Larger than",
                                    min = 1, max = 12, step = 0.5, value = 1))
          ),
          fluidRow(
              column(1),
              column(11,dataTableOutput(outputId = "Result.table"))),
          fluidRow(
              column(1),
              column(11,plotOutput(outputId = "Result.plot"))),
          fluidRow(
              column(1),
              column(11,plotOutput(outputId = "Result.MAplot")))
        ),
        tabPanel(
            title = "Gene",
            fluidRow(
                column(1),
                column(3,selectInput(inputId = "Gene.filter",
                                      label = "Pre-filter before the t-test:",
                                      choices = c("Nofilter" = "Nofilter","Filter by variance at 50%"="Filter")),
                column(4,sliderInput(inputId = "Gene.A.upper",
                                            label = "Filter with expression level(A) post ttest: Lower than",
                                            min = 1, max = 12, step = 0.5, value = 12)),
                column(4,sliderInput(inputId = "Gene.A.lower",
                                            label = "Filter with expression level(A) post ttest: Larger than",
                                            min = 1, max = 12, step = 0.5, value = 1))
            ),
            fluidRow(dataTableOutput(outputId = "queryResult")))
        ) 
    )
)


# SERVER part -------------------------------------------------------------


server <- function(input, output){
    
    output$Result.table <- renderDataTable({
        total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2)) %>% group_by(Case) %>% 
            filter(Method == input$Result.filter) %>%
            filter(!Case %in% c ("P6-Sphere > P6+fibroblast","P6-Sphere < P6+fibroblast", "P6-Sphere < P6+fibroblast nofilter", "P6-Sphere > P6+fibroblast nofilter")) %>% 
            filter(A < input$Result.A.upper) %>%
            filter(A > input$Result.A.lower) %>%
            summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                      n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                      n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                      n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE))
    })
    
    output$Result.plot <- renderPlot({
        total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2)) %>% ggplot() +
            geom_violin(aes(x = Case, y = A)) +
            facet_grid(Method ~ .) +
            geom_hline(yintercept = input$Result.A.upper) +
            geom_hline(yintercept = input$Result.A.lower) 
            
        
    })
    
    output$Result.MAplot <- renderPlot({
        total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2),
                                      M = estimate1 - estimate2,
                                      P = cut(adjusted.p, c(0,0.0001,0.001,0.01,0.05,1),c("p<0.0001","p<0.001","p<0.01","p<0.05","p>0.05"))) %>% ggplot() +
            geom_point(aes(x = A, y = M, colour=P), alpha = 0.5) +
            facet_grid(Method ~ Case) +
            geom_vline(xintercept = input$Result.A.upper) +
            geom_vline(xintercept = input$Result.A.lower)
    })
    
    output$queryResult <-  renderDataTable({
        searchHarmonizome(c("CD44","ALDH1A3","CD9","CDKN2A","DPP4","HSPB1","KIT","NANOG"))
    })
    
    
    
    
    }



shinyApp(ui = ui, server = server)
