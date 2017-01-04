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
# Page:Project ------------------------------------------------------------
    navbarPage(title = "Project",
        tabPanel(
            title = "Introduction",
            column(12,
                   includeMarkdown("DESCRIPTION.md")
                   )
          ),

# Page:Workflow -----------------------------------------------------------
        tabPanel(
            title = "Workflow",
            img(src="workflow.png")
           ),

# Page:Result Overview ----------------------------------------------------
        tabPanel(
          title = "Result Overview",
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
                                    min = 1, max = 15, step = 0.5, value = 15)),
              column(5,sliderInput(inputId = "Result.A.lower",
                                    label = "Filter with expression level(A) post ttest: Larger than",
                                    min = 0, max = 15, step = 0.5, value = 0))
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

# Page:DE Gene Overview ---------------------------------------------------------------


        tabPanel(
            title = "DE Gene Overview",
            fluidRow(
                column(1),
                column(3,selectInput(inputId = "Gene.filter",
                                      label = "Pre-filter before the t-test:",
                                      choices = c("Nofilter" = "Nofilter","Filter by variance at 50%"="Filter"))),
                column(4,sliderInput(inputId = "Gene.A.upper",
                                            label = "Filter with expression level(A) post ttest: Lower than",
                                            min = 1, max = 15, step = 0.5, value = 15)),
                column(4,sliderInput(inputId = "Gene.A.lower",
                                            label = "Filter with expression level(A) post ttest: Larger than",
                                            min = 0, max = 15, step = 0.5, value = 0))
               ),
            fluidRow(
                     column(1),
                     column(3,selectInput(inputId = "Gene.p",
                                          label = "Filter with p value",
                                          choices = c("0.05"  = 0.05,
                                                      "0.01"  = 0.01,
                                                      "0.001" = 0.001))),
                     column(4,sliderInput(inputId = "Gene.M.upper",
                                          label = "Filter with fold change(M) post ttest: Larger than",
                                          min = -10, max = 12, step = 0.5, value = 1)),
                     column(4,sliderInput(inputId = "Gene.M.lower",
                                          label = "Filter with fold change(M) post ttest: Lower than",
                                          min = -10, max = 12, step = 0.5, value = -1))
                ),
            fluidRow(
                column(1),
                column(6,selectInput(inputId = "Gene.case",
                                         label = "Multiple choose the inspect case",
                                         choices = c("P6-Sphere > p6",
                                                     "P6-Sphere < p6",
                                                     "P6+fibroblast > P6",
                                                     "P6+fibroblast < P6")
                                        ))
               ),
            fluidRow(column(12)),
            fluidRow(
                #column(3,helpText('The Number of DE probe')),
                column(1),
                column(5, plotOutput(outputId = "Gene.MAplot")),
                column(4, dataTableOutput('Gene.number')),
                column(2)
               ),
            fluidRow(
                  verbatimTextOutput(outputId = "Gene.select")
                 #column(3, plotOutput(outputId = "Gene.MAplot"))
                 #column(9,dataTableOutput(outputId = "Gene.query"))
                ),
            fluidRow(
                dataTableOutput(outputId = "Gene.query")
               )
            ),
        tabPanel(title = "Fibroblast caused DE gene related to stemness",
            fluidRow())
        ) 
    )






# SERVER part -------------------------------------------------------------


server <- function(input, output){
# Output$Result.table -----------------------------------------------------
    output$Result.table <- renderDataTable({
        total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2)) %>% group_by(Case) %>% 
            filter(Method == input$Result.filter) %>%
            filter(!Case %in% c ("P6-Sphere > P6+fibroblast",
                                 "P6-Sphere < P6+fibroblast",
                                 "P6-Sphere < P6+fibroblast",
                                 "P6-Sphere > P6+fibroblast")) %>% 
            filter(A < input$Result.A.upper) %>%
            filter(A > input$Result.A.lower) %>%
            summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                      n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                      n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                      n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE))
    })
# Output$Result.plot ------------------------------------------------------ 
    output$Result.plot <- renderPlot({
        total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2)) %>% ggplot() +
            geom_violin(aes(x = Case, y = A)) +
            facet_grid(Method ~ .) +
            geom_hline(yintercept = input$Result.A.upper) +
            geom_hline(yintercept = input$Result.A.lower) 
            
        
    })

# Output$Result.MAplot ----------------------------------------------------
    output$Result.MAplot <- renderPlot({
        total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2),
                                      M = estimate1 - estimate2,
                                      P = cut(adjusted.p, c(0,0.0001,0.001,0.01,0.05,1),c("p<0.0001","p<0.001","p<0.01","p<0.05","p>0.05"))) %>% ggplot() +
            geom_point(aes(x = A, y = M, colour=P), alpha = 0.5) +
            facet_grid(Method ~ Case) +
            geom_vline(xintercept = input$Result.A.upper) +
            geom_vline(xintercept = input$Result.A.lower)
    })


# Output$Gene.MAplot ------------------------------------------------------
    output$Gene.MAplot <- renderPlot({
        total_ttest_result %>%
            filter(Case %in% input$Gene.case) %>%
            filter(Method == input$Gene.filter) %>% 
            mutate(A = 0.5*(estimate1 + estimate2),
                                      M = estimate1 - estimate2,
                                      P = cut(adjusted.p, c(0,0.0001,0.001,0.01,0.05,1),c("p<0.0001","p<0.001","p<0.01","p<0.05","p>0.05"))) %>% ggplot() +
            geom_point(aes(x = A, y = M, colour=P), alpha = 0.5) +
            facet_grid(Method ~ Case) +
            geom_vline(xintercept = input$Gene.A.upper) +
            geom_vline(xintercept = input$Gene.A.lower) +
            geom_hline(yintercept = input$Gene.M.upper) +
            geom_hline(yintercept = input$Gene.M.lower) 
    })

# Reactive data -----------------------------------------------------------


    tmp <- reactive({
        total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2),
                                      M = estimate1 - estimate2) %>%  
            filter(Method == input$Gene.filter) %>%
            filter(Case == input$Gene.case) %>%
            filter(adjusted.p < input$Gene.p) %>%
            filter(A < input$Gene.A.upper) %>%
            filter(A > input$Gene.A.lower) %>%
            filter(M > input$Gene.M.upper | M < input$Gene.M.lower) 
    })
    
    gene.list <- reactive({
        unname(getSYMBOL(tmp()$gene,"hgu133plus2.db"))  
    })
    
    probe.list <- reactive({
        names(getSYMBOL(tmp()$gene,"hgu133plus2.db"))
    })

# Output$Gene.number ------------------------------------------------------  
    output$Gene.number <- renderDataTable({
        DT::datatable(data.frame(Probe = probe.list(),
                   Symbol = gene.list(),
                   Adjusted.p = tmp()$adjusted.p),
                   options = list(
                       lengthMenu = list(c(3, 15, -1), c("3", "15", "All")),
                       pageLength = 15
                   ))
       
    })
# Output$Gene.select ------------------------------------------------------
   output$Gene.select <- renderPrint({
      s <- input$Gene.number_rows_selected
      if (length(s)){
          cat('Thess rows were selected:\n\n')
          cat(s, sep=',')
      }
   })
        

# Output$Gene.query -------------------------------------------------------
    output$Gene.query <-  renderDataTable({
        s <- input$Gene.number_rows_selected
        if (length(s)){
            searchHarmonizome(gene.list()[s])
            }
        
    })

    }















shinyApp(ui = ui, server = server)
