library(shiny)
library(gProfileR)
library(plotly)
library(PANTHER.db)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(annotate)
library(hgu133plus2.db)
library(readr)
library(latticeExtra)
library(purrr)
library(rprojroot)
library(markdown)
library(tidyr)
library(httr)
library(rjson)
library(DT)
library(VennDiagram)
source("00_microarry_function.R")
load("total_ttest_result.Rdata")
load("hgu133plus2.RData")
load("total_probe_dataframe.Rdata")
load("annotated_entrez_symbol.Rdata")
load("hgu133plus2_annotation.RData")
load("Exprs_data.RData")
#pthOrganisms(PANTHER.db) <- "HUMAN"


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
              column(5,sliderInput(inputId = "Result.M.upper",
                                   label = "Filter with fold change(M) post ttest: Larger than",
                                   min = -10, max = 12, step = 0.5, value = 0)),
              column(5,sliderInput(inputId = "Result.M.lower",
                                   label = "Filter with fold change(M) post ttest: Lower than",
                                   min = -10, max = 12, step = 0.5, value = 0))
          ),
          fluidRow(
              column(5),
              column(1, actionButton(inputId = "Result.go",label = "Calculate!")),
              column(2,helpText('調整完條件後，按Calculate鈕！'))
          ),
          fluidRow(
              column(1),
              column(11,dataTableOutput(outputId = "Result.table"))),
          # fluidRow(
          #     column(1),
          #     column(11,plotOutput(outputId = "Result.plot"))),
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
                                          choices = c("0.001" = 0.001,
                                                      "0.01"  = 0.01,
                                                      "0.05"  = 0.05
                                                      ))),
                     column(4,sliderInput(inputId = "Gene.M.upper",
                                          label = "Filter with fold change(M) post ttest: Larger than",
                                          min = -10, max = 12, step = 0.5, value = 0)),
                     column(4,sliderInput(inputId = "Gene.M.lower",
                                          label = "Filter with fold change(M) post ttest: Lower than",
                                          min = -10, max = 12, step = 0.5, value = 0))
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
            fluidRow(
                column(5),
                column(1, actionButton(inputId = "Gene.go",label = "Calculate!")),
                column(2,helpText('調整完條件後，按Calculate鈕！'))
            ),
            fluidRow(column(12)),
            fluidRow(
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
                column(1),
                column(10,dataTableOutput(outputId = "Gene.gprofiler")),
                column(1)
            ),
            fluidRow(
                column(1),
                column(5,dataTableOutput(outputId = "Gene.panther")),
                column(5,dataTableOutput(outputId = "Gene.pathway")),
                column(1)
            ),
            fluidRow(
                dataTableOutput(outputId = "Gene.query")
               )
            
            ),

# Page:Fibroblast related ---------------------------------------------------


        tabPanel(title = "Fibroblast caused DE gene related to stemness",
                 fluidRow(
                     column(1),
                     column(3,selectInput(inputId = "Final.filter",
                                          label = "Pre-filter before the t-test:",
                                          choices = c("Nofilter" = "Nofilter","Filter by variance at 50%"="Filter"))),
                     column(4,sliderInput(inputId = "Final.A.upper",
                                          label = "Filter with expression level(A) post ttest: Lower than",
                                          min = 1, max = 15, step = 0.5, value = 15)),
                     column(4,sliderInput(inputId = "Final.A.lower",
                                          label = "Filter with expression level(A) post ttest: Larger than",
                                          min = 0, max = 15, step = 0.5, value = 0))
                 ),
                 fluidRow(
                     column(1),
                     column(3,selectInput(inputId = "Final.p",
                                          label = "Filter with p value",
                                          choices = c("0.001" = 0.001,
                                                      "0.01"  = 0.01,
                                                      "0.05"  = 0.05
                                          ))),
                     column(4,sliderInput(inputId = "Final.M.upper",
                                          label = "Filter with fold change(M) post ttest: Larger than",
                                          min = -10, max = 12, step = 0.5, value = 0)),
                     column(4,sliderInput(inputId = "Final.M.lower",
                                          label = "Filter with fold change(M) post ttest: Lower than",
                                          min = -10, max = 12, step = 0.5, value = 0))
                 ),
                 fluidRow(
                     column(5),
                     column(1, actionButton(inputId = "Final.go",label = "Calculate!")),
                     column(2,helpText('調整完條件後，按Calculate鈕！'))
                 ),
                 fluidRow(
                     column(1),
                     column(5, helpText('Both Up regulation gene in P6+fibroblast and P6-sphere compare to P6 ')),
                     column(5,helpText('Both Down regulation gene in P6+fibroblast and P6-sphere compare to P6')),
                     column(1)
                 ),
                fluidRow(
                column(1),
                column(5,plotOutput(outputId = "Final.venn.up")),
                column(5,plotOutput(outputId = "Final.venn.down")),
                column(1)
                
                ),
                fluidRow(
                    column(1),
                    column(5,dataTableOutput(outputId = "Final.up.result")),
                    column(5,dataTableOutput(outputId = "Final.down.result")),
                    column(1)
                    
                ),
                fluidRow(
                    plotlyOutput(outputId = "Final.3d")
                )
            )
        ) 
    )






# SERVER part -------------------------------------------------------------


server <- function(input, output){
# Output$Result.table -----------------------------------------------------
    output$Result.table <- renderDataTable({
        # tabke a dependency on action button
        input$Result.go
        # isolate to avoid dependency
        isolate(total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2),
                                      M = estimate1 - estimate2) %>% group_by(Case) %>% 
            filter(Method == input$Result.filter) %>%
            filter(!Case %in% c ("P6-Sphere > P6+fibroblast",
                                 "P6-Sphere < P6+fibroblast",
                                 "P6-Sphere < P6+fibroblast",
                                 "P6-Sphere > P6+fibroblast")) %>% 
            filter(A < input$Result.A.upper) %>%
            filter(A > input$Result.A.lower) %>%
            filter(M > input$Result.M.upper | M < input$Result.M.lower) %>%
            summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                      n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                      n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                      n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE)))
    })
# # Output$Result.plot ------------------------------------------------------ 
#     output$Result.plot <- renderPlot({
#         total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2)) %>% ggplot() +
#             geom_violin(aes(x = Case, y = A)) +
#             facet_grid(Method ~ .) +
#             geom_hline(yintercept = input$Result.A.upper) +
#             geom_hline(yintercept = input$Result.A.lower) 
#             
#         
#     })

# Output$Result.MAplot ----------------------------------------------------
    output$Result.MAplot <- renderPlot({
        # tabke a dependency on action button
        input$Result.go
        # isolate to avoid dependency
        
        isolate(total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2),
                                      M = estimate1 - estimate2,
                                      P = cut(adjusted.p, c(0,0.0001,0.001,0.01,0.05,1),c("p<0.0001","p<0.001","p<0.01","p<0.05","p>0.05"))) %>% ggplot() +
            geom_point(aes(x = A, y = M, colour=P), alpha = 0.5) +
            facet_grid(Method ~ Case) +
            geom_vline(xintercept = input$Result.A.upper) +
            geom_vline(xintercept = input$Result.A.lower) +
            geom_hline(yintercept = input$Result.M.upper) +
            geom_hline(yintercept = input$Result.M.lower))
            
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
        input$Gene.go
        isolate(total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2),
                                      M = estimate1 - estimate2) %>%  
            filter(Method == input$Gene.filter) %>%
            filter(Case == input$Gene.case) %>%
            filter(adjusted.p < input$Gene.p) %>%
            filter(A < input$Gene.A.upper) %>%
            filter(A > input$Gene.A.lower) %>%
            filter(M > input$Gene.M.upper | M < input$Gene.M.lower) )
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
    
# Output$Gene.profiler ----------------------------------------------------    
    output$Gene.gprofiler <- renderDataTable({
        
        gprofiler(query = probe.list(),
                  organism = "hsapiens",
                  ordered_query = T,
                  max_set_size = 1000,
                  min_isect_size = 0,
                  significant = T,
                  exclude_iea = T,
                  underrep = T,
                  correction_method = "fdr",
                  src_filter = c("GO:BP","REAC"),
                  hier_filtering = "moderate",
                  custom_bg = hgu133plus2.probe)
    })    
    
    
# Output$Gene.panther -----------------------------------------------------
    panther.result <- reactive({
        
        panther.probe <-probe.list()
        tmp <- annotated.entrez.symbol %>% filter(Probe %in% panther.probe)
        tmp.EntrezID <- tmp$EntrezID[!is.na(tmp$EntrezID)]
        k <- tmp.EntrezID
        choice <- c("CLASS_TERM")
        PANTHER.db::select(PANTHER.db,
               keys = k,
               columns = choice,
               keytype = "ENTREZ")   
    })
      
    panther.pathway <- reactive({
        
        panther.probe <-probe.list()
        tmp <- annotated.entrez.symbol %>% filter(Probe %in% panther.probe)
        tmp.EntrezID <- tmp$EntrezID[!is.na(tmp$EntrezID)]
        k <- tmp.EntrezID
        choice <- c("PATHWAY_TERM")
        PANTHER.db::select(PANTHER.db,
                           keys = k,
                           columns = choice,
                           keytype = "ENTREZ")   
    })
    
     output$Gene.panther <- renderDataTable({
         
         #panther.result() %>% filter(!is.na(CLASS_TERM)) %>% ggplot + geom_bar(aes(x = CLASS_TERM)) + coord_polar(theta = "x", direction=1)
         panther.result()$CLASS_TERM %>% table %>% sort() %>% as.data.frame %>% arrange(desc(Freq))
     })
     # Output$Gene.pathway -----------------------------------------------------    
     output$Gene.pathway <- renderDataTable({
         
         #panther.result() %>% filter(!is.na(CLASS_TERM)) %>% ggplot + geom_bar(aes(x = CLASS_TERM)) + coord_polar(theta = "x", direction=1)
         panther.pathway()$PATHWAY_TERM %>% table %>% sort() %>% as.data.frame %>% arrange(desc(Freq))
     })
     
     # Reactive for Final data -----------------------------------------------------------
     
     
     tmp.final <- reactive({
         input$Final.go
         isolate(total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2),
                                               M = estimate1 - estimate2) %>%  
                     filter(Method == input$Final.filter) %>%
                     filter(adjusted.p < input$Final.p) %>%
                     filter(A < input$Final.A.upper) %>%
                     filter(A > input$Final.A.lower) %>%
                     filter(M > input$Final.M.upper | M < input$Final.M.lower) )
     })
     
     final.gene.list <- reactive({
         unname(getSYMBOL(tmp.Final()$gene,"hgu133plus2.db"))  
     })
     
     final.probe.list <- reactive({
         names(getSYMBOL(tmp.Final()$gene,"hgu133plus2.db"))
     })
  
     up.CLF_CLS    <- reactive({
         names(getSYMBOL((tmp.final() %>% filter(Case == "P6+fibroblast > P6"))$gene,"hgu133plus2.db"))
     })
         
     up.Sphere_CLS <- reactive({
         names(getSYMBOL((tmp.final() %>% filter(Case == "P6-Sphere > p6"))$gene,"hgu133plus2.db"))
     }) 
     
     down.CLF_CLS    <- reactive({
         names(getSYMBOL((tmp.final() %>% filter(Case == "P6+fibroblast < P6"))$gene,"hgu133plus2.db"))
     })
     
     down.Sphere_CLS <- reactive({
         names(getSYMBOL((tmp.final() %>% filter(Case == "P6-Sphere < p6"))$gene,"hgu133plus2.db"))
     }) 
     # Output$Final.venn.up ----------------------------------------------------
       
     output$Final.venn.up <- renderPlot({
         area.1 <- length(up.CLF_CLS() )
         area.2 <- length(up.Sphere_CLS())
         area.12 <- length(intersect(up.CLF_CLS(), up.Sphere_CLS()))
         area.list <- list(area.1, area.2, area.12)
         draw.pairwise.venn(area1 = area.list[[1]], area2 = area.list[[2]], cross.area = area.list[[3]],
                            category = c("P6+fibroblast > P6", "P6-Sphere > p6"),
                            fill = c("Blue","Red"))
     })
       
         
     # Output$Final.venn.down --------------------------------------------------
     output$Final.venn.down <- renderPlot({
         area.1 <- length(down.CLF_CLS() )
         area.2 <- length(down.Sphere_CLS())
         area.12 <- length(intersect(down.CLF_CLS(), down.Sphere_CLS()))
         area.list <- list(area.1, area.2, area.12)
         draw.pairwise.venn(area1 = area.list[[1]], area2 = area.list[[2]], cross.area = area.list[[3]],
                            category = c("P6+fibroblast < P6", "P6-Sphere < p6"),
                            fill = c("Blue", "Red"))
     })
     
     # Output$Final.up.result --------------------------------------------------
     
     output$Final.up.result <- renderDataTable({
         
         left_join(annotated.entrez.symbol%>%select(-Symbol), hgu133plus2.probe.annotate[,c("Probe", "name", "description")]) %>%
         filter(Probe %in% intersect(up.CLF_CLS(), up.Sphere_CLS()))
     })
     
     
     
     # Output$Final.down.result ------------------------------------------------
     
     output$Final.down.result <- renderDataTable({
         
         left_join(annotated.entrez.symbol%>%select(-Symbol), hgu133plus2.probe.annotate[,c("Probe", "name", "description")]) %>%
             filter(Probe %in% intersect(down.CLF_CLS(), down.Sphere_CLS()))
     })
     
     # Output$Final.3d ---------------------------------------------------------
     exprssion <- log2(exprs(Exprs.data)+1) %>% as.data.frame()
     exprssion$Probe <- rownames(exprssion)
     scatter <- exprssion %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                                CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                                Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) %>%
                dplyr::select(Probe, CLF, CLS, Sphere)
     scatter <- left_join(scatter, hgu133plus2.probe.annotate[,c("Probe", "name")])
     usr.col <- brewer.pal(3, "Set1")
     
     scatter.data <- reactive({
         scatter$tag[scatter$Probe %in% up.CLF_CLS()] <- "P6+fibroblast > P6.Up"
         scatter$tag[scatter$Probe %in% up.Sphere_CLS()] <- "P6-Sphere > p6.Up"
         scatter$tag[scatter$Probe %in% down.CLF_CLS()] <- "P6+fibroblast < P6.Down"
         scatter$tag[scatter$Probe %in% down.Sphere_CLS()] <- "P6-Sphere < p6.Down"
         scatter$tag[scatter$Probe %in% intersect(up.CLF_CLS(), up.Sphere_CLS())] <- "Intersection.Up"
         scatter$tag[scatter$Probe %in% intersect(down.CLF_CLS(), down.Sphere_CLS())] <- "Intersection.Down"
         scatter$tag[is.na(scatter$tag)] <- "not"
         scatter
     })
     
     output$Final.3d <- renderPlotly({
         
         
         plot_ly(scatter.data(), x = ~CLF, y = ~CLS, z = ~Sphere,
                 color = ~tag, colors = usr.col, text= ~paste('Symbol:',name,'<br>Probe:',Probe)) %>%
             add_markers() %>%
             layout(scene = list(xaxis = list(title = 'CLF'),
                                 yaxis = list(title = 'CLS'),
                                 zaxis = list(title = 'Sphere')))
         
     })
     
     
     
}















 










shinyApp(ui = ui, server = server)
