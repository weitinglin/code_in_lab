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
load("updown_tukey_aov_result.Rdata")
load("norm_qn_log.Rdata")
load("hgu133plus2_annotation.RData")

# UI ----------------------------------------------------------------------


ui <- tagList(
  navbarPage(title = "Project",
             # _Introduction ---------------------------------------------------------
             tabPanel(
               title = "Introduction",
               column(12, includeMarkdown("DESCRIPTION.md"))
             ),
             # _Overview ---------------------------------------------------------
             tabPanel(
               title = "Overview",
               fluidRow(
                 column(1),
                 column(11,dataTableOutput(outputId = 'Overview.table'))
                 ),
               fluidRow(
                   column(6,plotlyOutput(outputId = 'Overview.MA.CLF_CLS')),
                   column(6,plotlyOutput(outputId = 'Overview.MA.Sphere_CLS'))
               )
              ),
              tabPanel(
                title = "DEgene",
                fluidRow(
                    column(2),
                    column(10,selectInput(inputId = "DEgene.cutoff",
                                     label = " p value cutoff",
                                     choices = c("0.0001" = 1,
                                                 "0.001"  = 2,
                                                 "0.05"  = 3,
                                                 "0.00001" = 4
                                     )))),
                fluidRow(column(6,plotOutput(outputId = 'DEgene.venn.up')),
                         column(6,plotOutput(outputId = 'DEgene.venn.down'))),
                fluidRow(column(1),
                         column(4, dataTableOutput(outputId = 'DEgene.up.table')),
                         column(2),
                         column(4, dataTableOutput(outputId = 'DEgene.down.table'))),
                         column(1)
                
              )
             )
)





# Server ------------------------------------------------------------------
tmp.up   <- bind_rows(up.tukey.aov.CLF_CLS, up.tukey.aov.Sphere_CLS) %>% mutate(type = "up")

tmp.down <- bind_rows(down.tukey.aov.CLF_CLS, down.tukey.aov.Sphere_CLS) %>% mutate(type = "down")

total.aov.tukey <- bind_rows(tmp.up, tmp.down)

#Norm
MA.aov.data <- norm %>% as.data.frame %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                                                 CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                                                 Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) 
MA.aov.data$Probe <- rownames(norm)
MA.aov.data <- MA.aov.data %>% dplyr::select(CLF,CLS,Sphere,Probe)

# MA plot

MA.plot.CLS_CLF <- total.aov.tukey %>% filter(case == "CLS-CLF")
MA.plot.CLS_CLF <- left_join(MA.plot.CLS_CLF, MA.aov.data[,c("CLS","CLF","Probe")], by="Probe")
MA.plot.CLS_CLF <- MA.plot.CLS_CLF %>% mutate(M = CLF-CLS, A = 0.5*(CLF+CLS),
                                              P= cut(BH.p.adj, c(0,0.0001,0.001,0.01,0.05,1),
                                                     c("p<0.0001","p<0.001","p<0.01","p<0.05","p>0.05")))

MA.plot.Sphere_CLS <- total.aov.tukey %>% filter(case == "Sphere-CLS")
MA.plot.Sphere_CLS <- left_join(MA.plot.Sphere_CLS, MA.aov.data[,c("Sphere","CLS","Probe")], by="Probe")
MA.plot.Sphere_CLS <- MA.plot.Sphere_CLS %>% mutate(M = Sphere-CLS, A=0.5*(Sphere + CLS),
                                                    P = cut(BH.p.adj, c(0,0.0001,0.001,0.01,0.05,1),
                                                            c("p<0.0001","p<0.001","p<0.01","p<0.05","p>0.05")))

usr.col <- brewer.pal(3, "Set1")

# _Preprocess --------------------------------------------------------------




server<- function(input, output){
  # _Overview.table--------------------------------------------------------------
  output$Overview.table <- renderDataTable({
    total.aov.tukey %>%
      group_by(type, case) %>%
      summarise(n_005 = sum(BH.p.adj < 0.05),
                n_001 = sum(BH.p.adj < 0.01),
                n_0001 = sum( BH.p.adj < 0.0001),
                n_00001 = sum( BH.p.adj < 0.00001))
    
  })
  # _Overview.MA--------------------------------------------------------------
  output$Overview.MA.CLF_CLS <- renderPlotly({
      plot_ly(MA.plot.CLS_CLF, x = ~A, y = ~M,
              color = ~P, colors = usr.col, text= ~paste('Symbol:', Symbol ,'<br>Probe:', Probe, '<br>Type:', type)) %>%
          add_markers() %>%
          layout(title = "P6-P6+fibroblast",scene = list(xaxis = list(title = 'A'),
                              yaxis = list(title = 'M')))
  })
  output$Overview.MA.Sphere_CLS <- renderPlotly({
      plot_ly(MA.plot.Sphere_CLS, x = ~A, y = ~M,
              color = ~P, colors = usr.col, text= ~paste('Symbol:', Symbol, '<br>Probe:', Probe, '<br>Type:', type)) %>%
          add_markers() %>%
          layout(title = "Sphere.p6-P6",scene = list(xaxis = list(title = 'A'),
                              yaxis = list(title = 'M')))
  })
  
  
  
  
  # _DEgene--------------------------------------------------------------
  
  input.cutoff.p <- reactive({
      
      if ( input$DEgene.cutoff == 1){
          return(0.0001)
      }else if(input$DEgene.cutoff == 2){
          return(0.001)
      }else if(input$DEgene.cutoff == 3){
          return(0.05)
      }else if(input$DEgene.cutoff == 4){
          return(0.00001)
      }
      
  })
  # _DEgene.venn.up--------------------------------------------------------------
  output$DEgene.venn.up <- renderPlot({
      cutoff.p <- input.cutoff.p()
      up.CLF_CLS    <- up.tukey.aov.CLF_CLS %>% filter(BH.p.adj < cutoff.p)
      up.Sphere_CLS <- up.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)
      
      area.1 <- length(up.CLF_CLS$Probe )
      area.2 <- length(up.Sphere_CLS$Probe)
      area.12 <- length(intersect(up.CLF_CLS$Probe, up.Sphere_CLS$Probe))
      area.list <- list(area.1, area.2, area.12)
      draw.pairwise.venn(area1 = area.list[[1]], area2 = area.list[[2]], cross.area = area.list[[3]],
                         category = c("P6+fibroblast > P6", "P6-Sphere > p6"),
                         fill = c("Blue","Red"), top="Up Regulation")
      
  })
  # _DEgene.venn.down--------------------------------------------------------------
  output$DEgene.venn.down <- renderPlot({
      cutoff.p <- input.cutoff.p()
      down.CLF_CLS    <- down.tukey.aov.CLF_CLS  %>% filter(BH.p.adj < cutoff.p)
      down.Sphere_CLS <- down.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)
      # The DE gene related to the 
      area.1 <- length(down.CLF_CLS$Probe )
      area.2 <- length(down.Sphere_CLS$Probe )
      area.12 <- length(intersect(down.CLF_CLS$Probe, down.Sphere_CLS$Probe ))
      area.list <- list(area.1, area.2, area.12)
      draw.pairwise.venn(area1 = area.list[[1]], area2 = area.list[[2]], cross.area = area.list[[3]],
                         category = c("P6+fibroblast < P6", "P6-Sphere < p6"),
                         fill = c("Blue", "Red"), top="Down Regulation")
  })
  
  
  # _DEgene.up.table--------------------------------------------------------------
  
  output$DEgene.up.table <- renderDataTable({
      cutoff.p <- input.cutoff.p()
      up.CLF_CLS    <- up.tukey.aov.CLF_CLS %>% filter(BH.p.adj < cutoff.p) %>% mutate(CLF_CLS.BH.p.adj=BH.p.adj)
      up.Sphere_CLS <- up.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p) %>% mutate(Sphere_CLS.BH.p.adj=BH.p.adj)
      up.probe <- intersect(up.CLF_CLS$Probe, up.Sphere_CLS$Probe )
      up.CLF_CLS <- up.CLF_CLS %>% filter(Probe %in% up.probe) %>% dplyr::select(Probe, CLF_CLS.BH.p.adj)
      up.Sphere_CLS <- up.Sphere_CLS %>% filter(Probe %in% up.probe) %>% dplyr::select(Probe, Sphere_CLS.BH.p.adj)
      up.merge <- left_join(up.CLF_CLS, up.Sphere_CLS, by = "Probe")
      tmp <- hgu133plus2.probe.annotate %>% dplyr::select(-alias.number,-alias,-target.number,-target,-namespace)%>% filter(Probe %in% up.probe)
      left_join(tmp, up.merge, by="Probe")
      
  })
  
  
  # _DEgene.down.table--------------------------------------------------------------
  output$DEgene.down.table <- renderDataTable({
      cutoff.p <- input.cutoff.p()
      down.CLF_CLS    <- down.tukey.aov.CLF_CLS  %>% filter(BH.p.adj < cutoff.p) %>% mutate(CLF_CLS.BH.p.adj=BH.p.adj)
      down.Sphere_CLS <- down.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p) %>% mutate(Sphere_CLS.BH.p.adj=BH.p.adj)
      down.probe <- intersect(down.CLF_CLS$Probe, down.Sphere_CLS$Probe )
      down.CLF_CLS <- down.CLF_CLS %>% filter(Probe %in% down.probe) %>% dplyr::select(Probe, CLF_CLS.BH.p.adj) 
      down.Sphere_CLS <- down.Sphere_CLS %>% filter(Probe %in% down.probe) %>% dplyr::select(Probe, Sphere_CLS.BH.p.adj)
      down.merge <- left_join(down.CLF_CLS, down.Sphere_CLS)
      tmp <- hgu133plus2.probe.annotate %>% dplyr::select(-alias.number,-alias,-target.number,-target,-namespace) %>%filter(Probe %in% down.probe)
      left_join(tmp, down.merge, by="Probe")
  })
  
  
  
  
}

shinyApp(ui = ui, server = server)