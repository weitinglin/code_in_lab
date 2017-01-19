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
load("test_annova_result.Rdata")
load("updown_tukey_aov_result.Rdata")
load("norm_qn_log.Rdata")


# UI ----------------------------------------------------------------------


ui <- tagList(
  navbarPage(title = "Project",
             # _Introduction ---------------------------------------------------------
             tabPanel(
               title = "Introduction",
               column(12, includeMarkdown("DESCRIPTION.md")),
               fluidRow()
             ),
             # _Overview ---------------------------------------------------------
             tabPanel(
               title = "Overview",
               fluidRow(
                 column(1),
                 column(4, dataTableOutput(outputId = 'Overview.table')),
                 column(3, plotOutput(outputId = 'Overview.MA.CLF_CLS')),
                 column(3, plotOutput(outputId = 'Overview.MA.Sphere_CLS'))
               )
             ))
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

MA.plot.CLS_CLF <- total.aov.tukey %>% filter(cas == "CLS-CLF")
MA.plot.CLS_CLF <- left_join(MA.plot.CLS_CLF, MA.aov.data[,c("CLS","CLF","Probe")], by="Probe")
MA.plot.CLS_CLF <- MA.plot.CLS_CLF %>% mutate(M=CLF-CLS, A=0.5(CLF+CLS), P= cut(adjusted.p, c(0,0.0001,0.001,0.01,0.05,1),c("p<0.0001","p<0.001","p<0.01","p<0.05","p>0.05")))
ggplot(data=MA.plot.CLS_CLF) + geom_point(aes(x = A, y = M, alpha=0.5))

MA.plot.Sphere_CLS <- total.aov.tukey %>% filter(cas == "Sphere-CLS")
MA.plot.Sphere_CLS <- left_join(MA.plot.Sphere_CLF, MA.aov.data[,c("Sphere","CLF","Probe")], by="Probe")
MA.plot.Sphere_CLS <- MA.plot.Sphere_CLS %>% mutate(M=Sphere-CLS, A=0.5(Sphere + CLS))

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
  output$Overview.MA.CLF_CLS <- renderPlot({
    
  })
  output$Overview.MA.Sphere_CLS <- renderPlot({
    
  })
}

shinyApp(ui = ui, server = server)