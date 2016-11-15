library(shiny)
library(VennDiagram)
source("/Users/Weitinglin/Documents/R_scripts/Lab/microarray201610/microarry_function.R")



ui <- fluidPage(
  titlePanel("Lung Cancer Microarray analysis"),
  
  sidebarLayout(
    
  sidebarPanel(
  #original distribution and post MASS background correction
  radioButtons(inputId = "origin.effect",
               label = "After the MASS, inspected effect by case or by sample",
               c("case"= "case", "sample"="sample")),
    
  #filtering method
  radioButtons(inputId = "prefilter",
               label = "prefiltering method",
               c("Variance Based Filtering"="filterbyvariance", "Filter with expression level"="filterbyexpressionlevel", "No Filter"="nofilter")),
  sliderInput(inputId = "variance_value", label="filtering by variance",min = 0,max=1,value = 0.1),
  numericInput(inputId = "ratio", label="filtering by expression:ratio",min = 0,max=1,value=1/3),
  numericInput(inputId = "exprlevel", label="filtering by expression:expression level",min = 0, max = 10, value=log2(100)),
  #differential expression method
  radioButtons(inputId = "demethod",
               label = "differential expression method",
               c("Wilcox"="wilcox","Limma"="limma")),

  # #p-value choose
  sliderInput(inputId = "pvalue", label ="parameters with the expression method(limma+bayes)", min=0, max=1, value=0.05),
  sliderInput(inputId = "prior", label ="parameters with the expression method(wilcox+bayes)", min = 0, max=1, value=0.2),
  sliderInput(inputId = "sig", label ="parameters with the expression method(wilcox+bayes)", min=0, max=1, value=0.05)
  ),
   mainPanel(
     tabsetPanel(
     id = "plot",
     tabPanel('Density',plotOutput(outputId = 'denDIS')),
     tabPanel('ClustHeatmap',plotOutput(outputId = 'originClusterHeatmap'),
              plotOutput(outputId = "postMASSclusterHeatmap")),
     tabPanel('venn',plotOutput(outputId = "venn"),
              plotOutput(outputId = "pvalue"), dataTableOutput(outputId = "TFinfo"))
       )
     )
  )
)






server <- function(input, output){
  
  data.path  <- file.path("/Users/Weitinglin/Documents/2016 實驗室資料處理/201510 microarray/raw data 20100114/set 8 CLS CLF Sphere")
  experiment.set <- c ( rep ("set8_withfibroblast" , 3 ) , rep ( "set8_withoutfibroblast" , 3 ), rep("set8_sphere", 3) )
  CLF_CLS.norm.data.path    <- file.path("/Users/Weitinglin/Documents/R_scripts/Lab/microarry20161107andshiny/CLF_CLS.norm.txt")
  CLF_sphere.norm.data.path <- file.path("/Users/Weitinglin/Documents/R_scripts/Lab/microarry20161107andshiny/CLF_sphere.norm.txt")
  CLS_sphere.norm.data.path <- file.path("/Users/Weitinglin/Documents/R_scripts/Lab/microarry20161107andshiny/CLS_sphere.norm.txt")
  probe.relatedTF.data.path <- file.path("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/proberelateTF.txt")
  celfile.set <- do_phenodata(data.path = data.path, experiment.set = experiment.set)
  print("Finish loading the celfile")

  #=====================================================================================
  postMASS <- read.delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarry20161107andshiny/postMASS.txt",
                         sep ="\t",
                         header = TRUE,
                         check.names = FALSE)
  print("load the MASS origin file")
  CLF_CLS.norm <- read.delim(CLF_CLS.norm.data.path,
                             sep ="\t",
                             header = TRUE,
                             check.names = FALSE)
  CLF_CLS.norm <- data.matrix(CLF_CLS.norm)
  CLF_sphere.norm <- read.delim(CLF_sphere.norm.data.path,
                                sep ="\t",
                                header = TRUE,
                                check.names = FALSE)
  CLF_sphere.norm <- data.matrix(CLF_sphere.norm)
  CLS_sphere.norm  <- read.delim(CLS_sphere.norm.data.path,
                                 sep ="\t",
                                 header = TRUE,
                                 check.names = FALSE)
  CLS_sphere.norm <- data.matrix(CLS_sphere.norm)
  print("Finish loading the post quantilenormalization norm")

  
  #====================================================================
  #packed back to affybash
  QN_CLF_CLS_eSet <-new("ExpressionSet",
                                    phenoData = phenoData(celfile.set[,1:6]),
                                    annotation   = annotation(celfile.set[,1:6]),
                                    protocolData = protocolData(celfile.set[,1:6]),
                                    experimentData = experimentData(celfile.set[,1:6]),
                                    exprs = CLF_CLS.norm )
  QN_CLF_Sphere_eSet <-new("ExpressionSet",
                                    phenoData = phenoData(celfile.set[,c(1:3,7:9)]),
                                    annotation   = annotation(celfile.set[,c(1:3,7:9)]),
                                    protocolData = protocolData(celfile.set[,c(1:3,7:9)]),
                                    experimentData = experimentData(celfile.set[,c(1:3,7:9)]),
                                    exprs = CLF_sphere.norm )
  QN_CLS_Sphere_eSet <-new("ExpressionSet",
                                    phenoData = phenoData(celfile.set[,4:9]),
                                    annotation   = annotation(celfile.set[,4:9]),
                                    protocolData = protocolData(celfile.set[,4:9]),
                                    experimentData = experimentData(celfile.set[,4:9]),
                                    exprs = CLS_sphere.norm )
  
  
   print("Finish create new expression set")
  #====================================================================
   table.probe.relate.TF <- read_delim(probe.relatedTF.data.path, delim = " ", col_types = c("ciiii"))
   
   
   
   
   #====================================================================
   
   prefilter.parameters <- reactive({
     if(input$prefilter == "filterbyvariance"){
         
           c(input$variance_value)
       }else if (input$prefilter == "filterbyexpressionlevel"){
           c(input$ratio, input$exprlevel)
       }#pOverA(1/3, log2(100))
       })
   
   eDATA.CLF_CLS    <-  reactive({prefilter_choose(eSet = QN_CLF_CLS_eSet, method = input$prefilter, parameters = prefilter.parameters() )})
   eDATA.CLF_Sphere <-  reactive({prefilter_choose(eSet = QN_CLF_Sphere_eSet, method = input$prefilter, parameters = prefilter.parameters() )})
   eDATA.CLS_Sphere <-  reactive({prefilter_choose(eSet = QN_CLS_Sphere_eSet, method = input$prefilter, parameters = prefilter.parameters() )})
   print("Finish the prefiltering")
   de.parameters <- reactive({ 
                                if(input$demethod=="limma"){
                                        c(input$pvalue)
                                                          }else if(input$demethod=="wilcox"){
                                         c(input$prior, input$sig)
                               }})
   
   top.CLF_CLS      <- reactive({tmp<- testing_choose(eSet = eDATA.CLF_CLS(), method=input$demethod, parameters = de.parameters(), group = 1)
                                  print(head(tmp))
                                  tmp })
   top.CLF_Sphere   <- reactive({tmp<- testing_choose(eSet = eDATA.CLF_Sphere(), method=input$demethod, parameters = de.parameters(), group = 2)
                                   print(head(tmp))
                                   tmp })
   top.CLS_Sphere   <- reactive({tmp<- testing_choose(eSet = eDATA.CLS_Sphere(), method=input$demethod, parameters = de.parameters(), group = 3)
                                  print(head(tmp))
                                   tmp                               })
   
   CLF_CLS.level    <- gl(n = 2, k =3, labels = c("CLF", "CLS"))
   CLF_sphere.level <- gl(n = 2, k =3, labels = c("CLF", "Sphere"))
   CLS_sphere.level <- gl(n = 2, k =3, labels = c("CLS", "Sphere"))
   
   mw.CLF_CLS       <- reactive({
                                MultiWilcoxonTest(eDATA.CLF_CLS(), CLF_CLS.level)
                                })
   mw.CLF_Sphere    <- reactive({
                               MultiWilcoxonTest(eDATA.CLF_Sphere(), CLF_sphere.level)
                               })
   mw.CLS_Sphere    <- reactive({ 
                               MultiWilcoxonTest(eDATA.CLS_Sphere(), CLS_sphere.level)
                                 })
   
   
   print("Finish the testing")
   
   CLF_CLS.probe    <- reactive({
                                print("CLF_CLS.prob")
                                print(head(rownames(top.CLF_CLS())))
                                rownames(top.CLF_CLS())})
   CLF_sphere.probe <- reactive({
                                print("CLF_sphere.prob")
                                print(head(rownames(top.CLF_Sphere())))
                                rownames(top.CLF_Sphere())})
   CLS_sphere.probe <- reactive({
                                print("CLS_sphere.prob")  
                                print(head(rownames(top.CLF_Sphere())))
                                rownames(top.CLS_Sphere())})
   
   area.list <- reactive({
   print("in the area list") 
   CLF_CLS.probe <- CLF_CLS.probe()
   CLF_sphere.probe <- CLF_sphere.probe()
   CLS_sphere.probe <- CLS_sphere.probe()
   area.1 <- length(CLF_CLS.probe)
   area.2 <- length(CLF_sphere.probe)
   area.3 <- length(CLS_sphere.probe)
   area.12 <- length(intersect(CLF_CLS.probe, CLF_sphere.probe))
   area.23 <- length(intersect(CLF_sphere.probe,CLS_sphere.probe))
   area.13 <- length(intersect(CLF_CLS.probe, CLS_sphere.probe))
   area.123 <- length(intersect(intersect(CLF_CLS.probe, CLF_sphere.probe),CLS_sphere.probe))
   list(area.1, area.2, area.3, area.12, area.23, area.13, area.123)
   })
   print("Finish the subsetting")
   
   
   area.23.probe <- reactive({
     intersect(intersect(CLF_sphere.probe(),CLS_sphere.probe()), intersect(CLF_CLS.probe(), CLS_sphere.probe()))
   })
   
   
   
   
   #====================================================================
   
   
output$denDIS <- renderPlot({
  
  total.expression <- log2(postMASS+1) %>% as.data.frame()
  CLF.expression <- total.expression[,1:3]%>% tidyr::gather("sample","value") %>% mutate(case="CLF")
  CLS.expression <- total.expression[,4:6]%>% tidyr::gather("sample","value") %>% mutate(case="CLS") 
  Sphere.expression <- total.expression[,7:9]%>% tidyr::gather("sample","value") %>% mutate(case="Sphere")
  total.expression <- bind_rows(CLF.expression, CLS.expression, Sphere.expression)
  ggplot(data = total.expression, aes(x = value)) + geom_density(aes(color = get(input$origin.effect) ))
})   

output$originClusterHeatmap <- renderPlot({
  make_clusterheatmap(celfile.set, expression=FALSE)
})

output$postMASSclusterHeatmap <- renderPlot({
  make_clusterheatmap(postMASS, expression=TRUE)
})

output$pvalue <- renderPlot({
  print("Make pvalue plot")
  if(input$demethod=="limma"){
    par(mfrow=c(3,1))
    print(str(top.CLF_CLS()))
    hist(top.CLF_CLS()$P.Value)
    hist(top.CLF_CLS()$P.Value)
    hist(top.CLF_CLS()$P.Value)
  }else if(input$demethod=="wilcox"){
    par(mfrow=c(3,2))
    hist(mw.CLF_CLS())
    plot(mw.CLF_CLS())
    hist(mw.CLF_Sphere())
    plot(mw.CLF_Sphere())
    hist(mw.CLS_Sphere())
    plot(mw.CLS_Sphere())
  }
  
  
})   
   
           
output$venn <- renderPlot({
  print("Begin to make vennDiagram")
  area.list <- area.list()
  draw.triple.venn(area1 = area.list[[1]], area2 = area.list[[2]], area3 = area.list[[3]], n12 = area.list[[4]], n23 = area.list[[5]], n13 = area.list[[6]],
                   n123 = area.list[[7]], category = c("CLF_CLS", "CLF_Sphere", "CLS_Sphere"), lty = "blank",
                   fill = c("skyblue", "pink1", "mediumorchid"))
})

output$TFinfo <- renderDataTable({
  table.probe.relate.TF %>% filter(Probe %in% area.23.probe())
})
  
}


shinyApp( ui = ui, server = server)