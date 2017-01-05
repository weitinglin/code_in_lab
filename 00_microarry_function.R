  
# Loading package ---------------------------------------------------------
#library(affy) #
#library(GEOquery) #
#library(hgu133plus2cdf) #
#bioconductor
#source("https://bioconductor.org/biocLite.R")

#Phenodata and load the celfile
#data.path  <- file.path("/Users/Weitinglin/Documents/2016 實驗室資料處理/201510 microarray/raw data 20100114/set 8 CLS CLF Sphere")
#experiment.set <- c ( rep ("set8_withfibroblast" , 3 ) , rep ( "set8_withoutfibroblast" , 3 ), rep("set8_sphere", 3) )

do_phenodata <- function(data.path, experiment.set){
  
 set<- list.files (data.path, pattern=".CEL")
 phenodata.set               <- matrix ( rep ( set, 2) , ncol = 2 )
 phenodata.set      <- as.data.frame ( phenodata.set )
 colnames ( phenodata.set )   <- c ( "Name" , "FileName" )
 phenodata.set$experiment.set <- experiment.set
 write.table(phenodata.set, paste(data.path,"/phenodata_set.txt", sep = ""),quote=F,sep="\t",row.names=F)
 print("save the file! in to a phenodata_set.txt")
 celfile.set <- ReadAffy ( celfile.path = data.path , phenoData = paste(data.path,"/phenodata_set.txt", sep = "") )
 return(celfile.set)
 }


#Quantile Normalization
quantile_normalization <- function(ecelfile.set, method="mean"){
  a <- ecelfile.set
  asort <- apply(a,2,sort)
  loc<-apply(a,2,order)  #原本probe的位置apply
  #step 2 :Takes the means of across rows
  if (method =="mean"){
  amean1<-matrix(rep(apply(asort,1,mean),dim(a)[2]),nrow=length(apply(asort,1,mean)))
  }else if(method == "median"){
    amean1<-matrix(rep(apply(asort,1,median),dim(a)[2]),nrow=length(apply(asort,1,median)))
  }
  #step 3 :Rearranging each column to make the same ordering as the original data
  norm<-matrix(0,dim(a)[1],dim(a)[2])
  
  for(i in 1:dim(amean1)[2]){
    #norm<-matrix(0,dim(amean)[1],1)
    for (j in 1:dim(amean1)[1]){
      #cat(loc[j,i],"\t",i,"\t",j,"\t",amean[j,i],"\n")
      norm[loc[j,i],i]<-amean1[j,i]
    }
  }
  norm <- log2(norm + 1)
  row.names(norm) <- row.names ( ecelfile.set) 
  colnames(norm)  <- colnames ( ecelfile.set)
  return(norm)
}


#pre-filtering
prefilter_choose <- function(eSet, method, parameters){
  if ( method == "filterbyvariance"){
    print("I'm in filter")
    print(parameters)
    return(nsFilter(eSet, remove.dupEntrez = FALSE, var.cutoff = parameters[1], var.func=IQR, var.filter = TRUE)$eset)
  } else if ( method == "filterbyexpressionlevel"){
    print(parameters)
    f1 <- pOverA(parameters[1], parameters[2]) #pOverA(1/3, log2(100))
    f2 <- function(x)(diff(range(x, na.rm = T)) > log2(1.5))
    ff <- filterfun(f1, f2)
    index <- genefilter(exprs(eSet), ff)
    return(eSet[index, ])
  } else if (method == "nofilter"){
    return(eSet)
  }
}

#testing
testing_choose <- function(eSet, method, parameters, group){
    CLF_CLS.level    <- gl(n = 2, k =3, labels = c("CLF", "CLS"))
    CLF_sphere.level <- gl(n = 2, k =3, labels = c("CLF", "Sphere"))
    CLS_sphere.level <- gl(n = 2, k =3, labels = c("CLS", "Sphere"))
    paired.level.list <- list(CLF_CLS.level, CLF_sphere.level, CLS_sphere.level)
    contrast.vector <- c("CLF - CLS", "CLF - Sphere", "CLS - Sphere")
    print(paste("Doing the group", group ))
   if (method == "limma") {
     print("Choose the method limma + bayes")
     design.list <- model.matrix(~0+paired.level.list[[group]])
     colnames(design.list) <- levels( paired.level.list[[group]] )
     print("HERE2")
     fit.list <- lmFit(eSet, design.list)
     print(paste("Finish the limma fit in group",group))
    contrast.list <- makeContrasts(contrasts = contrast.vector[group], levels = design.list)
    fit.list <- contrasts.fit(fit.list, contrast.list)
    fit.list <- eBayes(fit.list)
    print(paste("Finish the Bayes in group",group))
    top.list <- topTable(fit.list, coef=contrast.vector[group], number=nrow(fit.list), adjust.method = "BH", p.value = parameters[1]) #parameter=c(0.05)
    print(paste("Finish the testing in group",group))
    return(top.list)
  } else if ( method == "wilcox") {
    print("Choose the method wilcox rank-sum test + bayes")
    print(parameters)
    mw.list <- MultiWilcoxonTest(eSet, paired.level.list[[group]])
    mw <- mw.list
    print("HERE IN WILCOX")
    tidy.data <- data.frame(probe = names(mw@rank.sum.statistics), rank.sum.statistics = unname(mw@rank.sum.statistics))
    rownames(tidy.data) <- names(mw@rank.sum.statistics[])
    tidy.data <- tidy.data[unname(selectSignificant(mw, prior = parameters[1], significance = parameters[2])),]  # parameters = c(0.1, 0.05)
    print(head(tidy.data))
    return(tidy.data)
  }
}

#wilcox.test
testing_wilcox_old <- function(eSet){
  print("use homemade wilcox")
  expression <- exprs(eSet)
  wilcox.result <- data.frame()
  for (i in 1:nrow(expression )){
    wilcox.result <- wilcox.test(expression[i,1:3],expression[i,4:6], correct = FALSE, paired = FALSE) %>%
      tidy %>% mutate(gene = rownames(expression)[i]) %>%
      bind_rows(wilcox.result,.)
    
  }
  rownames(wilcox.result) <- wilcox.result$gene
  
  return(wilcox.result)
}

#whether a probe is a TF or not
# TF.BP.GO.path <- file.path("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/GO_BP_TF.txt")
# TF.CC.GO.path <- file.path("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/GO_CC_TF.txt")
# TF.MF.GO.path <- file.path("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/GO_MF_TF.txt")
# 
# TF.BP.GO <- read_delim(TF.BP.GO.path, col_names = c("GO_identifier", "description"),delim = "\t")
# TF.CC.GO <- read_delim(TF.CC.GO.path, col_names =c("GO_identifier", "description"), delim = "\t")
# TF.MF.GO <- read_delim(TF.MF.GO.path, col_names = c("GO_identifier", "description"), delim = "\t")
# 
# TF.BP.GO.id <- TF.BP.GO$GO_identifier 
# TF.CC.GO.id <- TF.CC.GO$GO_identifier 
# TF.MF.GO.id <- TF.MF.GO$GO_identifier


getprobeRelateToGoTF <- function(probelist){
probe.relatedTF <- data.frame()
for (i in 1:length(probelist)){
  index <- probelist[[i]]
  GO.list <- hgu133plus2GO[index] %>% toTable
  GO.BP.list <- GO.list %>% filter(Ontology == "BP")
  GO.CC.list <- GO.list %>% filter(Ontology == "BP")
  GO.MF.list <- GO.list %>% filter(Ontology == "BP")
  BP <- sum(GO.BP.list$go_id %in% TF.BP.GO.id)  
  CC <- sum(GO.CC.list$go_id %in% TF.CC.GO.id)
  GO <- sum(GO.MF.list$go_id %in% TF.MF.GO.id)
  result <- data.frame(Probe = index, BP = BP, CC = CC, GO = GO )
  probe.relatedTF <- bind_rows(probe.relatedTF, result)
}
return(probe.relatedTF)
}


#clustering heatmap function to celfile
make_clusterheatmap <- function(eSet, expression=FALSE){
  if (expression){
  dd <- dist2(log2(eSet+1))}
  else{
  dd <- dist2(log2(exprs(eSet)+1))}
  diag(dd) <- 0
  dd.row <- as.dendrogram(hclust(as.dist(dd)))
  row.ord <- order.dendrogram(dd.row)
  legend <- list(top = list(fun= dendrogramGrob, args = list(x = dd.row, side="top")))
  lp <- levelplot(dd[row.ord, row.ord], scales = list(x = list(rot=90)), xlab="", ylab="",legend=legend)
  plot(lp)
}

#testing_wilcox with foreach
testing_wilcox <- function(eSet, adj.method="BH"){
  registerDoParallel(cores=4)
  expression <- exprs(eSet)
  
  wilcox.result  <- foreach(i = 1:nrow(expression)) %dopar% {
    wilcox.test(expression[i,1:3],expression[i,4:6], correct = FALSE, paired = FALSE) %>%
      tidy %>% mutate(gene = rownames(expression)[i])
  }
  wilcox.result <- invoke(rbind, map(wilcox.result,data.frame))
  #do.call(rbind, lapply(wilcox.result,data.frame))
  wilcox.result <- as.data.frame(wilcox.result)
  rownames(wilcox.result) <- wilcox.result$gene
  wilcox.result <- wilcox.result %>%
    dplyr::select(gene,statistic,p.value,-method, -alternative)%>%
    mutate(adjusted.p = p.adjust(p.value, method=adj.method),
                                            adjusted.method = adj.method)
  return(wilcox.result)
}




# api ---------------------------------------------------------------------
#depent on purrr, dplyr, httr, rjson
searchHarmonizome <- function(genelist){
    dataframe <- data.frame(symbol = c(),
                            name = c(),
                            synonyms = c(),
                            description = c())
    url <- "http://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene/"
    checkExist <- function(request.result){
        if(is.null(request.result)){
            return("NULL")
        }else if(nchar(request.result) == 0){
            return("NULL")
        }else(return(request.result))
    }
    
    
    for ( item in genelist){
        raw <- paste0(url,item) %>% GET %>% httr::content(.,"text", encoding = "ISO-8859-1") %>% fromJSON
        
        gene.symbol       <- checkExist(raw$symbol)
        gene.name         <- checkExist(raw$name)
        gene.synonyms     <- checkExist(purrr::reduce(raw$synonyms, paste, sep=","))
        gene.description  <- checkExist(purrr::reduce(raw$description, paste, sep=","))
        
        
        
        request.dataframe <- data.frame(symbol = gene.symbol,
                                        name = gene.name,
                                        synonyms = gene.synonyms,
                                        description = gene.description)
        dataframe <- bind_rows(dataframe, request.dataframe)
    }
    return(dataframe)
    
}




