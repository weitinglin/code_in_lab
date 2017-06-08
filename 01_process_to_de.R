#microarray 2016/11/02
#20160316
# Mon Dec 12 18:07:35 2016 ------------------------------


#library in the source R code

source("/Users/Weitinglin/Documents/Repository/code_in_lab/00_microarry_function.R")






#=====================input the data================================
data.path  <- file.path("/Volumes/2016weiting/lab_data/validationdataset/Raw_data/GSE19804/GSE19804_RAW")

experiment.set <- c ( rep ("set8_withfibroblast" , 3 ) , rep ( "set8_withoutfibroblast" , 3 ), rep("set8_sphere", 3) )

set<- list.files (data.path, pattern=".CEL.gz")
phenodata.set   <- matrix ( rep ( set, 2) , ncol = 2 )
phenodata.set   <- as.data.frame ( phenodata.set )
colnames ( phenodata.set )   <- c ( "Name" , "FileName" )
phenodata.set$experiment.set <- experiment.set
write.table(phenodata.set, paste(data.path,"/phenodata_set.txt", sep = ""),quote=F,sep="\t",row.names=F)
print("save the file! in to a phenodata_set.txt")
celfile.set <- ReadAffy ( celfile.path = data.path , phenoData = paste(data.path,"/phenodata_set.txt", sep = "") )

#=====================preprocess the data ======
#use the mas5,in order not to log twice
Exprs.data <- mas5 ( celfile.set , normalize = FALSE, analysis = "absolute", sc = 500)

#location of the file
#/Users/Weitinglin/Documents/R_scripts/Lab/microarray/data/intermediate
save(Exprs.data, file="Exprs_data.RData")
load("Exprs_data.RData")
##expression profile
ecelfile.set <- exprs ( Exprs.data )



# expression boxplot ------------------------------------------------------
boxplot.data <- as.data.frame(log2(ecelfile.set+1))
boxplot.data$probe <- rownames(boxplot.data)
boxplot.data <- boxplot.data %>% tidyr::gather("Sample","value",-probe)
ggplot(data = boxplot.data, aes(x = Sample, y = value)) +
  geom_boxplot(aes(fill = Sample)) 


# cluster -----------------------------------------------------------------

##clustering 
library(latticeExtra)
dd <- dist2(log2(ecelfile.set+1))
dd <- dist2(norm)
dd <- as.matrix(dist(t(log2(Exprs.data.median+1)), method = "euclidean"))
#euclidean, maximum, manhattan, canberra, binary, mikowski
dd <- dist2(log2(exprs(Exprs.data )+1))
diag(dd) <- 0
dd.row <- as.dendrogram(hclust(as.dist(dd)))
row.ord <- order.dendrogram(dd.row)
legend <- list(top = list(fun= dendrogramGrob, args = list(x = dd.row, side="top")))
lp <- levelplot(dd[row.ord, row.ord], scales = list(x = list(rot=90)), xlab="", ylab="",legend=legend)
plot(lp)


# PCA ---------------------------------------------------------------------
t <- t(ecelfile.set)
t <- t(norm)
pca1 <- prcomp(t, scale. = TRUE)
library(scatterplot3d)
s3d <- scatterplot3d(x = pca1$x[,1], y = pca1$x[,2], z = pca1$x[,3], xlab="PC1", ylab="PC2", zlab="PC3")
s3d.coords <- s3d$xyz.convert(x = pca1$x[,1], y = pca1$x[,2], z = pca1$x[,3])
text(s3d.coords$x, s3d.coords$y, labels=row.names(pca1$x), cex = .5, pos=4)

library(plotly)
p <- plot_ly(test, x = ~PC1, y = ~PC2, z = ~PC3, color = ~sample, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))
#========================quantile normalization===========================
#Method Quantile Normalization+log2

norm <- quantile_normalization(ecelfile.set = Exprs.data, method = "median")


# boxplot after QN --------------------------------------------------------
boxplot.data <- as.data.frame(norm)
boxplot.data$probe <- rownames(boxplot.data)
boxplot.data <- boxplot.data %>% tidyr::gather("Sample","value",-probe)
ggplot(data = boxplot.data, aes(x = Sample, y = value)) +
  geom_boxplot(aes(fill = Sample)) 







norm<- data.matrix(norm)

#packed back to affybash
Quantile_Normalization_eSet <-new("ExpressionSet", 
                                  phenoData = phenoData(celfile.set), 
                                  annotation   = annotation(celfile.set),
                                  protocolData = protocolData(celfile.set),
                                  experimentData = experimentData(celfile.set),
                                  exprs = norm )
#GO exploration
probererlateTF <- read_delim("/home/weitinglin66/LungCSC_microarray/rawdata/proberelateTF.txt", delim = " ")

#=========================experiment set up=======================

#filter the genes that at least 2 out of 6 samples have an intensity of 100,
#at least max/min intensity is 1.5

#plot row standard deviations versus row means
library(vsn)
meanSdPlot(eDATA)


# EX01: MAS5,QN,LOG2, without filter, wilcoxon rank-sum  --------------------------------

#without filter
top.CLF_CLS      <- testing_wilcox(Quantile_Normalization_eSet[,1:6], adj.method = "BH")
top.CLF_Sphere   <- testing_wilcox(Quantile_Normalization_eSet[,c(1:3,7:9)], adj.method = "BH")
top.CLS_Sphere   <- testing_wilcox(Quantile_Normalization_eSet[,4:9], adj.method = "BH")

write.csv(top.CLF_CLS, file = "CLF_CLS_wilcox_adj_nofilter.csv")
write.csv(top.CLF_Sphere, file = "CLF_Sphere_wilcox_adj_nofilter.csv")
write.csv(top.CLS_Sphere, file = "CLS_Sphere_wilcox_adj_nofilter.csv")


CLF_CLS.probe    <- (top.CLF_CLS %>% filter(p.value < 0.05))$gene
CLF_sphere.probe <- (top.CLF_Sphere %>% filter(p.value < 0.05))$gene
CLS_sphere.probe <-  (top.CLS_Sphere %>% filter(p.value < 0.05))$gene

print("in the area list") 
area.1 <- length(CLF_CLS.probe)
area.2 <- length(CLF_sphere.probe)
area.3 <- length(CLS_sphere.probe)
area.12 <- length(intersect(CLF_CLS.probe, CLF_sphere.probe))
area.23 <- length(intersect(CLF_sphere.probe,CLS_sphere.probe))
area.13 <- length(intersect(CLF_CLS.probe, CLS_sphere.probe))
area.123 <- length(intersect(intersect(CLF_CLS.probe, CLF_sphere.probe),CLS_sphere.probe))
area.list <- list(area.1, area.2, area.3, area.12, area.23, area.13, area.123)
#intersect, union, setdiff
draw.triple.venn(area1 = area.list[[1]], area2 = area.list[[2]], area3 = area.list[[3]], n12 = area.list[[4]], n23 = area.list[[5]], n13 = area.list[[6]],
                 n123 = area.list[[7]], category = c("CLF_CLS", "CLF_Sphere", "CLS_Sphere"), lty = "blank",
                 fill = c("skyblue", "pink1", "mediumorchid"))

target.list <- setdiff(intersect(CLF_CLS.probe, CLS_sphere.probe), CLF_sphere.probe )

target.probe <- probererlateTF %>% filter(Probe %in% target.list)

getSYMBOL(target.probe$Probe, "hgu133plus2")
#GO exploration
probererlateTF <- read_delim("/home/weitinglin66/LungCSC_microarray/rawdata/proberelateTF.txt", delim = " ")



#withoutfilter but seperate normalization
CLF_CLS.norm.data.path    <- file.path("/home/weitinglin66/LungCSC_microarray/result/CLF_CLS.norm.txt")
CLF_sphere.norm.data.path <- file.path("/home/weitinglin66/LungCSC_microarray/result/CLF_sphere.norm.txt")
CLS_sphere.norm.data.path <- file.path("/home/weitinglin66/LungCSC_microarray/result/CLS_sphere.norm.txt")
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


sepQN.CLF_CLS      <- testing_wilcox(QN_CLF_CLS_eSet, adj.method = "BH")
sepQN.CLF_Sphere   <- testing_wilcox(QN_CLF_Sphere_eSet, adj.method = "BH")
sepQN.CLS_Sphere   <- testing_wilcox(QN_CLS_Sphere_eSet, adj.method = "BH")

write.csv(sepQN.CLF_CLS, file = "sepQNCLF_CLS_wilcox_adj_nofilter.csv")
write.csv(sepQN.CLF_Sphere, file = "sepQNCLF_Sphere_wilcox_adj_nofilter.csv")
write.csv(sepQN.CLS_Sphere, file = "sepQNCLS_Sphere_wilcox_adj_nofilter.csv")


CLF_CLS.probe    <- (sepQN.CLF_CLS %>% filter(p.value < 0.05))$gene
CLF_sphere.probe <- (sepQN.CLF_Sphere %>% filter(p.value < 0.05))$gene
CLS_sphere.probe <-  (sepQN.CLS_Sphere %>% filter(p.value < 0.05))$gene

sepQN.CLF_CLS    <- sepQN.CLF_CLS %>% mutate(case="CLF_CLS")
sepQN.CLF_Sphere <- sepQN.CLF_Sphere %>% mutate(case =" CLF_Sphere")
sepQN.CLS_Sphere <- sepQN.CLS_Sphere %>% mutate(case = "CLS_Sphere")
total <- bind_rows(sepQN.CLF_CLS, sepQN.CLF_Sphere, sepQN.CLS_Sphere)
(total %>% filter(case == "CLF_CLS") %>% filter( p.value < 0.1))$p.value %>% hist(.,main = "CLS_CLF",breaks=seq(0,0.1,by=0.01)) + abline(v=0.05, col="red")
(total %>% filter(case == "CLS_Sphere") %>% filter( p.value < 0.1))$p.value  %>% hist(., main ="CLS_Sphere",breaks=seq(0,0.1,by=0.01)) + abline(v=0.05, col="red")
(total %>% filter(case != "CLS_Sphere" & case != "CLF_CLS") %>% filter( p.value < 0.1))$p.value %>% hist(., main = "CLF_Sphere",breaks=seq(0,0.1,by=0.01)) + abline(v=0.05, col="red")

(total %>% filter(case == "CLF_CLS") %>% filter( adjusted.p < 0.5))$adjusted.p %>% hist(.,main = "CLS_CLF",breaks=seq(0,0.5,by=0.01)) + abline(v=0.05, col="red")
(total %>% filter(case == "CLS_Sphere") %>% filter( adjusted.p < 0.5))$adjusted.p  %>% hist(., main ="CLS_Sphere",breaks=seq(0,0.5,by=0.01)) + abline(v=0.05, col="red")
(total %>% filter(case != "CLS_Sphere" & case != "CLF_CLS") %>% filter( adjusted.p < 0.5))$adjusted.p %>% hist(., main = "CLF_Sphere",breaks=seq(0,0.5,by=0.01)) + abline(v=0.05, col="red")

ggplot(data=total, aes(x=p.value)) + geom_density(aes(color=case))


print("in the area list") 
area.1 <- length(CLF_CLS.probe)
area.2 <- length(CLF_sphere.probe)
area.3 <- length(CLS_sphere.probe)
area.12 <- length(intersect(CLF_CLS.probe, CLF_sphere.probe))
area.23 <- length(intersect(CLF_sphere.probe,CLS_sphere.probe))
area.13 <- length(intersect(CLF_CLS.probe, CLS_sphere.probe))
area.123 <- length(intersect(intersect(CLF_CLS.probe, CLF_sphere.probe),CLS_sphere.probe))
area.list <- list(area.1, area.2, area.3, area.12, area.23, area.13, area.123)
#intersect, union, setdiff
draw.triple.venn(area1 = area.list[[1]], area2 = area.list[[2]], area3 = area.list[[3]], n12 = area.list[[4]], n23 = area.list[[5]], n13 = area.list[[6]],
                 n123 = area.list[[7]], category = c("CLF_CLS", "CLF_Sphere", "CLS_Sphere"), lty = "blank",
                 fill = c("skyblue", "pink1", "mediumorchid"))

target.list <- setdiff(intersect(CLF_CLS.probe, CLS_sphere.probe), CLF_sphere.probe )

target.probe <- probererlateTF %>% filter(Probe %in% target.list)


# EX02: MAS5,QN,LOG2, with filtering by 50%, wilcoxon-rank-sum test -------
#filter for wilcox test on variance 50%
postfilteQN_CLF_CLS_eSet <- prefilter_choose(QN_CLF_CLS_eSet, method = "filterbyvariance", parameters = c(0.5))
postfilteQN_CLF_Sphere_eSet <- prefilter_choose(QN_CLF_Sphere_eSet, method = "filterbyvariance", parameters = c(0.5))
postfilteQN_CLS_Sphere_eSet <- prefilter_choose(QN_CLS_Sphere_eSet, method = "filterbyvariance", parameters = c(0.5))


postfilteCLF_CLS      <- testing_wilcox(postfilteQN_CLF_CLS_eSet, adj.method = "BH")
postfilteCLF_Sphere   <- testing_wilcox(postfilteQN_CLF_Sphere_eSet, adj.method = "BH")
postfilteCLS_Sphere   <- testing_wilcox(postfilteQN_CLS_Sphere_eSet, adj.method = "BH")

CLF_CLS.probe    <- (postfilteCLF_CLS %>% filter(p.value < 0.05))$gene
CLF_sphere.probe <- (postfilteCLF_Sphere %>% filter(p.value < 0.05))$gene
CLS_sphere.probe <-  (postfilteCLS_Sphere %>% filter(p.value < 0.05))$gene

print("in the area list") 
area.1 <- length(CLF_CLS.probe)
area.2 <- length(CLF_sphere.probe)
area.3 <- length(CLS_sphere.probe)
area.12 <- length(intersect(CLF_CLS.probe, CLF_sphere.probe))
area.23 <- length(intersect(CLF_sphere.probe,CLS_sphere.probe))
area.13 <- length(intersect(CLF_CLS.probe, CLS_sphere.probe))
area.123 <- length(intersect(intersect(CLF_CLS.probe, CLF_sphere.probe),CLS_sphere.probe))
area.list <- list(area.1, area.2, area.3, area.12, area.23, area.13, area.123)
#intersect, union, setdiff
draw.triple.venn(area1 = area.list[[1]], area2 = area.list[[2]], area3 = area.list[[3]], n12 = area.list[[4]], n23 = area.list[[5]], n13 = area.list[[6]],
                 n123 = area.list[[7]], category = c("CLF_CLS", "CLF_Sphere", "CLS_Sphere"), lty = "blank",
                 fill = c("skyblue", "pink1", "mediumorchid"))

target.list <- setdiff(intersect(CLF_CLS.probe, CLS_sphere.probe), CLF_sphere.probe )

target.probe <- probererlateTF %>% filter(Probe %in% target.list)





# EX03: MASS, QN, no LOG2, Wilcoxon-rank sum test -------------------------
# 2016.11.15
ecelfile.set <- exprs ( Exprs.data )



# do the quantile normalization 
a <- ecelfile.set
asort <- apply(a,2,sort)
loc<-apply(a,2,order)  #原本probe的位置apply
#step 2 :Takes the median of across rows
amean1<-matrix(rep(apply(asort,1,median),dim(a)[2]),nrow=length(apply(asort,1,median)))
#step 3 :Rearranging each column to make the same ordering as the original data
norm<-matrix(0,dim(a)[1],dim(a)[2])

for(i in 1:dim(amean1)[2]){
  #norm<-matrix(0,dim(amean)[1],1)
  for (j in 1:dim(amean1)[1]){
    #cat(loc[j,i],"\t",i,"\t",j,"\t",amean[j,i],"\n")
    norm[loc[j,i],i]<-amean1[j,i]
  }
}

rownames(norm) <- rownames(ecelfile.set)
colnames(norm) <- colnames(ecelfile.set)
#check the result of the QN
boxplot(log2(norm+1))



adj.method <- "BH"

registerDoParallel(cores=8)
#function
expression <- norm[1:1000,1:6]

wilcoxon_test <- function(expression, adj.method = "BD"){
  gene.list <- rownames(expression)  
  wilcox.result  <- foreach(i = 1:nrow(expression)) %dopar% {
    wilcox.test(expression[i,1:3],expression[i,4:6], correct = FALSE, paired = FALSE) %>%
      tidy %>% mutate(gene = gene.list[i])
  }
  wilcox.result <- invoke(rbind, map(wilcox.result,data.frame))
  
  wilcox.result <- as.data.frame(wilcox.result)
  rownames(wilcox.result) <- wilcox.result$gene
  wilcox.result <- wilcox.result %>%
    dplyr::select(gene,statistic,p.value,-method, -alternative)%>%
    mutate(adjusted.p = p.adjust(p.value, method=adj.method),
           adjusted.method = adj.method)
  return(wilcox.result)
}
#for CLS-CLF
nologQN.CLF_CLS      <- wilcoxon_test(expression = norm[,1:6], adj.method = "BH")
nologQN.CLF_Sphere   <- wilcoxon_test(expression = norm[,c(1:3,7:9)], adj.method = "BH")
nologQN.CLS_Sphere   <- wilcoxon_test(expression = norm[,4:9], adj.method = "BH")

write.table(nologQN.CLF_CLS,"nologQNwilcoxCLF_CLS.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(nologQN.CLF_Sphere,"nologQNwilcoxCLF_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(nologQN.CLS_Sphere,"nologQNwilcoxCLS_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)

#check the p value
CLF_CLS.probe    <- (nologQN.CLF_CLS %>% filter(p.value < 0.05))$gene
CLF_sphere.probe <- (nologQN.CLF_Sphere %>% filter(p.value < 0.05))$gene
CLS_sphere.probe <-  (nologQN.CLS_Sphere %>% filter(p.value < 0.05))$gene

print("in the area list") 
#draw the venndiagram
draw_venngram(CLF_CLS.probe, CLF_sphere.probe, CLS_sphere.probe)

target.list <- setdiff(intersect(CLF_CLS.probe, CLS_sphere.probe), CLF_sphere.probe )

target.probe <- probererlateTF %>% filter(Probe %in% target.list)



# EX04: MAS5, QN, LOG2,without filter, two sample t.test ------------------
ecelfile.set <- exprs ( Exprs.data )

#QN + log2

norm <- quantile_normalization(ecelfile.set = Exprs.data)

#two.sample t.test

adj.method <- "BH"

registerDoParallel(cores=8)
#function


t_test <- function(expression, adj.method = "BH"){
  gene.list <- rownames(expression)  
  t.result  <- foreach(i = 1:nrow(expression)) %dopar% {
    t.test(expression[i,1:3],expression[i,4:6], paired = FALSE,  conf.level = 0.95) %>%
      tidy %>% mutate(gene = gene.list[i])
  }
  t.result <- invoke(rbind, map(t.result,data.frame))
  
  t.result <- as.data.frame(t.result)
  rownames(t.result) <- t.result$gene
  t.result <- t.result %>%
    dplyr::select(-method, -alternative)%>%
    mutate(adjusted.p = p.adjust(p.value, method=adj.method),
           adjusted.method = adj.method)
  return(t.result)
}
#“holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”

#for CLS-CLF
t.CLF_CLS      <- t_test(expression = norm[,1:6], adj.method = "BH")
t.CLF_Sphere   <- t_test(expression = norm[,c(1:3,7:9)], adj.method = "BH")
t.CLS_Sphere   <- t_test(expression = norm[,4:9], adj.method = "BH")

write.table(t.CLF_CLS,"QNttestCLF_CLS.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(t.CLF_Sphere,"QNttestCLF_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(t.CLS_Sphere,"QNttestCLS_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)

hist(t.CLF_CLS$adjusted.p) + abline(v = 0.05, col = "red")
hist(t.CLF_Sphere$adjusted.p) + abline(v = 0.05, col = "red")
hist(t.CLS_Sphere$adjusted.p) + abline(v = 0.05, col = "red")

#pvalue
CLF_CLS.probe    <- (t.CLF_CLS %>% filter(p.value < 0.05))$gene
CLF_sphere.probe <- (t.CLF_Sphere %>% filter(p.value < 0.05))$gene
CLS_sphere.probe <-  (t.CLS_Sphere %>% filter(p.value < 0.05))$gene
#adjusted.p 0.05
CLF_CLS.probe    <- (t.CLF_CLS %>% filter(adjusted.p < 0.05))$gene
CLF_sphere.probe <- (t.CLF_Sphere %>% filter(adjusted.p < 0.05))$gene
CLS_sphere.probe <-  (t.CLS_Sphere %>% filter(adjusted.p < 0.05))$gene


draw_venngram(CLF_CLS.probe, CLF_sphere.probe, CLS_sphere.probe)

length(CLF_CLS.probe)
length(CLF_sphere.probe)
length(CLS_sphere.probe)
#adjusted.p 0.01
CLF_CLS.probe    <- (t.CLF_CLS %>% filter(adjusted.p < 0.01))$gene
CLF_sphere.probe <- (t.CLF_Sphere %>% filter(adjusted.p < 0.01))$gene
CLS_sphere.probe <-  (t.CLS_Sphere %>% filter(adjusted.p < 0.01))$gene


draw_venngram(CLF_CLS.probe, CLF_sphere.probe, CLS_sphere.probe)

length(CLF_CLS.probe)
length(CLF_sphere.probe)
length(CLS_sphere.probe)

#adjusted.p 0.001
CLF_CLS.probe    <- (t.CLF_CLS %>% filter(adjusted.p < 0.001))$gene
CLF_sphere.probe <- (t.CLF_Sphere %>% filter(adjusted.p < 0.001))$gene
CLS_sphere.probe <-  (t.CLS_Sphere %>% filter(adjusted.p < 0.001))$gene


draw_venngram(CLF_CLS.probe, CLF_sphere.probe, CLS_sphere.probe)

length(CLF_CLS.probe)
length(CLF_sphere.probe)
length(CLS_sphere.probe)

#figure out the TF related probeset

target.list <- setdiff(intersect(CLF_CLS.probe, CLS_sphere.probe), CLF_sphere.probe )
length(target.list)
target.probe <- probererlateTF %>% filter(Probe %in% target.list)
dim(target.probe)
target.probe$Total %>% table()


# EX05: MAS5, QN, LOG2, with filter by 0.5, two sample t.test -------------

ecelfile.set <- exprs ( Exprs.data )

#QN + log2

norm <- quantile_normalization(ecelfile.set = Exprs.data)

#filter by variance 
Quantile_Normalization_eSet <-new("ExpressionSet", 
                                  phenoData = phenoData(celfile.set), 
                                  annotation   = annotation(celfile.set),
                                  protocolData = protocolData(celfile.set),
                                  experimentData = experimentData(celfile.set),
                                  exprs = norm )
QN_CLF_CLS_eSet     <- Quantile_Normalization_eSet[,1:6]
QN_CLF_Sphere_eSet  <- Quantile_Normalization_eSet[,c(1:3,7:9)]
QN_CLS_Sphere_eSet  <- Quantile_Normalization_eSet[,c(4:9)]

filterQN_CLF_CLS   <- prefilter_choose(QN_CLF_CLS_eSet , method = "filterbyvariance", parameters = c(0.5))
filterQN_CLF_Sphere<- prefilter_choose(QN_CLF_Sphere_eSet , method = "filterbyvariance", parameters = c(0.5))
filterQN_CLS_Sphere<- prefilter_choose(QN_CLS_Sphere_eSet , method = "filterbyvariance", parameters = c(0.5))

#for CLS-CLF
f.t.CLF_CLS      <- t_test(expression = exprs(filterQN_CLF_CLS), adj.method = "BH")
f.t.CLF_Sphere   <- t_test(expression = exprs(filterQN_CLF_Sphere), adj.method = "BH")
f.t.CLS_Sphere   <- t_test(expression = exprs(filterQN_CLS_Sphere), adj.method = "BH")

write.table(f.t.CLF_CLS,"filter50QNttestCLF_CLS.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(f.t.CLF_Sphere,"filter50QNttestCLF_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(f.t.CLS_Sphere,"filter50QNttestCLS_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)

hist(f.t.CLF_CLS$adjusted.p) + abline(v = 0.05, col = "red")
hist(f.t.CLF_Sphere$adjusted.p) + abline(v = 0.05, col = "red")
hist(f.t.CLS_Sphere$adjusted.p) + abline(v = 0.05, col = "red")

#pvalue
CLF_CLS.probe    <- (f.t.CLF_CLS %>% filter(p.value < 0.05))$gene
CLF_sphere.probe <- (f.t.CLF_Sphere %>% filter(p.value < 0.05))$gene
CLS_sphere.probe <-  (f.t.CLS_Sphere %>% filter(p.value < 0.05))$gene
#adjusted.p 0.05
CLF_CLS.probe    <- (f.t.CLF_CLS %>% filter(adjusted.p < 0.05))$gene
CLF_sphere.probe <- (f.t.CLF_Sphere %>% filter(adjusted.p < 0.05))$gene
CLS_sphere.probe <-  (f.t.CLS_Sphere %>% filter(adjusted.p < 0.05))$gene


draw_venngram(CLF_CLS.probe, CLF_sphere.probe, CLS_sphere.probe)

length(CLF_CLS.probe)
length(CLF_sphere.probe)
length(CLS_sphere.probe)
#adjusted.p 0.01
CLF_CLS.probe    <- (f.t.CLF_CLS %>% filter(adjusted.p < 0.01))$gene
CLF_sphere.probe <- (f.t.CLF_Sphere %>% filter(adjusted.p < 0.01))$gene
CLS_sphere.probe <-  (f.t.CLS_Sphere %>% filter(adjusted.p < 0.01))$gene


draw_venngram(CLF_CLS.probe, CLF_sphere.probe, CLS_sphere.probe)

length(CLF_CLS.probe)
length(CLF_sphere.probe)
length(CLS_sphere.probe)

#adjusted.p 0.001
CLF_CLS.probe    <- (f.t.CLF_CLS %>% filter(adjusted.p < 0.001))$gene
CLF_sphere.probe <- (f.t.CLF_Sphere %>% filter(adjusted.p < 0.001))$gene
CLS_sphere.probe <-  (f.t.CLS_Sphere %>% filter(adjusted.p < 0.001))$gene


draw_venngram(CLF_CLS.probe, CLF_sphere.probe, CLS_sphere.probe)

length(CLF_CLS.probe)
length(CLF_sphere.probe)
length(CLS_sphere.probe)

#figure out the TF related probeset

target.list <- setdiff(intersect(CLF_CLS.probe, CLS_sphere.probe), CLF_sphere.probe )
length(target.list)
target.probe <- probererlateTF %>% filter(Probe %in% target.list)
dim(target.probe)
target.probe$Total %>% table()
target.probe$Probe %>% getSYMBOL(., "hgu133plus2")

target.list %>% getSYMBOL(., "hgu133plus2")


# **************************************** --------------------------------


# EX06:Consider Direction, T-test,no filter -------------------------------
ecelfile.set <- exprs ( Exprs.data )

#QN + log2

norm <- quantile_normalization(ecelfile.set = Exprs.data)

#two.sample t.test

adj.method <- "BH"

registerDoParallel(cores=8)
#function


t_test <- function(expression, adj.method = "BH", hypothesis = "greater"){
  gene.list <- rownames(expression)  
  t.result  <- foreach(i = 1:nrow(expression)) %dopar% {
    t.test(expression[i,1:3],expression[i,4:6],alternative = hypothesis, paired = FALSE,  conf.level = 0.95) %>%
      tidy %>% mutate(gene = gene.list[i])
  }
  t.result <- invoke(rbind, map(t.result,data.frame))
  
  t.result <- as.data.frame(t.result)
  rownames(t.result) <- t.result$gene
  t.result <- t.result %>%
    dplyr::select(-method, -alternative)%>%
    mutate(adjusted.p = p.adjust(p.value, method=adj.method),
           adjusted.method = adj.method)
  return(t.result)
}
#“holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”

#for CLS-CLF
# _greater ----------------------------------------------------------------
t.greater.CLF_CLS      <- t_test(expression = norm[,1:6], adj.method = "BH", hypothesis = "greater")
t.greater.Sphere_CLF   <- t_test(expression = norm[,c(7:9, 1:3)], adj.method = "BH", hypothesis = "greater")
t.greater.Sphere_CLS   <- t_test(expression = norm[,c(7:9, 4:6)], adj.method = "BH", hypothesis = "greater")

# _less -------------------------------------------------------------------
t.less.CLF_CLS      <- t_test(expression = norm[,1:6], adj.method = "BH", hypothesis = "less")
#the comparisom have change
t.less.Sphere_CLF   <- t_test(expression = norm[,c(7:9, 1:3)], adj.method = "BH", hypothesis = "less")
t.less.Sphere_CLS   <- t_test(expression = norm[,c(7:9, 4:6)], adj.method = "BH", hypothesis = "less")


# _write table ------------------------------------------------------------

write.table(t.greater.CLF_CLS,"greaterQNttestCLF_CLS.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(t.greater.Sphere_CLF,"greaterQNttestCLF_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(t.greater.Sphere_CLS,"greaterQNttestCLS_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)

write.table(t.less.CLF_CLS,"lessQNttestCLF_CLS.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(t.less.Sphere_CLF,"lessQNttestCLF_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(t.less.Sphere_CLS,"lessQNttestCLS_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)



# EX07:Consider Direction, T-test,filter by 0.5 ---------------------------
ecelfile.set <- exprs ( Exprs.data )

#QN + log2

norm <- quantile_normalization(ecelfile.set = Exprs.data)

#filter by variance 
Quantile_Normalization_eSet <-new("ExpressionSet", 
                                  phenoData = phenoData(celfile.set), 
                                  annotation   = annotation(celfile.set),
                                  protocolData = protocolData(celfile.set),
                                  experimentData = experimentData(celfile.set),
                                  exprs = norm )
QN_CLF_CLS_eSet     <- Quantile_Normalization_eSet[,1:6]
QN_Sphere_CLF_eSet  <- Quantile_Normalization_eSet[,c(7:9, 1:3)]
QN_Sphere_CLS_eSet  <- Quantile_Normalization_eSet[,c(7:9, 4:6)]

filterQN_CLF_CLS   <- prefilter_choose(QN_CLF_CLS_eSet , method = "filterbyvariance", parameters = c(0.5))
filterQN_Sphere_CLF<- prefilter_choose(QN_Sphere_CLF_eSet , method = "filterbyvariance", parameters = c(0.5))
filterQN_Sphere_CLS<- prefilter_choose(QN_Sphere_CLS_eSet , method = "filterbyvariance", parameters = c(0.5))


# _greater ----------------------------------------------------------------
f.t.greater.CLF_CLS      <- t_test(expression = exprs(filterQN_CLF_CLS), adj.method = "BH", hypothesis = "greater")
f.t.greater.Sphere_CLF   <- t_test(expression = exprs(filterQN_Sphere_CLF), adj.method = "BH", hypothesis = "greater")
f.t.greater.Sphere_CLS   <- t_test(expression = exprs(filterQN_Sphere_CLS), adj.method = "BH", hypothesis = "greater")

# _less -------------------------------------------------------------------
f.t.less.CLF_CLS      <- t_test(expression = exprs(filterQN_CLF_CLS), adj.method = "BH", hypothesis = "less")
f.t.less.Sphere_CLF   <- t_test(expression = exprs(filterQN_Sphere_CLF), adj.method = "BH", hypothesis = "less")
f.t.less.Sphere_CLS   <- t_test(expression = exprs(filterQN_Sphere_CLS), adj.method = "BH", hypothesis = "less")


# _write table ------------------------------------------------------------

write.table(f.t.greater.CLF_CLS,"greaterfilter50QNttestCLF_CLS.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(f.t.greater.Sphere_CLF,"greaterfilter50QNttestCLF_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(f.t.greater.Sphere_CLS,"greaterfilter50QNttestCLS_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)

write.table(f.t.less.CLF_CLS,"lessfilter50QNttestCLF_CLS.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(f.t.less.Sphere_CLF,"lessfilter50QNttestCLF_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)
write.table(f.t.less.Sphere_CLS,"lessfilter50QNttestCLS_Sphere.txt", quote=TRUE, sep="\t",row.names=TRUE, col.names = TRUE)



# EX07 MAS5, QN, LOG2, Characteristic Direction ---------------------------

GeoDE.norm <- norm
GeoDE.norm <- as.data.frame(GeoDE.norm)
GeoDE.norm$Probenames <- rownames(GeoDE.norm) 
GeoDE.norm.Sphere_CLS <- GeoDE.norm[,c(10,4:9)]
GeoDE.norm.CLF_CLS    <- GeoDE.norm[c(10,4:6,1:3)]
library(GeoDE)
# Load the example data
data(example_expression_data)
data(example_sampleclass)
data(example_gammas)
# Examine the expression data
head(example_expression_data)
# Examine the corresponding sample class factor
example_sampleclass

# Run the analysis
chdir_analysis_example <- chdirAnalysis(GeoDE.norm.Sphere_CLS,example_sampleclass,example_gammas
                                        ,CalculateSig=TRUE,nnull=10)
# Examine the results with the first value of the shrinkage parameter (gamma)
# show the first few of the most important genes.
lapply(chdir_analysis_example$results, function(x) x[1:10])
# We can also extract the results of the \code{chdirSig} function
# for example chdir_analysis_example$chdirprops[[1]] gives the whole
# characteristic direction vector for each value of gamma:
lapply(chdir_analysis_example$chdirprops[[1]],head)
# and the estimated number of significant genes can be recovered with
chdir_analysis_example$chdirprops$number_sig_genes


# EX08 MAS5, QN, LOG2, ANNOVA ---------------------------------------------
ecelfile.set <- exprs ( Exprs.data )

#QN + log2

norm <- quantile_normalization(ecelfile.set = ecelfile.set, method="median")

#annova
annova_test <- function(norm){
registerDoParallel(cores=16)
ex.design <- factor(c(rep("CLF",3), rep("CLS",3), rep("Sphere",3)))

probe.name <- rownames(norm)

test.annova.result  <- foreach(i = 1:nrow(norm)) %dopar% {
    aov(norm[i,] ~ ex.design) 
}



# _Post-hoc ---------------------------------------------------------------


# __TukeyHSD --------------------------------------------------------------



test.tukey.result  <- foreach(i = 1:nrow(norm)) %dopar% {
    TukeyHSD(test.annova.result[[i]]) %>%
        tidy %>%
        filter(comparison != "Sphere-CLF") %>%
        mutate(annova.p.value=(test.annova.result[[i]] %>%
                                   tidy() %>%
                                   select(p.value))[[1]][1], probe = probe.name[i])
}

test.tukey.result <- invoke(rbind, map(test.tukey.result,data.frame)) %>% select(-term)
return(test.tukey.result)
}

annova.tukey.result <- annova_test(norm)

# not use the doparrell
ex.design <- factor(c(rep("CLF",3), rep("CLS",3), rep("Sphere",3)))
probe.name <- rownames(norm)

test.annova.result <- list()
test.tukey.result  <- list()
Sys.time()
# Step 1
for (i in 1:length(norm)){
    test.annova.result[[i]] <- aov(norm[i,] ~ ex.design)
}
Sys.time()
save(test.annova.result, file="test_annova_result.Rdata")
# Step 2
for ( i in 1:length(norm)){
    test.tukey.result[[i]] <- TukeyHSD(test.annova.result[[i]])
    mem_used()
    print(i)
}

Sys.time()
# Step 3
test.tukey.result  <- foreach(i = 1:nrow(norm)) %dopar% {
  TukeyHSD(test.annova.result[[i]]) 
}
save(test.tukey.result, file="test_tukey_result.Rdata")
Sys.time()

# Step 4
test.tukey.result.dataframe  <- foreach(i = 1:nrow(norm)) %dopar% {
  test.tukey.result[[i]]  %>% with(ex.design) %>% tidy()
}

# Step 5
test.tukey.result.p.value  <- foreach(i = 1:nrow(norm)) %dopar% {
  test.annova.result[[i]] %>% tidy %>% slice_(1) %>% with(p.value)
}
rm(test.annova.result)
test.tukey.result.p.value <- unlist(test.tukey.result.p.value)

# Step 6 
test.tukey.result.dataframe.mutate  <- foreach(i = 1:nrow(norm)) %dopar% {
  test.tukey.result.dataframe[[i]]%>% mutate(annova.p.value = test.tukey.result.p.value[[i]],
                                             Probe = probe.name[[i]]) %>% rename(case=`.rownames`)
}

# Step 7
test.tukey.final.result <- invoke(rbind, map(test.tukey.result.dataframe.mutate,data.frame)) 
Sys.time()

# Step 8 annotate
test.tukey.final.result <- left_join(test.tukey.final.result, annotated.entrez.symbol, by="Probe")

# __Holm ------------------------------------------------------------------
registerDoParallel(cores=4)

system.time(test.aov.holm.result  <- foreach(i = 1:nrow(norm)) %dopar% {
    pairwise.t.test(norm[i,],ex.design,p.adjust="holm")
})


system.time(tidy.test.aov.holm.result         <- foreach(i = 1:nrow(norm)) %dopar% {
    test.aov.holm.result[[i]] %>% broom::tidy() %>% dplyr::slice(c(1,3))
})

system.time(tidy.aov.holm.result <- foreach( i = 1:nrow(norm)) %dopar% {
tidy.test.aov.holm.result[[i]] %>% mutate(Probe=probe.name[[i]], case=paste0(group1,"-",group2))
})

system.time(aov.holm.final <- invoke(rbind, map(tidy.aov.holm.result, data.frame)))


# __BH ------------------------------------------------------------------------
system.time(test.aov.bamjamin.result  <- foreach(i = 1:nrow(norm)) %dopar% {
    pairwise.t.test(norm[i,],ex.design,p.adjust="BH")
})


system.time(tidy.test.aov.bamjamin.result <- foreach(i = 1:nrow(norm)) %dopar% {
    test.aov.bamjamin.result[[i]] %>% broom::tidy() %>% dplyr::slice(c(1,3))
})

system.time(tidy.aov.bamjamin.result <- foreach(i = 1:nrow(norm)) %dopar% {
    tidy.test.aov.bamjamin.result[[i]] %>% mutate(Probe=probe.name[[i]], case=paste0(group1,"-",group2))
})

aov.bamjamin.final <- invoke(rbind, map(tidy.aov.bamjamin.result, data.frame))
# __Bonferroni ----------------------------------------------------------------

system.time(test.aov.bonferroni.result  <- foreach(i = 1:nrow(norm)) %dopar% {
    pairwise.t.test(norm[i,],ex.design,p.adjust="bonferroni")
})


tidy.test.aov.bonferroni.result <- foreach(i = 1:nrow(norm)) %dopar% {
    test.aov.bonferroni.result[[i]] %>% broom::tidy() %>% dplyr::slice(c(1,3))
}

system.time(tidy.aov.bonferroni.result <- foreach(i = 1:nrow(norm)) %dopar% {
    tidy.test.aov.bonferroni.result[[i]] %>% mutate(Probe=probe.name[[i]], case=paste0(group1,"-",group2))
}
)
system.time(aov.bonferroni.final <- invoke(rbind, map(tidy.aov.bonferroni.result, data.frame)))

save(aov.holm.final, aov.bamjamin.final, aov.bonferroni.final, file="pairttest_holm_bh_bf_dataframe.Rdata")
# Appendix: Annotation  -------------------------------------------------------------
# the library dependence
library(hgu133plus2.db)
library(readxl)
library(annotate)
# getSYMBOL(x, data)
# getLL(x, data)
# getEG(x, data)
# getGO(x, data)
# getPMID(x, data)
# getGOdesc(x, which)
# lookUp(x, data, what, load = FALSE)
# getUniqAnnItem()
# Load the Gene Ontology information  

hgu133plus2.probe <- rownames(ecelfile.set)
save(hgu133plus2.probe, file = "hgu133plus2.RData")
load("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/data/intermediate/hgu133plus2.RData")


hgu133plus2.probe.annotate <- gconvert(hgu133plus2.probe,
                                       organism = "hsapiens",
                                       target="AFFY_HG_U133_PLUS_2",
                                       mthreshold=1,
                                       filter_na = TRUE,
                                       df =T)
#turn upper to lower
hgu133plus2.probe.annotate <- hgu133plus2.probe.annotate %>% mutate(Probe=tolower(alias))

#add the information of entrezid
Symbol_Entrezid <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = hgu133plus2.probe.annotate$name ,
                                         columns = c("ENTREZID"),
                                         keytype = "SYMBOL")
colnames(Symbol_Entrezid) <- c("name","ENTREZID")
hgu133plus2.probe.annotate <- left_join(hgu133plus2.probe.annotate, Symbol_Entrezid, by="name") %>% unique

#filter some duplicated annotatione 
dup.prob <- hgu133plus2.probe.annotate.v1$Probe[duplicated(hgu133plus2.probe.annotate.v1$Probe)]
dup.annotated.probe <- hgu133plus2.probe.annotate.v1 %>% filter(Probe %in% dup.prob)
dup.symbol <- dup.annotated.probe$name %>% unique()
right.entrzid <- annotated.entrez.symbol %>% filter(Symbol %in% dup.symbol) %>% with(EntrezID)
not.used.entrezid <- dup.annotated.probe %>% filter(!(ENTREZID %in% right.entrzid)) %>% with(ENTREZID)

hgu133plus2.probe.annotate.v1 <- hgu133plus2.probe.annotate.v1 %>% filter(!(ENTREZID %in% not.used.entrezid ))
dup.prob <- hgu133plus2.probe.annotate.v1$Probe[duplicated(hgu133plus2.probe.annotate.v1$Probe)]


probe_6421 <-  annotated.entrez.symbol %>% filter(Symbol == "SFPQ", EntrezID == "6421") %>% with(Probe)
probe_654780 <-  annotated.entrez.symbol %>% filter(Symbol == "SFPQ", EntrezID == "654780") %>% with(Probe)

hgu133plus2.probe.annotate.v1$ENTREZID[hgu133plus2.probe.annotate.v1$Probe %in% probe_6421] <- "6421"
hgu133plus2.probe.annotate.v1$ENTREZID[hgu133plus2.probe.annotate.v1$Probe %in% probe_654780] <- "654780"

hgu133plus2.probe.annotate.v1 <- hgu133plus2.probe.annotate.v1 %>% unique

#save to RData
save(hgu133plus2.probe.annotate, file="hgu133plus2_annotation_v1.RData")


# unnotated.probe <- total.probe.data.frame[is.na(total.probe.data.frame$Gene),]
# unnotated.probe <- unnotated.probe$Probe
# getGO(unnotated.probe,"hgu133plus2")

probel.relatedTF <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/data/Affymetrix_reference/proberelateTF.txt",delim=" ")

cancer_stem_marker <- read_excel("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/data/CSCdbdata/marker_list.xls")
#57
cancer_stem_related_marker <- read_excel("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/data/CSCdbdata/annotation.xls")
#13,308
study.interest <- c("NANOG","IGF1R","SDF1","POU5F2","POU5F1P3","POU5F1P4","POU5F1B","HGF","POU5F1","IGF2","LIF","TGF","SOX2","THY1","JAK1",
                    "CXCR4","IGF2","THY-1","ABCG2","ABCG","ABCG1","ALDH1","CD133","SDF2",
                    "CD44","PROM1","MYC","BMI1","EZH2","KIT","PAR","IGFBP","IGFBP1","IGFBP2","GCSF",
                    "CHBD","ANGPT","NT","IL3","IL4","IL5","IL7","CCL18","CCL11","FGF2","EPHA1","TGFBR1","TGFBR2","LIFR","SMAD2","YAP1","TCF21")
Cancer.stem.marker <- unique(cancer_stem_marker$symbol)
Cancer.stem.related <- unique(cancer_stem_related_marker$symbol)

total.probe <- hgu133plus2.probe
total.probe.data.frame <- data.frame(Probe = total.probe)
total.probe.data.frame <- left_join(total.probe.data.frame, probel.relatedTF)
# get gene symbom annotation with hgu133plus2.db
#total.probe.data.frame <- total.probe.data.frame %>% mutate(Gene = getSYMBOL(Probe,"hgu133plus2"))
total.probe.data.frame <- total.probe.data.frame %>% left_join(.,dplyr::select(hgu133plus2.probe.annotate,Probe,name,description))

total.probe.data.frame$CSmarker[total.probe.data.frame$name %in% Cancer.stem.marker] <- 1
total.probe.data.frame$CSrelated[total.probe.data.frame$name %in% Cancer.stem.related] <- 1


save(total.probe.data.frame, file = "total_probe_dataframe.Rdata")



# Appendix: Annotation with PANTHER ---------------------------------------
load("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/data/intermediate/hgu133plus2.RData")

probe <- hgu133plus2.probe
get(probe[1], env=hgu133plus2ENTREZID)
entrezID <- unlist(lookUp(probe, "hgu133plus2.db", "ENTREZID"))
symbol <- getSYMBOL(probe, "hgu133plus2.db")
annotated.entrez.symbol <- data.frame(Probe = probe,
                                      EntrezID = entrezID,
                                      Symbol = symbol)

save(annotated.entrez.symbol, file= "/Users/Weitinglin/Documents/Repository/code_in_lab/annotated_entrez_symbol.Rdata")
# Post Filter Assessmnet --------------------------------------------------
nofilter <- exprs(QN_CLF_CLS_eSet )
nofilter <- as.data.frame(nofilter)
nm1 <- colnames(nofilter)

rowIQR <- function(x){
  
  d <- dim(x)[[1]]
  n <- rep(0,time=d)
  for (i in 1:d){
    n[[i]] <- IQR(x[i,])
  }
  return(n)
  
}


nofilter <- nofilter %>% mutate(IQR=rowIQR(as.matrix(.[nm1]))) 
filter <- nofilter[nofilter$IQR > 0.825841347,]
filter <- filter[filter$IQR > 0.825841347,]
quantile(nofilter$IQR)
length(nofilter$IQR)
filter   <- exprs(filterQN_CLF_CLS)
filter <- as.data.frame(filter)
filter <- filter %>% mutate(IQR=rowIQR(as.matrix(.[nm1]))) 
quantile(filter$IQR)
length(filter$IQR)

combine  <- bind_rows(nofilter, filter) 

test1 <- combine %>% filter(type == "nofilter" && Sample == "CHW_CLF 13-6-1.CEL")
test2 <-  combine %>% filter(type == "filter") %>% filter(Sample == "CHW_CLF 13-6-1.CEL") 

ggplot(data = combine, aes(x = Sample, y = value)) +
  geom_boxplot(aes(fill = Sample)) + 
  geom_jitter(width=0.2, alpha=0.2) + 
  facet_grid(type ~.)

ggplot(data = combine, aes(x = value)) +
  geom_density(aes(color = Sample)) + facet_grid(type ~.)




# ==================== *******************************************************----------------------------------------------------

# EX09 --------------------------------------------------------------------
#=====================input the data================================
data.path  <- file.path("/Users/Weitinglin/Documents/2016 實驗室資料處理/201510 microarray/raw data 20100114/set 9 CLS Sphere")
experiment.set <- c ( rep ("set9_CLS_P0" , 3 ) , rep ( "set9_CLS_Sphere" , 3 ) )


celfile.set <- do_phenodata(data.path = data.path, experiment.set = experiment.set)

#=====================preprocess the data ======
# use the mas5
Set9.Exprs.data <- mas5 ( celfile.set , normalize = FALSE, analysis = "absolute", sc = 500)

# save the temporary file
save(Set9.Exprs.data, file="Set9_Exprs_data.RData")

# expression profile
ecelfile.set <- exprs ( Set9.Exprs.data )

# quantile normalization
norm <- quantile_normalization(ecelfile.set = ecelfile.set,  method="median")

# t-test


registerDoParallel(cores=8)
#function

t_test <- function(tmp, hypothesis = "two.sides", adj.method = "BH"){

gene.list <- rownames(tmp)
hypothesis <-  "two.sided"

set.t.bamjamin.result  <- foreach(i = 1:nrow(tmp)) %dopar% {
    t.test(tmp[i,1:3], tmp[i,4:6], alternative = hypothesis, paired = FALSE,  conf.level = 0.95) 
}


tmp.tidy.set.t.bamjamin.result <- foreach(i = 1:nrow(tmp)) %dopar% {
    set.t.bamjamin.result[[i]] %>% broom::tidy() 
}

tidy.set.t.bamjamin.result <- foreach(i = 1:nrow(tmp)) %dopar% {
    tmp.tidy.set.t.bamjamin.result[[i]] %>%  mutate(gene = gene.list[i])
}

set.t.bamjamin.final.result<- invoke(rbind, map(tidy.set.t.bamjamin.result, data.frame))
set.t.bamjamin.final.result <- set.t.bamjamin.final.result %>% mutate(adjusted.p = p.adjust(p.value, method=adj.method),
                                    adjusted.method = adj.method)
return(set.t.bamjamin.final.result)
}
# EX10 --------------------------------------------------------------------

data.path  <- file.path("/Users/Weitinglin/Documents/2016 實驗室資料處理/201510 microarray/raw data 20100114/set 10 CLS Sphere  TimeSeries")
experiment.set <- c ( rep ("set10_CLS_P0" , 3 ) , rep ( "set10_CLS_Sphere" , 3 ) , rep ( "set10_CLS_Timeseies" , 3 ) )
setwd(data.path)

celfile.set <- do_phenodata(data.path = data.path, experiment.set = experiment.set)

#=====================preprocess the data ======
# use the mas5
Set10.Exprs.data <- mas5 ( celfile.set , normalize = FALSE, analysis = "absolute", sc = 500)

# save the temporary file
setwd("/Users/Weitinglin/Documents/Repository/code_in_lab")
save(Set10.Exprs.data, file="Set10_Exprs_data.RData")

# expression profile
ecelfile.set <- exprs ( Set10.Exprs.data )

# quantile normalization
norm <- quantile_normalization(ecelfile.set = ecelfile.set,  method="median")

# t.test
set10.t.result <- t_test(norm, adj.method = "BH", hypothesis = "two.sided")

# GSEA analysis -----------------------------------------------------------
load("~/Documents/Repository/code_in_lab/norm_qn_log.Rdata")
norm.tmp         <- norm[,c(1,2,3,4,5,6)]
norm.name        <- rownames(norm.tmp)
norm.description <- rep("na", nrow(norm.tmp))

norm.tmp <- as.data.frame(norm.tmp)
norm.tmp$NAME <- norm.name
norm.tmp$DESCRIPTION <- norm.description
colnames(norm.tmp) <- c("Sphere.1","Sphere.2","Sphere.3", "Sphere.1","Sphere.2","Sphere.3", "Name", "Description")
norm.tmp <- norm.tmp[,c(7,8,1:6)]
write.table(norm.tmp, file="set10_SphereCLS_norm_qn.gct", quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

# EX11 --------------------------------------------------------------------

data.path  <- file.path("/Users/Weitinglin/Documents/2016 實驗室資料處理/201510 microarray/raw data 20100114/set 11 CAF Sphere")
experiment.set <- c ( rep ("set11_Sphere" , 3 ) , rep ( "set11_CLF" , 3 ))
setwd(data.path)

celfile.set <- do_phenodata(data.path = data.path, experiment.set = experiment.set)

#=====================preprocess the data ======
# use the mas5
Set11.Exprs.data <- mas5 ( celfile.set , normalize = FALSE, analysis = "absolute", sc = 500)

# save the temporary file
setwd("/Users/Weitinglin/Documents/Repository/code_in_lab")
save(Set11.Exprs.data, file="Set11_Exprs_data.RData")

# expression profile
ecelfile.set <- exprs ( Set11.Exprs.data )

# quantile normalization
norm <- quantile_normalization(ecelfile.set = ecelfile.set,  method="median")

# t.test
set11.t.result <- t_test(norm, adj.method = "BH", hypothesis = "two.sided")




save(set11.t.result, file="Set11_ttestBH.Rdata")






# EX12 --------------------------------------------------------------------
data.path  <- file.path("/Users/Weitinglin/Documents/2016 實驗室資料處理/201510 microarray/raw data 20100114/set 12 Sphere  Timeseries")
experiment.set <- c ( rep ("set12_CLS_Sphere" , 3 ) , rep ( "set12_CLF" , 3 ))
setwd(data.path)

celfile.set <- do_phenodata(data.path = data.path, experiment.set = experiment.set)

#=====================preprocess the data ======
# use the mas5
Set11.Exprs.data <- mas5 ( celfile.set , normalize = FALSE, analysis = "absolute", sc = 500)

# save the temporary file
setwd("/Users/Weitinglin/Documents/Repository/code_in_lab")
save(Set11.Exprs.data, file="Set11_Exprs_data.RData")


