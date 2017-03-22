library(ReactomePA)
library(reactome.db)
library(dplyr)
library(topGO)
library(purrr)
library(magrittr)
library(readr)
library(stringr)
library(AnnotationDbi)



# Load the data -----------------------------------------------------------
load("/Users/Weitinglin/Documents/Repository/code_in_lab/updown_tukey_aov_result.Rdata")

## Upregulation
cutoff.p <- 0.05
up.CLF_CLS    <- up.tukey.aov.CLF_CLS %>% filter(BH.p.adj < cutoff.p)
up.Sphere_CLS <- up.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)


# Downregulation
cutoff.p <- 0.05
down.CLF_CLS    <- down.tukey.aov.CLF_CLS  %>% filter(BH.p.adj < cutoff.p)
down.Sphere_CLS <- down.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)





# _down.stream analysis ---------------------------------------------------
# Hierachial cluster use heatmaply
union.up       <- union(up.CLF_CLS$EntrezID, up.Sphere_CLS$EntrezID)
intersect.up   <- intersect(up.CLF_CLS$EntrezID, up.Sphere_CLS$EntrezID)

union.down     <- union(down.CLF_CLS$EntrezID, down.Sphere_CLS$EntrezID)
intersect.down <- intersect(down.CLF_CLS$EntrezID, down.Sphere_CLS$EntrezID)


unique(intersect.up,intersect.down) %>% length

# SPIA analysis:KEGG -----------------------------------------------------------
library(SPIA)
data(colorectalcancer)
res=spia(de=DE_Colorectal,
         all=ALL_Colorectal,
         organism="hsa"
         ,nB=2000,
         plots=FALSE,
         beta=NULL,
         combine="fisher")


calculateLog2foldchange <- function(x){
   set1 <- 2^x[,1:3]
   set2 <- 2^x[,4:6]
   probe.name <- rownames(x)
   set1.rowmean <- rowMeans(set1)
   set2.rowmean <- rowMeans(set2)
   foldchange <- set1.rowmean/set2.rowmean
   foldchange <- foldchange + 1
   log2foldchange <- log2(foldchange)
   result <- data.frame(Probe=probe.name, logFC=log2foldchange)
}

meanDuplicateEntrezID <- function(vector){
    SPIA.vector <- vector    
    tmp.data.frame <- data.frame(ENTREZID=names(SPIA.vector), logFC=unname(SPIA.vector))
    tmp.data.frame <- tmp.data.frame %>% group_by(ENTREZID) %>% summarise(logFC=mean(logFC))
    new.vector <- tmp.data.frame$logFC
    names(new.vector) <- tmp.data.frame$ENTREZID
    return(new.vector)
}

CLF_CLS.subset    <- norm[,c(1:6)] 
Sphere_CLS.subset <- norm[,c(7:9,4:6)]
load("~/Documents/Repository/code_in_lab/hgu133plus2_annotation_v1.RData")
Probe_EntrezID_Symbol <- hgu133plus2.probe.annotate[,c(5,8,9)]

tukey.aov.CLF_CLS    <- bind_rows(up.tukey.aov.CLF_CLS, down.tukey.aov.CLF_CLS)
tukey.aov.CLF_CLS    <- left_join(tukey.aov.CLF_CLS, Probe_EntrezID_Symbol, by="Probe")


tukey.aov.Sphere_CLS <- bind_rows(up.tukey.aov.Sphere_CLS, down.tukey.aov.Sphere_CLS)
tukey.aov.Sphere_CLS <- left_join(tukey.aov.Sphere_CLS, Probe_EntrezID_Symbol, by="Probe")




hash.table <- vector()
for ( i in 1:length(tukey.aov.CLF_CLS$Probe)){
    hash.table[[tukey.aov.CLF_CLS$Probe[i]]] <- tukey.aov.CLF_CLS$ENTREZID[i]
}


# CLF_CLS
cutoff.p <- 0.0001
CLF_CLS.log2change         <- calculateLog2foldchange(CLF_CLS.subset)
CLF_CLS.sig <- tukey.aov.CLF_CLS %>% filter(BH.p.adj < cutoff.p) %>% with(Probe)

SPIA.vector <- CLF_CLS.log2change$logFC
names(SPIA.vector) <- hash.table[CLF_CLS.log2change$Probe]
SPIA.vector <- SPIA.vector[!is.na(names(SPIA.vector))]
mean.SPIA.vector <- meanDuplicateEntrezID(SPIA.vector)
DE.SPIA.vector <-  mean.SPIA.vector[names(mean.SPIA.vector) %in% (hash.table[CLF_CLS.sig] %>% unname %>% unique)]


CLF_CLS.spia.result <- spia(de=DE.SPIA.vector,
                            all=names(mean.SPIA.vector),
                            organism = "hsa",
                            nB=2000,
                            plots = FALSE,
                            combine = "fisher")


# Sphere_CLS
cutoff.p <- 0.0001
Sphere_CLS.log2change   <- calculateLog2foldchange(Sphere_CLS.subset)
Sphere_CLS.sig <- tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p) %>% with(Probe)


SPIA.vector <- Sphere_CLS.log2change$logFC
names(SPIA.vector) <- hash.table[Sphere_CLS.log2change$Probe]
SPIA.vector <- SPIA.vector[!is.na(names(SPIA.vector))]
mean.SPIA.vector <- meanDuplicateEntrezID(SPIA.vector)
DE.SPIA.vector <-  mean.SPIA.vector[names(mean.SPIA.vector) %in% (hash.table[Sphere_CLS.sig] %>% unname %>% unique)]

CLS_Sphere.spia.result <- spia(de=DE.SPIA.vector,
                            all=names(mean.SPIA.vector),
                            organism = "hsa",
                            nB=2000,
                            plots = FALSE,
                            combine = "fisher")




# use time series data ----------------------------------------------------
load("~/Documents/Repository/code_in_lab/hgu133plus2_annotation_v1.RData")
load("/Users/Weitinglin/Documents/Repository/code_in_lab/updown_tukey_aov_result.Rdata")
probe.name <- rownames(norm)
norm.time.series <- as.data.frame(norm.time.series)
norm.time.series <- norm.time.series %>% mutate(p6=(`CHW_CLS 1-2 p.6-1.CEL`+`CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-1.CEL`)/3) %>% dplyr::select(c(1,2,6,7))
norm.time.series <- norm.time.series %>% rename(p14=`CHW_CLS 1-2 p.14.CEL`, p3=`CHW_CLS 1-2 p.3.CEL`, p8=`CHW_CLS 1-2 p.8.CEL`)
norm.time.series <- norm.time.series[,c(2,4,3,1)]
time.series.logFC <- norm.time.series %>% transmute(p6_p3logFC=p6-p3, p8_p6logFC=p8-p6, p14_p8logFC=p14-p8 )



getSPIAvector <- function(x,probe.name){
    SPIA.vector <- x[,1]
    names(SPIA.vector) <- hash.table[probe.name]
    SPIA.vector <- SPIA.vector[!is.na(names(SPIA.vector))]
    mean.SPIA.vector <- meanDuplicateEntrezID(SPIA.vector)
    return(mean.SPIA.vector)
}
cutoff.p <- 0.0001
Sphere_CLS.sig <- tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p) %>% with(Probe)

set_p6_p3 <- time.series.logFC %>% dplyr::select(p6_p3logFC)
set_p6_p3.SPIA.vector  <-  getSPIAvector(set_p6_p3)
set_p6_p3.DE.SPIA.vector <-  set_p6_p3.SPIA.vector[names(set_p6_p3.SPIA.vector) %in% (hash.table[Sphere_CLS.sig] %>% unname %>% unique)]

set_p8_p6 <- time.series.logFC %>% dplyr::select(p8_p6logFC)
set_p8_p6.SPIA.vector  <- getSPIAvector(set_p8_p6)  
set_p8_p6.DE.SPIA.vector <-  set_p8_p6.SPIA.vector[names(set_p8_p6.SPIA.vector) %in% (hash.table[Sphere_CLS.sig] %>% unname %>% unique)]

set_p14_p8 <- time.series.logFC %>% dplyr::select(p14_p8logFC)
set_p14_p8.SPIA.vector <- getSPIAvector(set_p14_p8)
set_p14_p8.DE.SPIA.vector <-  set_p14_p8.SPIA.vector[names(set_p14_p8.SPIA.vector) %in% (hash.table[Sphere_CLS.sig] %>% unname %>% unique)]


set_p6_p3.spia.result <- spia(de=set_p6_p3.DE.SPIA.vector,
                               all=names(set_p6_p3.SPIA.vector),
                               organism = "hsa",
                               nB=2000,
                               plots = FALSE,
                               combine = "fisher")

set_p8_p6.spia.result <- spia(de=set_p8_p6.DE.SPIA.vector,
                              all=names(set_p8_p6.SPIA.vector),
                              organism = "hsa",
                              nB=2000,
                              plots = FALSE,
                              combine = "fisher")

set_p14_p8.spia.result <- spia(de=set_p14_p8.DE.SPIA.vector,
                              all=names(set_p14_p8.SPIA.vector),
                              organism = "hsa",
                              nB=2000,
                              plots = FALSE,
                              combine = "fisher")
save(set_p6_p3.spia.result, set_p8_p6.spia.result, set_p14_p8.spia.result, file="time_series_spia.Rdata")
# Lactase -----------------------------------------------------------------

lactose_biosynthetic <- c("SLC2A1",	"B4GALT1", 	"LALBA")



# Reactome pathway analysis -----------------------------------------------


# Seperate
#*******************

# UP
up.reactome.CLF_CLS     <- enrichPathway(gene=up.CLF_CLS$EntrezID, pvalueCutoff = 0.05, readable = T)
up.reactome.Sphere_CLS  <- enrichPathway(gene=up.Sphere_CLS$EntrezID, pvalueCutoff = 0.05, readable = T)

# DOWN
down.reactome.CLF_CLS   <- enrichPathway(gene=down.CLF_CLS$EntrezID, pvalueCutoff = 0.05, readable = T)
down.reactome.Sphere_CLS<- enrichPathway(gene=down.Sphere_CLS$EntrezID, pvalueCutoff = 0.05, readable = T)



# Intersect
#*******************

# UP
intersect.up.reactome.pathway <- enrichPathway(gene=intersect.up, pvalueCutoff=0.05, readable=T)
intersect.up.reactome.result <- intersect.up.reactome.pathway@result %>% mutate(case = "Up", type = "intersect")

# DOWN
intersect.down.reactome.pathway <- enrichPathway(gene=intersect.down, pvalueCutoff=0.05, readable=T)
intersect.down.reactome.result <- intersect.down.reactome.pathway@result %>% mutate(case = "Down", type = "intersect")

# Union
#********************


# UP

uion.up.reactome.pathway <- enrichPathway(gene=union.up, pvalueCutoff=0.05, readable=T)
uion.up.reactome.result <- uion.up.reactome.pathway@result %>% mutate(case = "Up", type = "union")
# DOWN
uion.down.reactome.pathway <- enrichPathway(gene=union.down, pvalueCutoff=0.05, readable=T)
uion.down.reactome.result <- uion.down.reactome.pathway@result %>% mutate(case = "Down", type = "union")

reactome.result.total <- bind_rows(intersect.up.reactome.result,
                                  intersect.down.reactome.result,
                                  uion.up.reactome.result,
                                  uion.down.reactome.result)
write.csv(reactome.result.total, file = "reactom_result_tukey_10_5.csv")



# GO term analysis --------------------------------------------------------

# data preparation
all.entrez.hgu133plus2 <- annotated.entrez.symbol["EntrezID"] %>% filter(!is.na(EntrezID)) %>% unique()

up.CLF_CLS.p <- up.CLF_CLS %>% filter(EntrezID %in% intersect.up, !is.na(EntrezID))  %>% with(BH.p.adj)
up.Sphere_CLS.p <- up.Sphere_CLS %>% filter(EntrezID %in% intersect.up, !is.na(EntrezID))  %>% with(BH.p.adj)

tmp <- up.CLF_CLS.p
names(tmp) <- intersect.up[!is.na(intersect.up)]

total.up  <- up.tukey.aov.CLF_CLS %>% filter(!is.na(EntrezID)) %>% with(BH.p.adj)
total.up <-  up.tukey.aov.CLF_CLS %>% filter(!is.na(EntrezID)) %>% with(EntrezID)

sampleGOdata <- new("topGOdata",
                    description = "Simple test",
                    ontology = "BP",
                    allGenes = total.up,
                    geneSel = tmp,
                    nodeSize = 10)


# use GO.db
library(GO.db)
library(org.Hs.eg.db)


# _up:CLF_CLS --------------------------------------------------------------


# Prepare the data structure topGO needed
all.up.GO.CLF_CLS.entrezID  <-  up.tukey.aov.CLF_CLS %>% with(EntrezID)
all.up.GO.CLF_CLS <- AnnotationDbi::select(org.Hs.eg.db,
                                           keys = all.up.GO.CLF_CLS.entrezID [!is.na(all.up.GO.CLF_CLS.entrezID )] %>% as.character,
                                           columns = "GO",
                                           keytype = "ENTREZID")

# turn the data.frame into list-like data structure
all.up.GO.CLF_CLS.list <- all.up.GO.CLF_CLS %>% filter(EVIDENCE != "IEA") %>% split(.$ONTOLOGY)
all.up.GO.BP.CLF_CLS   <- all.up.GO.CLF_CLS.list$BP %>% split(.$ENTREZID)  %>% map(with, GO)
all.up.GO.BP.CLF_CLS.GO2geneID <- inverseList(all.up.GO.BP.CLF_CLS)

# get the genelist
geneNames  <- names(all.up.GO.BP.CLF_CLS)
geneList <- factor(as.integer(geneNames %in% (up.CLF_CLS$EntrezID[!is.na(up.CLF_CLS$EntrezID)] %>% as.character)))
names(geneList) <- geneNames

# construct topGOdata
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
               annot = annFUN.gene2GO, gene2GO = all.up.GO.BP.CLF_CLS)

# calculate the result
resultFisher <- runTest(GOdata , algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

# query the result
allRes  <- GenTable(GOdata, classicFisher= resultFisher,
                    classicKS = resultKS, elimKS = resultKS.elim,
                    orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)
allRes  <- allRes %>% transmute(GO.ID=GO.ID,
                     Term=Term,
                     Annotated=Annotated,
                     Significant=Significant,
                     Expected=Expected,
                     classicFisher=as.numeric(classicFisher),
                     classicKS=as.numeric(classicKS),
                     elimKS=as.numeric(elimKS))

allRes.up.CLF_CLS <- allRes %>% filter(classicFisher < 0.05)
write.csv(allRes.up.CLF_CLS, file = "topGO_top200_BP_up_CLF_CLS.csv")

# _up:Sphere_CLS ----------------------------------------------------------
# Prepare the data structure topGO needed
all.up.GO.Sphere_CLS.entrezID  <-  up.tukey.aov.Sphere_CLS %>% with(EntrezID)
all.up.GO.Sphere_CLS <- AnnotationDbi::select(org.Hs.eg.db,
                                           keys = all.up.GO.Sphere_CLS.entrezID [!is.na(all.up.GO.Sphere_CLS.entrezID )] %>% as.character,
                                           columns = "GO",
                                           keytype = "ENTREZID")

# turn the data.frame into list-like data structure
all.up.GO.Sphere_CLS.list <- all.up.GO.Sphere_CLS %>% filter(EVIDENCE != "IEA") %>% split(.$ONTOLOGY, .$ENTREZID)
all.up.GO.BP.Sphere_CLS   <- all.up.GO.Sphere_CLS.list$BP
all.up.GO.BP.Sphere_CLS   <- all.up.GO.BP.Sphere_CLS %>% split(.$ENTREZID) %>% map(with, GO)
all.up.GO.BP.Sphere_CLS.GO2geneID <- inverseList(all.up.GO.BP.Sphere_CLS)

# get the genelist
geneNames  <- names(all.up.GO.BP.Sphere_CLS)
geneList <- factor(as.integer(geneNames %in% (up.Sphere_CLS$EntrezID[!is.na(up.Sphere_CLS$EntrezID)] %>% as.character)))
names(geneList) <- geneNames

# construct topGOdata
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = all.up.GO.BP.Sphere_CLS)

# calculate the result
resultFisher <- runTest(GOdata , algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

# query the result
allRes  <- GenTable(GOdata, classicFisher= resultFisher,
                    classicKS = resultKS, elimKS = resultKS.elim,
                    orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)
allRes  <- allRes %>% transmute(GO.ID=GO.ID,
                                Term=Term,
                                Annotated=Annotated,
                                Significant=Significant,
                                Expected=Expected,
                                classicFisher=as.numeric(classicFisher),
                                classicKS=as.numeric(classicKS),
                                elimKS=as.numeric(elimKS))

allRes.up.Sphere_CLS <- allRes %>% filter(classicFisher < 0.05)
write.csv(allRes.up.Sphere_CLS, file = "topGO_top200_BP_up_Sphere_CLS.csv")

# _down:CLF_CLS -----------------------------------------------------------
# Prepare the data structure topGO needed
all.down.GO.CLF_CLS.entrezID  <-  down.tukey.aov.CLF_CLS %>% with(EntrezID)
all.down.GO.CLF_CLS <- AnnotationDbi::select(org.Hs.eg.db,
                                             keys = all.down.GO.CLF_CLS.entrezID [!is.na(all.down.GO.CLF_CLS.entrezID )] %>% as.character,
                                             columns = "GO",
                                             keytype = "ENTREZID")

# turn the data.frame into list-like data structure
all.down.GO.CLF_CLS.list <- all.down.GO.CLF_CLS %>% filter(EVIDENCE != "IEA") %>% split(.$ONTOLOGY, .$ENTREZID)
all.down.GO.BP.CLF_CLS   <- all.down.GO.CLF_CLS.list$BP
all.down.GO.BP.CLF_CLS   <- all.down.GO.BP.CLF_CLS %>% split(.$ENTREZID) %>% map(with, GO)
all.down.GO.BP.CLF_CLS.GO2geneID <- inverseList(all.down.GO.BP.CLF_CLS)

# get the genelist
geneNames  <- names(all.down.GO.BP.CLF_CLS)
geneList <- factor(as.integer(geneNames %in% (down.CLF_CLS$EntrezID[!is.na(down.CLF_CLS$EntrezID)] %>% as.character)))
names(geneList) <- geneNames

# construct topGOdata
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = all.down.GO.BP.CLF_CLS)

# calculate the result
resultFisher <- runTest(GOdata , algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

# query the result
allRes  <- GenTable(GOdata, classicFisher= resultFisher,
                    classicKS = resultKS, elimKS = resultKS.elim,
                    orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)
allRes  <- allRes %>% transmute(GO.ID=GO.ID,
                                Term=Term,
                                Annotated=Annotated,
                                Significant=Significant,
                                Expected=Expected,
                                classicFisher=as.numeric(classicFisher),
                                classicKS=as.numeric(classicKS),
                                elimKS=as.numeric(elimKS))

allRes.down.CLF_CLS <- allRes %>% filter(classicFisher < 0.05)
write.csv(allRes.down.CLF_CLS, file = "topGO_top200_BP_down_CLF_CLS.csv")
# _down:Sphere_CLS --------------------------------------------------------

# Prepare the data structure topGO needed
all.down.GO.Sphere_CLS.entrezID  <-  down.tukey.aov.Sphere_CLS %>% with(EntrezID)
all.down.GO.Sphere_CLS <- AnnotationDbi::select(org.Hs.eg.db,
                                                keys = all.down.GO.Sphere_CLS.entrezID [!is.na(all.down.GO.Sphere_CLS.entrezID )] %>% as.character,
                                                columns = "GO",
                                                keytype = "ENTREZID")

# turn the data.frame into list-like data structure
all.down.GO.Sphere_CLS.list <- all.down.GO.Sphere_CLS %>% filter(EVIDENCE != "IEA") %>% split(.$ONTOLOGY, .$ENTREZID)
all.down.GO.BP.Sphere_CLS   <- all.down.GO.Sphere_CLS.list$BP
all.down.GO.BP.Sphere_CLS   <- all.down.GO.BP.Sphere_CLS %>% split(.$ENTREZID) %>% map(with, GO)
all.down.GO.BP.Sphere_CLS.GO2geneID <- inverseList(all.down.GO.BP.Sphere_CLS)

# get the genelist
geneNames  <- names(all.down.GO.BP.Sphere_CLS)
geneList <- factor(as.integer(geneNames %in% (down.Sphere_CLS$EntrezID[!is.na(down.Sphere_CLS$EntrezID)] %>% as.character)))
names(geneList) <- geneNames

# construct topGOdata
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = all.down.GO.BP.Sphere_CLS)

# calculate the result
resultFisher <- runTest(GOdata , algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

# query the result
allRes  <- GenTable(GOdata, classicFisher= resultFisher,
                    classicKS = resultKS, elimKS = resultKS.elim,
                    orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)
allRes  <- allRes %>% transmute(GO.ID=GO.ID,
                                Term=Term,
                                Annotated=Annotated,
                                Significant=Significant,
                                Expected=Expected,
                                classicFisher=as.numeric(classicFisher),
                                classicKS=as.numeric(classicKS),
                                elimKS=as.numeric(elimKS))

allRes.down.Sphere_CLS <- allRes %>% filter(classicFisher < 0.05)
write.csv(allRes.down.Sphere_CLS, file = "topGO_top200_BP_down_Sphere_CLS.csv")



# ＿up:CLF_CLS/Sphere_CLS union --------------------------------------------
# Prepare the data structure topGO needed
all.union.up.entrezID <- union(up.tukey.aov.CLF_CLS$EntrezID[!is.na(up.tukey.aov.CLF_CLS$EntrezID)],up.tukey.aov.Sphere_CLS$EntrezID[!is.na(up.tukey.aov.Sphere_CLS$EntrezID)])

getEntrezIDToGOList <- function(EntrezID){

all.union.up.entrezID <- EntrezID
all.union.up.GO <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = all.union.up.entrezID %>% as.character,
                                        columns = "GO",
                                        keytype = "ENTREZID")

# turn the data.frame into list-like data structure
all.union.up.GO.list <- all.union.up.GO %>% filter(EVIDENCE != "IEA") %>% split(.$ONTOLOGY, .$ENTREZID)
all.union.up.GO.BP   <- all.union.up.GO.list$BP
all.union.up.GO.BP   <- all.union.up.GO.BP %>% split(.$ENTREZID) %>% map(with, GO)
return(all.union.up.GO.BP)
}
all.union.up.GO.BP.GO2geneID <- inverseList(all.union.up.GO.BP)

# get the genelist
getTopGO_Object <- function(BP.GO.list, sub.set.entrezID){

all.union.up.GO.BP <- BP.GO.list
union.up <- sub.set.entrezID 
    
geneNames  <- names(all.union.up.GO.BP)
geneList <- factor(as.integer(geneNames %in% (union.up %>% as.character)))
names(geneList) <- geneNames

# construct topGOdata
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = all.union.up.GO.BP)
return(GOdata)
}
# calculate the result

getGOanalysisResultTop200 <- function(GOdata){
resultFisher <- runTest(GOdata , algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

# query the result
allRes  <- GenTable(GOdata, classicFisher= resultFisher,
                    classicKS = resultKS, elimKS = resultKS.elim,
                    orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)
allRes  <- allRes %>% transmute(GO.ID=GO.ID,
                                Term=Term,
                                Annotated=Annotated,
                                Significant=Significant,
                                Expected=Expected,
                                classicFisher=as.numeric(classicFisher),
                                classicKS=as.numeric(classicKS),
                                elimKS=as.numeric(elimKS))

allRes.up.union <- allRes %>% filter(classicFisher < 0.05)
return(allRes.up.union)}

write.csv(allRes.up.union, file = "topGO_top200_BP_up_union.csv")

# _down:CLF_CLS/Sphere_CLS union ------------------------------------------
all.union.down.entrezID <- union(down.tukey.aov.CLF_CLS$EntrezID[!is.na(down.tukey.aov.CLF_CLS$EntrezID)],down.tukey.aov.Sphere_CLS$EntrezID[!is.na(down.tukey.aov.Sphere_CLS$EntrezID)])

#step 1
all.union.down.GO.BP <- getEntrezIDToGOList(EntrezID = all.union.down.entrezID)

#step 2
GOdata <- getTopGO_Object(BP.GO.list = all.union.down.GO.BP, sub.set.entrezID = union.down)

#step 3
allRes.down.union <- getGOanalysisResultTop200(GOdata = GOdata)
write.csv(allRes.down.union, file="topGO_top200_BP_down_union.csv")

# _up:CLF_CLS/Sphere_CLS intersect ----------------------------------------
all.intersect.up.entrezID <- intersect(up.tukey.aov.CLF_CLS$EntrezID[!is.na(up.tukey.aov.CLF_CLS$EntrezID)],up.tukey.aov.Sphere_CLS$EntrezID[!is.na(up.tukey.aov.Sphere_CLS$EntrezID)])

#step 1
all.intersect.up.GO.BP <- getEntrezIDToGOList(EntrezID = all.intersect.up.entrezID)

#step 2
GOdata <- getTopGO_Object(BP.GO.list = all.intersect.up.GO.BP, sub.set.entrezID = intersect.up)

#step 3
allRes.up.intersect <- getGOanalysisResultTop200(GOdata = GOdata)
write.csv(allRes.up.intersect, file="topGO_top200_BP_up_intersect.csv")

# _down:CLF_CLS/Sphere_CLS intersect --------------------------------------
all.intersect.down.entrezID <- intersect(down.tukey.aov.CLF_CLS$EntrezID[!is.na(down.tukey.aov.CLF_CLS$EntrezID)],down.tukey.aov.Sphere_CLS$EntrezID[!is.na(down.tukey.aov.Sphere_CLS$EntrezID)])

#step 1
all.intersect.down.GO.BP <- getEntrezIDToGOList(EntrezID = all.intersect.down.entrezID)

#step 2
GOdata <- getTopGO_Object(BP.GO.list = all.intersect.down.GO.BP, sub.set.entrezID = intersect.down)

#step 3
allRes.down.intersect <- getGOanalysisResultTop200(GOdata = GOdata)
write.csv(allRes.down.intersect, file="topGO_top200_BP_down_intersect.csv")



# use the pathview  -------------------------------------------------------

library(pathview)
data("gse16873.d")
pv.out <- pathview(gene.data = gse16873.d[,1], pathway.id = "04110", species = "hsa", out.suffix = "gse16873")

data("demo.paths")
i <- 1

#classic version
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873", kegg.native = T)

list.files(pattern="hsa04110", full.names=T)
str(pv.out)

head(pv.out$plot.data.gene)

# double the output size in order to get less computing time
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.2layer", kegg.native = T,
                   same.layer = F)

# graphviz style
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species =  "hsa", out.suffix = "gse16873", kegg.native = F, 
                   sign.pos = demo.paths$spos[i] )
dim(pv.out$plot.data.gene)

head(pv.out$plot.data.gene)



# use g:profiler to get kegg pathway and visualization with pathvi --------
library(gProfileR)
gprofiler(query,
          organism = "hsapiens",
          ordered_query = F,
          significant =T,
          exclude_iea = F,
          underrep = F,
          evcodes = F,
          region_query =F,
          max_p_value = 1,
          min_set_size = 0,
          max_set_size = 0,
          min_isect_size = 0,
          correction_method = "analytical",
          hier_filtering = "none",
          domain_size = "annotated",
          custom_bg = "",
          numeric_ns = "",
          png_fn = NULL,
          include_graph = F,
          src_filter = NULL)


# up:CLF_CLS --------------------------------------------------------------
up.query <- up.CLF_CLS$Symbol[!is.na(up.CLF_CLS$Symbol)] %>% unique()
up.query <- up.query %>% as.vector()
query.background <- up.tukey.aov.CLF_CLS$Symbol[!is.na(up.tukey.aov.CLF_CLS$Symbol)] %>% unique()
query.background <- query.backgroud  %>% as.vector()
query.background <- query.backgroud  %>% as.character()
big.background   <- hgu133plus2.probe.annotate$name[!is.na(hgu133plus2.probe.annotate$name)] %>% unique() %>% as.character()
up.gprofiler.CLF_CLS.kegg <- gprofiler(query = up.query,
                                       organism = "hsapiens",
                                       ordered_query = F,
                                       significant =T,
                                       exclude_iea = T,
                                       underrep = F,
                                       evcodes = F,
                                       region_query =F,
                                       max_p_value = 10,
                                       # min_set_size = 0,
                                       # max_set_size = 1000,
                                       # min_isect_size = 0,
                                       correction_method = "fdr",
                                       hier_filtering = "moderate", #"none", "moderate", "strong"
                                       domain_size = "annotated",
                                       custom_bg = query.background ,
                                       numeric_ns = "",
                                       png_fn = NULL,
                                       include_graph = F,
                                       src_filter = c("KEGG"))




# Create Network with Reactome -----------------------------------------------------------
reactome_005_up_result   <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/pathway_reactome/result_up.csv")  
reactome_005_up_map      <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/pathway_reactome/mapping_up.csv") 
reactome_005_down_result <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/pathway_reactome/result_down.csv") 
reactome_005_down_map    <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/pathway_reactome/mapping_down.csv") 


# construct the edge information

reverse_aggregate_reactome <- function(result){

tmp.one <- result$`Submitted entities found` %>% str_split(";")
tmp.two <- result$`Pathway name` 
tmp.thr <- data.frame("Pathway name" = rep(tmp.two, sapply(tmp.one,length)), "EntrezID" = unlist(tmp.one), check.names = FALSE)
tmp.new <- left_join(tmp.thr, reactome_005_up_result, by="Pathway name")
return(tmp.new)
}

reactome_005_up_result   <- reverse_aggregate_reactome(result=reactome_005_up_result)
reactome_005_down_result <- reverse_aggregate_reactome(result=reactome_005_down_result)


# construct the node information
load("~/Documents/Repository/code_in_lab/hgu133plus2_annotation.RData")
library(org.Hs.eg.db)
library(AnnotationDbi)
hgu133plus2.probe <- hgu133plus2.probe.annotate[,c(8,5,6)]
Symbol_Entrezid <- AnnotationDbi::select(org.Hs.eg.db,
                      keys = hgu133plus2.probe$name ,
                      columns = c("ENTREZID"),
                      keytype = "SYMBOL")
colnames(Symbol_Entrezid) <- c("name","ENTREZID")
hgu133plus2.probe.annotate <- left_join(hgu133plus2.probe.annotate, Symbol_Entrezid, by="name")
hgu133plus2.probe <- left_join(hgu133plus2.probe, Symbol_Entrezid, by="name")
hgu133plus2.probe <- unique(hgu133plus2.probe)
colnames(hgu133plus2.probe) <- c("Probe", "Name", "Description", "EntrezID")

colnames(reactome_005_up_map)   <- c("EntrezID","ReactomeIdentifier","Resource")
colnames(reactome_005_down_map) <- c("EntrezID","ReactomeIdentifier","Resource")
reactome_005_up_map$EntrezID <- as.character(reactome_005_up_map$EntrezID)
reactome_005_down_map$EntrezID <- as.character(reactome_005_down_map$EntrezID)

reactome_005_up_map <- left_join(reactome_005_up_map, hgu133plus2.probe, by="EntrezID")
reactome_005_up_map <- reactome_005_up_map %>% dplyr::select(-c(Probe, ReactomeIdentifier))
reactome_005_down_map <- left_join(reactome_005_down_map, hgu133plus2.probe, by="EntrezID")
reactome_005_down_map <- reactome_005_down_map %>% dplyr::select(-c(Probe, ReactomeIdentifier))


# Create Network from the data --------------------------------------------


# Load the data -----------------------------------------------------------


## Upregulation
cutoff.p <- 0.05
up.CLF_CLS    <- up.tukey.aov.CLF_CLS %>% filter(BH.p.adj < cutoff.p)
up.Sphere_CLS <- up.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)


# Downregulation
cutoff.p <- 0.05
down.CLF_CLS    <- down.tukey.aov.CLF_CLS  %>% filter(BH.p.adj < cutoff.p)
down.Sphere_CLS <- down.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)





# _down.stream analysis ---------------------------------------------------
# Hierachial cluster use heatmaply
intersect.up   <- intersect(up.CLF_CLS$EntrezID, up.Sphere_CLS$EntrezID)
intersect.down <- intersect(down.CLF_CLS$EntrezID, down.Sphere_CLS$EntrezID)




#panthner
library(PANTHER.db)
pthOrganisms(PANTHER.db) <- "HUMAN"
search_panther_database <- function(EntrezID){
    
class.term <- AnnotationDbi::select(PANTHER.db,
                                    keys=unique(EntrezID),
                                    columns=c("CLASS_TERM"),
                                    keytype = "ENTREZ") 
class.term <- class.term %>% dplyr::rename(PANTHER_TERM=CLASS_TERM) %>% mutate(TERM_TYPE="CLASS_TERM")
component.term <- AnnotationDbi::select(PANTHER.db,
                                    keys=unique(EntrezID),
                                    columns=c("COMPONENT_TERM"),
                                    keytype = "ENTREZ")
component.term <- component.term %>% dplyr::rename(PANTHER_TERM=COMPONENT_TERM) %>% mutate(TERM_TYPE="COMPONENT_TERM")
family.term <- AnnotationDbi::select(PANTHER.db,
                                     keys=unique(EntrezID),
                                    columns=c("FAMILY_TERM"),
                                    keytype = "ENTREZ")
family.term <- family.term %>% dplyr::rename(PANTHER_TERM=FAMILY_TERM) %>% mutate(TERM_TYPE="FAMILY_TERM")
pathway.term <- AnnotationDbi::select(PANTHER.db,
                                     keys=unique(EntrezID),
                                     columns=c("PATHWAY_TERM"),
                                     keytype = "ENTREZ")
pathway.term <- pathway.term %>% dplyr::rename(PANTHER_TERM=PATHWAY_TERM) %>% mutate(TERM_TYPE="PATHWAY_TERM")
total <- bind_rows(class.term, component.term) %>%
         bind_rows(.,family.term) %>%
         bind_rows(.,pathway.term)
total <- total %>% filter(!is.na(PANTHER_TERM))
return(total)
}

hgu133plus2.005.up.panther   <-  search_panther_database(EntrezID = intersect.up) 
node.up.term   <- hgu133plus2.005.up.panther %>% mutate(REGULATION_TYPE="UP")

hgu133plus2.005.down.panther <-  search_panther_database(EntrezID = intersect.down)
node.down.term   <- hgu133plus2.005.down.panther %>% mutate(REGULATION_TYPE="DOWN")

node.up.Symbol_Entrezid <- AnnotationDbi::select(org.Hs.eg.db,
                                            keys = unique(node.up.term$ENTREZ) ,
                                            columns = c("SYMBOL"),
                                            keytype = "ENTREZID")
node.up.Symbol_Entrezid <- node.up.Symbol_Entrezid %>% mutate(EXPRESSION_TYPE="UP")

node.down.Symbol_Entrezid <- AnnotationDbi::select(org.Hs.eg.db,
                                                 keys = unique(node.down.term$ENTREZ) ,
                                                 columns = c("SYMBOL"),
                                                 keytype = "ENTREZID")
node.down.Symbol_Entrezid <- node.down.Symbol_Entrezid %>% mutate(EXPRESSION_TYPE="DOWN")

total.node <- bind_rows(node.up.Symbol_Entrezid, node.down.Symbol_Entrezid)

write_delim(total.node, path = "/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/first/reactom_node.txt", delim = "\t")
write_delim(total.edge.node, path = "/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/first/reactom_edge.txt", delim = "\t")
# create association attribute

TF.entrez <- bind_rows(node.up.term, node.down.term) %>% filter(str_detect(PANTHER_TERM,"transcription factor"), TERM_TYPE == "CLASS_TERM") 


TF.Symbol_Entrezid <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = TF.entrez$ENTREZ ,
                                         columns = c("SYMBOL"),
                                         keytype = "ENTREZID")
TF.Symbol_Entrezid <- TF.Symbol_Entrezid %>% unique


# a gene symbol can have multiple uniprot

# Query with PISCQUIC -----------------------------------------------------
library(PSICQUIC)
psicquic <- PSICQUIC()
providers(psicquic)


SearchPSICQUIC <- function(symbol_list){
    total.result <- data.frame()
    for (item in symbol_list){
      one.result <-   interactions(psicquic, id=item, species="9606", provider=c("BioGrid", "IntAct", "Reactome-FIs"))
      total.result <- bind_rows(total.result,one.result[,c(1,2,7,9,12,13,14)])
      }
    return(total.result)
}

TF.interaction.protein <- SearchPSICQUIC(TF.Symbol_Entrezid$SYMBOL)


#translate the id
entrez.TF.interaction.protein     <- TF.interaction.protein %>% filter(str_detect(A,"entrez"))
uniprot.TF.interaction.protein    <- TF.interaction.protein %>% filter(str_detect(A, "uniprotkb"))
not.entuni.TF.interaction.protein <- TF.interaction.protein %>% filter(!str_detect(A,"entrez"),!str_detect(A, "uniprotkb"))


TurnEnrezUniprotTextOff <- function(data, type){

    if (type == "entrez"){
    total.entrez <- c(data$A, data$B)
    total.entrez <- total.entrez %>% unique %>% str_replace_all("entrez gene\\/locuslink:","")
    }
    if (type == "uniprot"){
    total.entrez <- c(data$A, data$B)
    total.entrez <- total.entrez %>% unique %>% str_replace_all("uniprotkb:","")
    }
    return(total.entrez)
}

id.entrez.TF  <- TurnEnrezUniprotTextOff(data=entrez.TF.interaction.protein , type="entrez")
id.uniprot.TF <- TurnEnrezUniprotTextOff(data=uniprot.TF.interaction.protein , type="uniprot")


#use org.HG
id.symbol.uniprot <- AnnotationDbi::select(org.Hs.eg.db,
                            keys = id.uniprot.TF,
                            columns = c("SYMBOL"),
                            keytype = "UNIPROT")

id.symbol.from.entrez  <- AnnotationDbi::select(org.Hs.eg.db,
                                                 keys = id.entrez.TF,
                                                 columns = c("SYMBOL"),
                                                 keytype = "ENTREZID")



id.symbol.uniprot <- id.symbol.uniprot %>% filter(!is.na(SYMBOL))
#remove prefix
uniprot.TF.interaction.protein <- uniprot.TF.interaction.protein %>% mutate(A=str_replace_all(A,"uniprotkb:",""),B=str_replace_all(B,"uniprotkb:","")) 
entrez.TF.interaction.protein <- entrez.TF.interaction.protein %>% mutate(A=str_replace_all(A,"entrez gene\\/locuslink:",""),B=str_replace_all(B,"entrez gene\\/locuslink:","")) 


#filter without annotation id
uniprot.symbol <- vector()

for ( i in 1:nrow(id.entrez.from.uniprot)){
    uniprot.symbol[id.symbol.uniprot$UNIPROT[i]] <- id.symbol.uniprot$SYMBOL[i]    
}

uniprot.TF.interaction.protein <- uniprot.TF.interaction.protein %>% mutate(A=unname(uniprot.symbol[A]), B=unname(uniprot.symbol[B])) 
uniprot.TF.interaction.protein <- uniprot.TF.interaction.protein %>% filter(!is.na(A),!is.na(B))

# annotation the entrez with symbol
entrez.symbol <- vector()
for ( i in 1:nrow(id.entrez.from.entrez)){
    entrez.symbol[id.symbol.from.entrez$ENTREZID[i]] <- id.symbol.from.entrez$SYMBOL[i]    
}

entrez.TF.interaction.protein <- entrez.TF.interaction.protein %>% mutate(A=unname(entrez.symbol[A]), B=unname(entrez.symbol[B])) 
entrez.TF.interaction.protein <- entrez.TF.interaction.protein %>% filter(!is.na(A),!is.na(B))

total.clean.TF.interaction.ptotein <- bind_rows(entrez.TF.interaction.protein, uniprot.TF.interaction.protein)
# gather the total protein list
total.005.entrez.updown <- c(node.up.term$ENTREZ,node.down.term$ENTREZ) %>% unique
total.005.symbol.from.entrez  <- AnnotationDbi::select(org.Hs.eg.db,
                                                keys = total.005.entrez.updown,
                                                columns = c("SYMBOL"),
                                                keytype = "ENTREZID")





# edge information
tmp <- total.clean.TF.interaction.ptotein %>% filter(A %in% total.005.symbol.from.entrez$SYMBOL & B %in% total.005.symbol.from.entrez$SYMBOL)
tmp.1 <- tmp %>% filter(A %in% TF.Symbol_Entrezid$SYMBOL)
tmp.2 <- tmp %>% filter(B %in% TF.Symbol_Entrezid$SYMBOL & !(A %in% TF.Symbol_Entrezid$SYMBOL))
tmp.2 <- tmp.2 %>% dplyr::rename(A=B,B=A)
tmp.3 <- tmp %>% filter(A %in% TF.Symbol_Entrezid$SYMBOL & B %in% TF.Symbol_Entrezid$SYMBOL)

total.interaction.edge.inform <- bind_rows(tmp.1,tmp.2,tmp.3)

# node information
total.node <- total.node %>% mutate(TF="") 
total.node$TF[total.node$ENTREZID %in% TF.Symbol_Entrezid$ENTREZID ] <- "TF"
total.node$TF[!(total.node$ENTREZID %in% TF.Symbol_Entrezid$ENTREZID )] <- "Not TF"

#MI definiation
#http://irefindex.org/wiki/index.php?title=Mapping_of_terms_to_MI_term_ids_-_iRefIndex_6.0

tmp.interaction <- total.interaction.edge.inform
tmp.interaction <- tmp.interaction[,c(1,2,5)]
tmp.interaction <- unique(tmp.interaction)

#tmp.interaction <- tmp.interaction %>% filter(str_detect(type,"direct interaction") | str_detect(type,"physical interaction") |str_detect(type, "association"))
tmp.interaction <- tmp.interaction %>% filter(A %in% total.node$SYMBOL) %>% filter(B %in% total.node$SYMBOL)

#filter some interaction
# 1. self interaction



write_delim(tmp.interaction, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/second/interaction_edge.txt", delim = "\t")


#tag with the finding p value

TagPvalue <- function(pvalue_list){
    empty.list <- list()
    
    for (item in pvalue_list){
        ## Upregulation
            cutoff.p <- item
            up.CLF_CLS    <- up.tukey.aov.CLF_CLS %>% filter(BH.p.adj < cutoff.p)
            up.Sphere_CLS <- up.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)
            # Downregulation
            down.CLF_CLS    <- down.tukey.aov.CLF_CLS  %>% filter(BH.p.adj < cutoff.p)
            down.Sphere_CLS <- down.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)

            # Hierachial cluster use heatmaply
            intersect.up  <- intersect(up.CLF_CLS$EntrezID, up.Sphere_CLS$EntrezID)
            intersect.down <- intersect(down.CLF_CLS$EntrezID, down.Sphere_CLS$EntrezID)
            empty.list[as.character(item)] <- list(c(unique(c(intersect.up, intersect.down))))
            print(head(empty.list))
     }
  return(empty.list)
}


tag.p.list <- TagPvalue(c(0.05,0.01,0.001,0.0001))

tag.p.value.data.frame <- data.frame(ENTREZID = tag.p.list$`0.05`)
tag.p.value.data.frame <- tag.p.value.data.frame %>% mutate(Pvalue="0.05")
tag.p.value.data.frame$Pvalue[tag.p.value.data.frame$ENTREZID %in% tag.p.list$`0.01`] <- "0.01"
tag.p.value.data.frame$Pvalue[tag.p.value.data.frame$ENTREZID %in% tag.p.list$`0.001`] <- "0.001"
tag.p.value.data.frame$Pvalue[tag.p.value.data.frame$ENTREZID %in% tag.p.list$`1e-04`] <- "0.0001"
tag.p.value.data.frame <- tag.p.value.data.frame[!is.na(tag.p.value.data.frame$ENTREZID),]
total.node.tag.p.value <- left_join(total.node, tag.p.value.data.frame, by="ENTREZID")


write_delim(total.node.tag.p.value, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/second/interaction_node.txt", delim = "\t")
total.node <- read_delim("/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/second/interaction_node.txt", delim = "\t")


# Add the information GO/PANTHER ------------------------------------------
# add more information on the node (GO/PANTHER)

total.entrezid <- total.node.tag.p.value$ENTREZID
library(PANTHER.db)
pthOrganisms(PANTHER.db) <- "HUMAN"
total.entrezid.panther <- search_panther_database(EntrezID = total.entrezid) 


total.entrezid.up    <- total.node.tag.p.value %>% filter(EXPRESSION_TYPE=="UP") %>% with(ENTREZID)
total.entrezid.down  <- total.node.tag.p.value %>% filter(EXPRESSION_TYPE=="DOWN") %>% with(ENTREZID)

total.entrezid.up.GO <- AnnotationDbi::select(org.Hs.eg.db,
                      keys = total.entrezid.up,
                      columns = c("GO","ONTOLOGY","EVIDENCE"),
                      keytype = "ENTREZID")



total.entrezid.down.GO <- AnnotationDbi::select(org.Hs.eg.db,
                                              keys = total.entrezid.down,
                                              columns = c("GO","ONTOLOGY","EVIDENCE"),
                                              keytype = "ENTREZID")

total.entrezid.up.GO <- total.entrezid.up.GO %>% filter(EVIDENCE != "IEA")
total.entrezid.down.GO <- total.entrezid.down.GO %>% filter(EVIDENCE != "IEA")

total.GOID <- unique(c(total.entrezid.up.GO$GO,total.entrezid.down.GO$GO))


library(GO.db)
columns(GO.db)
keytypes(GO.db)

total.GO.annotation <- AnnotationDbi::select(GO.db,
                      keys = total.GOID,
                      columns = c("GOID","DEFINITION","ONTOLOGY","TERM"),
                      keytypes = "GOID")



# Projection the network into the pathway or GO term ----------------------

# ECM related
ECM.GO.related.term <- c("GO:0030198","GO:0043062","GO:0030199","GO:0061448","GO:0022617","GO:0032963","GO:0044259","GO:0006898","GO:0000904")
total.entrezid.up.related.ECM   <- total.entrezid.up.GO %>% filter(GO %in% ECM.GO.related.term) %>% with(ENTREZID)
total.entrezid.down.related.ECM <- total.entrezid.down.GO %>% filter(GO %in% ECM.GO.related.term) %>% with(ENTREZID)

SYMBOL.related.ECM <- total.node.tag.p.value %>% filter(ENTREZID %in% c(total.entrezid.up.related.ECM, total.entrezid.down.related.ECM)) %>% with(SYMBOL)
sub.edge.related.ECM <- tmp.interaction %>% filter(A %in% SYMBOL.related.ECM | B %in% SYMBOL.related.ECM)
sub.node.related.ECM <- total.node.tag.p.value %>% filter(SYMBOL %in% unique(c(sub.edge.related.ECM$A, sub.edge.related.ECM$B)))

write_delim(sub.node.related.ECM, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/ECM_related_subnetwork/interaction_node.txt", delim = "\t")
write_delim(sub.edge.related.ECM, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/ECM_related_subnetwork/interaction_edge.txt", delim = "\t")


# 
Metabolism.GO.related.term <- c("GO:0043604", "GO:0043603", "GO:0043603", "GO:0043043", "GO:0006518", "GO:0031396", "GO:0016567", "GO:1900040", "GO:0043039", "GO:1903322" ,"GO:0045184", "GO:0051248")
total.entrezid.up.related.Metabolism   <- total.entrezid.up.GO %>% filter(GO %in% Metabolism.GO.related.term) %>% with(ENTREZID)
total.entrezid.down.related.Metabolism <- total.entrezid.down.GO %>% filter(GO %in% Metabolism.GO.related.term) %>% with(ENTREZID)

SYMBOL.related.Metabolism <- total.node.tag.p.value %>% filter(ENTREZID %in% c(total.entrezid.up.related.Metabolism, total.entrezid.down.related.Metabolism)) %>% with(SYMBOL)
sub.edge.related.Metabolism <- tmp.interaction %>% filter(A %in% SYMBOL.related.Metabolism | B %in% SYMBOL.related.Metabolism)
sub.node.related.Metabolism <- total.node.tag.p.value %>% filter(SYMBOL %in% unique(c(sub.edge.related.Metabolism$A, sub.edge.related.Metabolism$B)))

write_delim(sub.node.related.Metabolism, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/Metabolism_related_subnetwork/interaction_node.txt", delim = "\t")
write_delim(sub.edge.related.Metabolism, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/Metabolism_related_subnetwork/interaction_edge.txt", delim = "\t")
# Metabolism related





