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


# GSEA analysis -----------------------------------------------------------
load("~/Documents/Repository/code_in_lab/norm_qn_log.Rdata")
norm.tmp         <- norm
norm.name        <- rownames(norm.tmp)
norm.description <- rep("na", nrow(norm.tmp))

norm.tmp <- as.data.frame(norm.tmp)
norm.tmp$NAME <- norm.name
norm.tmp$DESCRIPTION <- norm.description
colnames(norm.tmp) <- c("CLF.1","CLF.2", "CLF.3", "CLS.1","CLS.2","CLS.3","Sphere.1","Sphere.2","Sphere.3", "Name", "Description")
norm.tmp <- norm.tmp[,c(10,11,1,2,3,4,5,6,7,8,9)]
write.table(norm.tmp, file="norm_qn.gct", quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

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

quantile_normalization <- function(ecelfile.set, method="median"){
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

# Custom the SPIA pathway for analysis ------------------------------------


# download xml
hs.kegg.pathway <- read_delim("http://rest.kegg.jp/list/pathway/hsa", delim = "\t", col_names = FALSE)
colnames(hs.kegg.pathway) <- c("KEGG_Pathway_ID","Name")

index <- c(301,286,283,282,215,213,212,193,186,183,182,181,180,179,174,173,172,171,170,169,168,163,161,160,159,157,149,148,147,146,145,144,143,142,141,137,135,134,132,129,128,127,126,125,124,123,122,121,120,119,118,117,116,112,111,110,109,108,107,106,105,104,103,102,101,100)
id    <- hs.kegg.pathway$KEGG_Pathway_ID[index] %>% str_replace_all("path:","")
name  <-  hs.kegg.pathway$Name[index] %>% str_replace_all("[:space:]|\\/|\\(.+\\)","")

KMLname <- list(as.list(id),as.list(name)) %>% pmap(paste, sep="_") %>% unlist
KMLdata <- as.list(hs.kegg.pathway$KEGG_Pathway_ID[index] %>% str_replace_all("path:","") )%>%map(~ paste0("http://rest.kegg.jp/get/",.,"/kgml")) %>% unlist

fs <- c()
tmpdir <- tempdir()
setwd(tempdir())
setwd(getwd())

for (i in 1:length(KMLdata)) {
    path <- paste0(KMLname[i], ".xml")
    fs <- c(fs, path)
    kegg.api.result <- read_xml(KMLdata[i])
    write_xml(kegg.api.result, path, options = "format")
}
zip(zipfile="output.zip", files=fs)

kml.dir <- "/Users/Weitinglin/Documents/R_scripts/Lab/kegg/output"
kml.spia.dir <- "/Users/Weitinglin/Documents/R_scripts/Lab/kegg/spiaXML/"

# just parse the signal pathway, not metabolite related pathway
makeSPIAdata(kgml.path=kml.dir, organism="hsa", out.path=kml.spia.dir)

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
                            combine = "fisher",
                            data.dir = kml.spia.dir)


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
                            combine = "fisher",
                            data.dir = kml.spia.dir)

res                      <- spia(de=DE.SPIA.vector,
                               all=names(mean.SPIA.vector),
                               organism = "hsa",
                               nB=10,
                               plots = FALSE,
                               combine = "fisher",
                               data.dir = kml.spia.dir)


# graphite ----------------------------------------------------------------
library(graphite)
humanReactome <- pathways("hsapiens", "reactome")
humanbiocarta <- pathways("hsapiens", "humancyc")
biocarta
names(humanReactome)[1:10]
pathwayDatabases()
prepareSPIA(humanReactome[1:1639], "reactomeTest")
prepareSPIA(humanbiocarta[1:2], "biocartaEx")
res <- runSPIA(de=DE.SPIA.vector, all=names(mean.SPIA.vector), "reactomeTest")
res <- runSPIA(de=DE.SPIA.vector, all=names(mean.SPIA.vector), "biocartaEx")

p <-  humanReactome[["Glucose metabolism"]]
p@title
head(nodes(p))
head(edges(p))
pSymbol <- convertIdentifiers(p, "SYMBOL")
head(nodes(pSymbol))


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

library(org.Hs.eg.db)
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

# let the probe simultaneous have up/down type to tag as mixed expression
mixed.updown.type <- intersect(node.up.Symbol_Entrezid$ENTREZID, node.down.Symbol_Entrezid$ENTREZID)

node.up.Symbol_Entrezid$EXPRESSION_TYPE[node.up.Symbol_Entrezid$ENTREZID %in% mixed.updown.type] <- "MIXED"
node.down.Symbol_Entrezid$EXPRESSION_TYPE[node.down.Symbol_Entrezid$ENTREZID %in% mixed.updown.type] <- "MIXED"

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
      one.result <-   interactions(psicquic, id=item, species="9606", provider=c("BioGrid", "IntAct", "Reactome-FIs", "UniProt"))
      total.result <- bind_rows(total.result,one.result[,c(1,2,7,9,12,13,14)])
      }
    return(total.result)
}

TF.interaction.protein <- SearchPSICQUIC(TF.Symbol_Entrezid$SYMBOL)
save(TF.interaction.protein, file = "TF_interaction_protein.RData")



CleanInteraction <- function(TF.interaction.protein){

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
id.symbol.from.uniprot <- AnnotationDbi::select(org.Hs.eg.db,
                            keys = id.uniprot.TF,
                            columns = c("SYMBOL"),
                            keytype = "UNIPROT")

id.symbol.from.entrez  <- AnnotationDbi::select(org.Hs.eg.db,
                                                 keys = id.entrez.TF,
                                                 columns = c("SYMBOL"),
                                                 keytype = "ENTREZID")



id.symbol.from.uniprot <- id.symbol.from.uniprot %>% filter(!is.na(SYMBOL))
#remove prefix
uniprot.TF.interaction.protein <- uniprot.TF.interaction.protein %>% mutate(A=str_replace_all(A,"uniprotkb:",""),B=str_replace_all(B,"uniprotkb:","")) 
entrez.TF.interaction.protein <- entrez.TF.interaction.protein %>% mutate(A=str_replace_all(A,"entrez gene\\/locuslink:",""),B=str_replace_all(B,"entrez gene\\/locuslink:","")) 


#filter without annotation id
uniprot.symbol <- vector()

for ( i in 1:nrow(id.symbol.from.uniprot)){
    uniprot.symbol[id.symbol.from.uniprot$UNIPROT[i]] <- id.symbol.from.uniprot$SYMBOL[i]    
}

uniprot.TF.interaction.protein <- uniprot.TF.interaction.protein %>% mutate(A=unname(uniprot.symbol[A]), B=unname(uniprot.symbol[B])) 
uniprot.TF.interaction.protein <- uniprot.TF.interaction.protein %>% filter(!is.na(A),!is.na(B))

# annotation the entrez with symbol
entrez.symbol <- vector()
for ( i in 1:nrow(id.symbol.from.entrez)){
    entrez.symbol[id.symbol.from.entrez$ENTREZID[i]] <- id.symbol.from.entrez$SYMBOL[i]    
}

entrez.TF.interaction.protein <- entrez.TF.interaction.protein %>% mutate(A=unname(entrez.symbol[A]), B=unname(entrez.symbol[B])) 
entrez.TF.interaction.protein <- entrez.TF.interaction.protein %>% filter(!is.na(A),!is.na(B))



total.clean.TF.interaction.protein <- bind_rows(entrez.TF.interaction.protein, uniprot.TF.interaction.protein)
return(total.clean.TF.interaction.protein)
}

total.clean.TF.interaction.protein <- CleanInteraction(TF.interaction.protein = TF.interaction.protein)






# gather the total protein list
total.005.entrez.updown <- c(node.up.term$ENTREZ, node.down.term$ENTREZ) %>% unique
total.005.symbol.from.entrez  <- AnnotationDbi::select(org.Hs.eg.db,
                                                keys = total.005.entrez.updown,
                                                columns = c("SYMBOL"),
                                                keytype = "ENTREZID")





# edge information
tmp <- total.clean.TF.interaction.protein %>% filter(A %in% total.005.symbol.from.entrez$SYMBOL & B %in% total.005.symbol.from.entrez$SYMBOL)




# let TF always on the left
tmp.1 <- tmp %>% filter(A %in% TF.Symbol_Entrezid$SYMBOL)
tmp.2 <- tmp %>% filter(B %in% TF.Symbol_Entrezid$SYMBOL & !(A %in% TF.Symbol_Entrezid$SYMBOL))
tmp.2 <- tmp.2 %>% dplyr::rename(A=B,B=A)
tmp.3 <- tmp %>% filter(A %in% TF.Symbol_Entrezid$SYMBOL & B %in% TF.Symbol_Entrezid$SYMBOL)

total.interaction.edge.inform <- bind_rows(tmp.1,tmp.2,tmp.3)
total.interaction.edge.inform <- tmp %>% unique
save(total.interaction.edge.inform, file="total_interaction_edge_information.RData")


# Secondary interaction retrieve ------------------------------------------
tmp.node <- total.interaction.edge.inform %>% dplyr::select(A,B) %>% unlist %>% unique()
tmp.node.removeTF <- tmp.node[!(tmp.node %in% TF.Symbol_Entrezid$SYMBOL)]
save(tmp.node.removeTF, file = "tmp_node_removeTF.RData")
load("~/Documents/Repository/code_in_lab/Total_Clean_SecondaryTF_interaction_protein.RData")

tmp <- total.clean.secondaryTF.interaction.protein %>% dplyr::select(A,B) %>% unique

# only gene exist in the o.o5 p expression 

tmp <- tmp %>% filter(A %in% total.005.symbol.from.entrez$SYMBOL & B %in% total.005.symbol.from.entrez$SYMBOL)
tmp <- tmp %>% arrange(A,B)
#MI definiation
#http://irefindex.org/wiki/index.php?title=Mapping_of_terms_to_MI_term_ids_-_iRefIndex_6.0
inter.gene.list <- foreach(i = 1:nrow(tmp)) %dopar% {
    tmp[i,1:2] %>% unlist
}

inter.gene.list <- inter.gene.list %>% map(unique) %>% keep(~length(.x) >1) 

start.i <- c()
end.i   <- c()
#26032
for (i in 1:260){
    start.i[i] <- 1+100*(i-1)
    end.i[i]   <- 100*i
}
start.i[261] <- 26001
end.i[261] <- 26032

# norm.time.list <- list()
# for (i in 1:2){
#     tmp.gene.list <- inter.gene.list[start.i[i]:end.i[i]] %>% unlist %>% unique()
#     norm.time.list[[i]] <- norm.time %>% filter(name %in% tmp.gene.list)
# }

tmp.gene.list <- foreach(i = 1:261) %dopar%{
    inter.gene.list[start.i[i]:end.i[i]] %>% unlist %>% unique()
}

norm.time.list <- foreach(i = 1:261) %dopar%{
    norm.time %>% filter(name %in% tmp.gene.list[[i]])
}

total.secondary.interaction_correlation.data.frame <- data.frame()
for (i in 1:261){
    tmp.norm.time <- norm.time.list[[i]]
    for ( j in start.i[i]:end.i[i]){
        norm.time.data  <- tmp.norm.time %>% filter(name %in% inter.gene.list[[j]])
        tmp.correlation <- cor(unlist(norm.time.data[1,3:6]), unlist(norm.time.data[2,3:6]))
        tmp.interaction_correlation.data.frame <- inter.gene.list[[j]] %>% t %>% data.frame %>% mutate(correlation=tmp.correlation)
        total.secondary.interaction_correlation.data.frame <- bind_rows(total.secondary.interaction_correlation.data.frame, tmp.interaction_correlation.data.frame)
    }
}
save(total.secondary.interaction_correlation.data.frame, file="total_secondary_interaction_correlation_data.RData")

# export the interaction into cytoscape
total.secondary.interaction_correlation.data.frame <- total.secondary.interaction_correlation.data.frame %>% filter(!is.na(correlation))
total.secondary.interaction_correlation.data.frame  <- total.secondary.interaction_correlation.data.frame  %>% rename(A=X1,B=X2)
tmp.total.secondary.interaction_correlation.data.frame <- total.secondary.interaction_correlation.data.frame %>% filter(correlation > 0.5 | correlation < -0.5)
write_delim(tmp.total.secondary.interaction_correlation.data.frame, path="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/secondary TF interaction/interaction_edge.txt", delim = "\t")



# Add the correlation relationship ----------------------------------------
load("~/Documents/Repository/code_in_lab/norm_qn_log_time.Rdata")
load("~/Documents/Repository/code_in_lab/hgu133plus2_annotation_v1.RData")
load("~/Documents/Repository/code_in_lab/TF_interaction_protein.RData")
annotation.use <- hgu133plus2.probe.annotate[,c(5,8,9)]
rowname.probe <- rownames(norm.time.series )
norm.time <- norm.time.series
colnames(norm.time) <- c("p14","p3","p6-1","p6-2","p6-3","p8")
norm.time <- norm.time %>% as.data.frame() %>% mutate(Probe=rowname.probe) %>% left_join(.,annotation.use, by="Probe")
norm.time <- norm.time %>% mutate(p6=(`p6-1`+`p6-2`+`p6-3`)/3) %>% dplyr::select(-`p6-1`, -`p6-2`, -`p6-3`)
norm.time <- norm.time %>% dplyr::group_by(name, ENTREZID) %>% dplyr::summarise(p3=median(p3), p6=median(p6), p8=median(p8), p14=median(p14))

# from the ppi interaction
library(doParallel)
library(foreach)

# for ( i in 1:nrow(total.interaction.edge.inform)){
# 
# inter.gene <- total.interaction.edge.inform[i,1:2] %>% unlist
# cor.data <- norm.time %>% filter(name %in% inter.gene)
# total.interaction.edge.correlation.vector[i] <- cor(unlist(cor.data[1,3:6]), unlist(cor.data[2,3:6]))
# 
# 
# }

registerDoParallel(cores=4)
inter.gene.list <- foreach(i = 1:nrow(total.interaction.edge.inform)) %dopar% {
    total.interaction.edge.inform[i,1:2] %>% unlist
}

inter.gene.list <- inter.gene.list %>% map(unique) %>% keep(~length(.x) >1) 

start.i <- c()
end.i   <- c()
#7477
for (i in 1:74){
    start.i[i] <- 1+100*(i-1)
    end.i[i]   <- 100*i
}
start.i[75] <- 7401
end.i[75] <- 7477

# norm.time.list <- list()
# for (i in 1:2){
#     tmp.gene.list <- inter.gene.list[start.i[i]:end.i[i]] %>% unlist %>% unique()
#     norm.time.list[[i]] <- norm.time %>% filter(name %in% tmp.gene.list)
# }

tmp.gene.list <- foreach(i = 1:75) %dopar%{
    inter.gene.list[start.i[i]:end.i[i]] %>% unlist %>% unique()
}

norm.time.list <- foreach(i = 1:75) %dopar%{
    norm.time %>% filter(name %in% tmp.gene.list[[i]])
}

total.interaction_correlation.data.frame <- data.frame()
for (i in 1:75){
    tmp.norm.time <- norm.time.list[[i]]
    for ( j in start.i[i]:end.i[i]){
       norm.time.data  <- tmp.norm.time %>% filter(name %in% inter.gene.list[[j]])
       tmp.correlation <- cor(unlist(norm.time.data[1,3:6]), unlist(norm.time.data[2,3:6]))
       tmp.interaction_correlation.data.frame <- inter.gene.list[[j]] %>% t %>% data.frame %>% mutate(correlation=tmp.correlation)
       total.interaction_correlation.data.frame <- bind_rows(total.interaction_correlation.data.frame, tmp.interaction_correlation.data.frame)
    }
}
save(total.interaction_correlation.data.frame, file="total_interaction_correlation_data.RData")


final.interaction_correlation.data.frame <- total.interaction_correlation.data.frame %>% filter(!is.na(correlation)) %>% unique
final.interaction_correlation.data.frame <- final.interaction_correlation.data.frame %>% dplyr::rename(A=X1,B=X2)



write_delim(final.interaction_correlation.data.frame, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/second/interaction_edge.txt", delim = "\t")









tmp.interaction <- total.interaction.edge.inform
tmp.interaction <- tmp.interaction[,c(1,2,5)]
tmp.interaction <- unique(tmp.interaction)

#tmp.interaction <- tmp.interaction %>% filter(str_detect(type,"direct interaction") | str_detect(type,"physical interaction") |str_detect(type, "association"))
tmp.interaction <- tmp.interaction %>% filter(A %in% total.node$SYMBOL) %>% filter(B %in% total.node$SYMBOL)





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




# node information
total.node.tag.p.value <- total.node.tag.p.value %>% mutate(TF="") 
total.node.tag.p.value$TF[total.node.tag.p.value$ENTREZID %in% TF.Symbol_Entrezid$ENTREZID ] <- "TF"
total.node.tag.p.value$TF[!(total.node.tag.p.value$ENTREZID %in% TF.Symbol_Entrezid$ENTREZID )] <- "Not TF"


filter.final.interaction_correlation.data.frame <- final.interaction_correlation.data.frame %>% filter(correlation > 0.5 | correlation < -0.5)
tmp.total.symbol <- filter.final.interaction_correlation.data.frame[,c(1,2)] %>% unlist %>% unique 

total.node.tag.p.value.in.interaction <- total.node.tag.p.value %>% filter(SYMBOL %in% tmp.total.symbol)

write_delim(total.node.tag.p.value.in.interaction, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/third/interaction_node.txt", delim = "\t")
write_delim(filter.final.interaction_correlation.data.frame,  path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/third/interaction_edge.txt", delim = "\t")














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


# ECM and RNA splicing and Metabolism
ECM.RNA.Metabolism.related.term <- c(ECM.GO.related.term, Metabolism.GO.related.term, "GO:0043039", "GO:0006418", "GO:0043603", "GO:0006518", "GO:0043604", "GO:0043043", "GO:0006412", "GO:0008380", "GO:0006396", "GO:0000375", "GO:0071704", "GO:0000377", "GO:0000398", "GO:0006397", "GO:0043484")
total.entrezid.up.related.ECM.RNA.Metabolism   <- total.entrezid.up.GO %>% filter(GO %in% ECM.GO.related.term) %>% with(ENTREZID)
total.entrezid.down.related.ECM.RNA.Metabolism <- total.entrezid.down.GO %>% filter(GO %in% ECM.GO.related.term) %>% with(ENTREZID)

SYMBOL.related.ECM.RNA.Metabolism   <- total.node.tag.p.value %>% filter(ENTREZID %in% c(total.entrezid.up.related.ECM.RNA.Metabolism , total.entrezid.down.related.ECM.RNA.Metabolism )) %>% with(SYMBOL)
sub.edge.related.ECM.RNA.Metabolism <- final.interaction_correlation.data.frame %>% filter(A %in% SYMBOL.related.ECM.RNA.Metabolism | B %in% SYMBOL.related.ECM.RNA.Metabolism)
sub.node.related.ECM.RNA.Metabolism <- total.node.tag.p.value %>% filter(SYMBOL %in% unique(c(sub.edge.related.ECM.RNA.Metabolism$A, sub.edge.related.ECM.RNA.Metabolism$B)))

write_delim(sub.node.related.ECM.RNA.Metabolism, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/ECM_RNA_Metabolism/interaction_node.txt", delim = "\t")
write_delim(sub.edge.related.ECM.RNA.Metabolism, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/ECM_RNA_Metabolism/interaction_edge.txt", delim = "\t")


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
Metabolism.GO.related.term <- c("GO:0008152","GO:0044237","GO:0071704",'GO:0044248',"GO:0019538","GO:0009056","GO:0043604", "GO:0043603", "GO:0043603", "GO:0043043", "GO:0006518", "GO:0031396", "GO:0016567", "GO:1900040", "GO:0043039", "GO:1903322" ,"GO:0045184", "GO:0051248")
total.entrezid.up.related.Metabolism   <- total.entrezid.up.GO %>% filter(GO %in% Metabolism.GO.related.term) %>% with(ENTREZID)
total.entrezid.down.related.Metabolism <- total.entrezid.down.GO %>% filter(GO %in% Metabolism.GO.related.term) %>% with(ENTREZID)

SYMBOL.related.Metabolism <- total.node.tag.p.value %>% filter(ENTREZID %in% c(total.entrezid.up.related.Metabolism, total.entrezid.down.related.Metabolism)) %>% with(SYMBOL)
sub.edge.related.Metabolism <- tmp.interaction %>% filter(A %in% SYMBOL.related.Metabolism | B %in% SYMBOL.related.Metabolism)
sub.node.related.Metabolism <- total.node.tag.p.value %>% filter(SYMBOL %in% unique(c(sub.edge.related.Metabolism$A, sub.edge.related.Metabolism$B)))

write_delim(sub.node.related.Metabolism, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/Metabolism_related_subnetwork/interaction_node.txt", delim = "\t")
write_delim(sub.edge.related.Metabolism, path ="/Users/Weitinglin/Documents/2016 實驗室資料處理/Cytoscape/Metabolism_related_subnetwork/interaction_edge.txt", delim = "\t")
# Metabolism related


# ****************************************************** ------------------
#load("~/Documents/Repository/code_in_lab/set9_t_bh_final.Rdata")
load("~/Documents/Repository/code_in_lab/set10_aov_tukey_final.Rdata")
load("~/Documents/Repository/code_in_lab/hgu133plus2_annotation_v1.RData")

hgu133plus2.annotate.set <- hgu133plus2.probe.annotate %>% dplyr::select(name, Probe, ENTREZID)
# remove the AFFX 62 probes


# Protein Protein interaction ---------------------------------------------
load("~/Documents/Repository/code_in_lab/set10_005_total_symbol_interaction.Rdata")
load("~/Documents/Repository/code_in_lab/set10_005_cortest.Rdata")
# Set 10 ------------------------------------------------------------------

# function
calculateTopGO <- function(up.set10.Sphere_CLS.all.GO, up.set10.Sphere_CLS.005, number){
    
    up.set10.Sphere_CLS.all.GO.list <- up.set10.Sphere_CLS.all.GO %>% filter(EVIDENCE != "IEA") %>% split(.$ONTOLOGY)
    up.set10.Sphere_CLS.all.GO.BP   <- up.set10.Sphere_CLS.all.GO.list$BP %>% split(.$ENTREZID) %>% map(with, GO)
    up.set10.Sphere_CLS.all.GO.BP.GO2geneID <- inverseList(up.set10.Sphere_CLS.all.GO.BP)
    
    
    set10.geneNames <- names(up.set10.Sphere_CLS.all.GO.BP)
    set10.geneList <- factor(as.integer(set10.geneNames %in% as.character(up.set10.Sphere_CLS.005$ENTREZID[!is.na(up.set10.Sphere_CLS.005$ENTREZID)])))
    names(set10.geneList) <- set10.geneNames
    
    GOdata <- new("topGOdata", ontology = "BP", allGenes = set10.geneList,
                  annot = annFUN.gene2GO, gene2GO = up.set10.Sphere_CLS.all.GO.BP)
    
    # calculate the result
    resultFisher <- runTest(GOdata , algorithm = "classic", statistic = "fisher")
    resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    
    allRes  <- GenTable(GOdata, classicFisher= resultFisher,
                        classicKS = resultKS, elimKS = resultKS.elim,
                        orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = number)
    
    allRes  <- allRes %>% transmute(GO.ID=GO.ID,
                                    Term=Term,
                                    Annotated=Annotated,
                                    Significant=Significant,
                                    Expected=Expected,
                                    classicFisher=as.numeric(classicFisher),
                                    classicKS=as.numeric(classicKS),
                                    elimKS=as.numeric(elimKS))
    return(allRes)
}


# _sepDesign --------------------------------------------------------------

set10.Sphere_CLS        <- set10.test.tukey.fine.result %>%
    filter(case == "Sphere-CLS") %>%
    filter(!str_detect(Probe, "AFFX")) %>%
    mutate(BH.p.adj=p.adjust(p.adj, method="BH"))

set10.Sphere_CLS <- left_join(set10.Sphere_CLS, hgu133plus2.annotate.set, by="Probe")

set10.Timeseries_Sphere <- set10.test.tukey.fine.result %>%
    filter(case == "Timeseries-Sphere") %>%
    filter(!str_detect(Probe, "AFFX")) %>%
    mutate(BH.p.adj=p.adjust(p.adj, method="BH"))

set10.Timeseries_Sphere <- left_join(set10.Timeseries_Sphere, hgu133plus2.annotate.set, by="Probe")

set10.Timeseries_CLS <- set10.test.tukey.fine.result %>%
    filter(case == "Timeseries-CLS") %>%
    filter(!str_detect(Probe, "AFFX")) %>%
    mutate(BH.p.adj=p.adjust(p.adj, method="BH"))

set10.Timeseries_CLS <- left_join(set10.Timeseries_CLS, hgu133plus2.annotate.set, by="Probe")


# _sepUpandDown ------------------------------------------------------------



up.set10.Sphere_CLS     <- set10.Sphere_CLS %>% filter(diff > 0)
down.set10.Sphere_CLS   <- set10.Sphere_CLS %>% filter(diff < 0)

up.set10.Timeseries_Sphere <- set10.Timeseries_Sphere %>% filter(diff < 0)
down.set10.Timeseries_Sphere <- set10.Timeseries_Sphere %>% filter(diff > 0)



# _filterPvalue ------------------------------------------------------------

up.set10.Sphere_CLS.005              <-     up.set10.Sphere_CLS %>% filter(BH.p.adj < 0.05)
down.set10.Sphere_CLS.005            <-   down.set10.Sphere_CLS %>% filter(BH.p.adj < 0.05)

up.set10.Timeseries_Sphere.005       <-   up.set10.Timeseries_Sphere %>% filter(BH.p.adj < 0.05)
down.set10.Timeseries_Sphere.005     <-   down.set10.Timeseries_Sphere %>% filter(BH.p.adj < 0.05)

# _gene list ---------------------------------------------------------------

up.set10.intersect    <- intersect(up.set10.Sphere_CLS.005$ENTREZID, up.set10.Timeseries_Sphere.005$ENTREZID )
down.set10.intersect  <- intersect(down.set10.Sphere_CLS.005$ENTREZID, down.set10.Timeseries_Sphere.005$ENTREZID)


total.set10.Sphere_CLS.005 <- c(up.set10.Sphere_CLS.005$name,down.set10.Sphere_CLS.005$name)
total.set10.Sphere_CLS.005 <- c(up.set10.Sphere_CLS.005$name)
total.set10.Sphere_CLS.005 <- c(down.set10.Sphere_CLS.005$name)

total.set10.Sphere_CLS.005 <- total.set10.Sphere_CLS.005[!is.na(total.set10.Sphere_CLS.005)]
# GO analysis -------------------------------------------------------------
library(topGO)
library(GO.db)
library(org.Hs.eg.db)


# set10 -------------------------------------------------------------------


# _Prepare the data structure  topGO need

up.set10.Sphere_CLS.all.entrezID   <-  up.set10.Sphere_CLS   %>%  with(ENTREZID)

up.set10.Sphere_CLS.all.GO <- AnnotationDbi::select(org.Hs.eg.db,
                                                    keys = as.character(up.set10.Sphere_CLS.all.entrezID[!is.na(up.set10.Sphere_CLS.all.entrezID)]),
                                                    columns = "GO",
                                                    keytype = "ENTREZID")


down.set10.Sphere_CLS.all.entrezID   <-  down.set10.Sphere_CLS  %>%  with(ENTREZID)

down.set10.Sphere_CLS.all.GO <- AnnotationDbi::select(org.Hs.eg.db,
                                                      keys = as.character(down.set10.Sphere_CLS.all.entrezID[!is.na(down.set10.Sphere_CLS.all.entrezID)]),
                                                      columns = "GO",
                                                      keytype = "ENTREZID")

up.set10.Timeseries_Sphere.all.entrezID   <-   up.set10.Timeseries_Sphere %>% with(ENTREZID) 

down.set10.Timeseries_Sphere.all.entrezID <-    down.set10.Timeseries_Sphere %>%  with(ENTREZID)


up.set10.Timeseries_Sphere.all.GO <- AnnotationDbi::select(org.Hs.eg.db,
                                                                 keys = as.character(up.set10.Timeseries_Sphere.all.entrezID[!is.na(up.set10.Timeseries_Sphere.all.entrezID)]),
                                                                 columns = "GO",
                                                                 keytype = "ENTREZID")
down.set10.Timeseries_Sphere.all.GO <- AnnotationDbi::select(org.Hs.eg.db,
                                                             keys = as.character(down.set10.Timeseries_Sphere.all.entrezID[!is.na(down.set10.Timeseries_Sphere.all.entrezID)]),
                                                             columns = "GO",
                                                             keytype = "ENTREZID")

# _Turn data.frame into list-like data structure ---------------------------

up.set10.Sphere_CLS.GOanalysis.result   <- calculateTopGO(up.set10.Sphere_CLS.all.GO, up.set10.Sphere_CLS.005, 200)
down.set10.Sphere_CLS.GOanalysis.result <- calculateTopGO(down.set10.Sphere_CLS.all.GO, up.set10.Sphere_CLS.005, 200)

up.set10.Timeseries_Sphere.GOanalysis.result <- calculateTopGO(up.set10.Timeseries_Sphere.all.GO, up.set10.Timeseries_Sphere.005, 200)
down.set10.Timeseries_Sphere.GOanalysis.result <- calculateTopGO(down.set10.Timeseries_Sphere.all.GO, down.set10.Timeseries_Sphere.005, 200)

write.csv(up.set10.Sphere_CLS.GOanalysis.result, file = "up_set10_Sphere_CLS_GO_005.csv")
write.csv(down.set10.Sphere_CLS.GOanalysis.result, file = "down_set10_Sphere_CLS_GO_005.csv")
write.csv(up.set10.Timeseries_Sphere.GOanalysis.result, file = "up_set10_Timeseries_Sphere_GO_005.csv")
write.csv(down.set10.Timeseries_Sphere.GOanalysis.result, file = "down_set10_Timeseries_Sphere_GO_005.csv")


index.set10.up    <- intersect(up.set10.Sphere_CLS.GOanalysis.result$GO.ID, up.set10.Timeseries_Sphere.GOanalysis.result$GO.ID)
index.set10.down  <- intersect(down.set10.Sphere_CLS.GOanalysis.result$GO.ID, down.set10.Timeseries_Sphere.GOanalysis.result$GO.ID)

up.set10.intersect.GO.term   <- up.set10.Sphere_CLS.GOanalysis.result %>% filter(GO.ID %in% index.set10.up)
down.set10.intersect.GO.term <- down.set10.Sphere_CLS.GOanalysis.result %>%  filter(GO.ID %in% index.set10.down)


total.intersect.index <- intersect(down.set10.intersect.GO.term$GO.ID, up.set10.intersect.GO.term$GO.ID)
total.intersect.GO.term <- down.set10.intersect.GO.term %>% filter(GO.ID %in% total.intersect.index)

filter.up.set10.intersect.GO.term <- up.set10.intersect.GO.term %>%  filter(!(GO.ID %in% total.intersect.index))
filter.down.set10.intersect.GO.term <- down.set10.intersect.GO.term %>%  filter(!(GO.ID %in% total.intersect.index))


# summarise the GO analysis result ----------------------------------------

up_set10_Sphere_CLS_005_topGO_result <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/up_set10_Sphere_CLS_GO_005.csv", delim=",")
down_set10_Sphere_CLS_005_topGO_result <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/down_set10_Sphere_CLS_GO_005.csv", delim=",")

up_set10_Timeseries_Sphere_005_topGO_result <- read_delim(file = "/Users/Weitinglin/Documents/Repository/code_in_lab/up_set10_Timeseries_Sphere_GO_005.csv", delim=",")
down_set10_Timeseries_Sphere_005_topGO_result  <-  read_delim(file= "/Users/Weitinglin/Documents/Repository/code_in_lab/down_set10_Timeseries_Sphere_GO_005.csv", delim=",")


up.set10.TopGO.intersect.GO.ID    <- intersect(up_set10_Sphere_CLS_005_topGO_result$GO.ID, up_set10_Timeseries_Sphere_005_topGO_result$GO.ID)

down.set10.TopGO.intersect.GO.ID  <- intersect(down_set10_Sphere_CLS_005_topGO_result$GO.ID, down_set10_Timeseries_Sphere_005_topGO_result$GO.ID)


up_set10_Sphere_CLS_005_topGO_result.filter        <- up_set10_Sphere_CLS_005_topGO_result %>% filter(GO.ID %in% up.set10.TopGO.intersect.GO.ID) %>% arrange(classicFisher, desc(Annotated)) %>% select(-X1)
up_set10_Timeseries_Sphere_005_topGO_result.filter <- up_set10_Timeseries_Sphere_005_topGO_result %>% filter(GO.ID %in% up.set10.TopGO.intersect.GO.ID) %>%  arrange(classicFisher, desc(Annotated)) %>%  select(-X1)

down_set10_Sphere_CLS_005_topGO_result.filter        <- down_set10_Sphere_CLS_005_topGO_result %>% filter(GO.ID %in% down.set10.TopGO.intersect.GO.ID) %>% arrange(classicFisher, desc(Annotated)) %>% select(-X1)
down_set10_Timeseries_Sphere_005_topGO_result.filter <- down_set10_Timeseries_Sphere_005_topGO_result %>% filter(GO.ID %in% down.set10.TopGO.intersect.GO.ID) %>% arrange(classicFisher, desc(Annotated)) %>% select(-X1)


# Use gprofiler -----------------------------------------------------------



up.set10.Sphere_CLS.001             <-     up.set10.Sphere_CLS %>% filter(BH.p.adj < 0.01)
down.set10.Sphere_CLS.001            <-   down.set10.Sphere_CLS %>% filter(BH.p.adj < 0.01)

up.set10.Timeseries_Sphere.001       <-   up.set10.Timeseries_Sphere %>% filter(BH.p.adj < 0.01)
down.set10.Timeseries_Sphere.001     <-   down.set10.Timeseries_Sphere %>% filter(BH.p.adj < 0.001)
library(gProfileR)

up.set10.Sphere_CLS.probe <- up.set10.Sphere_CLS.001$Probe
up.set10.Timeseries_Sphere.probe <- up.set10.Timeseries_Sphere.001$Probe 
down.set10.Sphere_CLS.probe <- down.set10.Sphere_CLS.001$Probe
down.set10.Timeseries_Sphere.probe<- down.set10.Timeseries_Sphere.0001$Probe



hgu133plus2.probe <- as.character(annotated.entrez.symbol$Probe)

up.set10.Sphere_CLS_001_gprofiler.result <- gprofiler(  query = up.set10.Sphere_CLS.probe,
                                                        organism = "hsapiens",
                                                        ordered_query = T,
                                                        max_set_size = 1000,
                                                        min_isect_size = 0,
                                                        significant = T,
                                                        underrep = T,
                                                        exclude_iea = T,
                                                        correction_method = "fdr",
                                                        hier_filtering = "moderate",
                                                        custom_bg = hgu133plus2.probe,
                                                        src_filter = c("GO:BP","GO:CC"))




down.set10.Sphere_CLS_001_gprofiler.result <- gprofiler(  query = down.set10.Sphere_CLS.probe,
                                                          organism = "hsapiens",
                                                          ordered_query = T,
                                                          max_set_size = 1000,
                                                          min_isect_size = 0,
                                                          significant = T,
                                                          underrep = T,
                                                          exclude_iea = T,
                                                          correction_method = "fdr",
                                                          hier_filtering = "moderate",
                                                          custom_bg = hgu133plus2.probe,
                                                          src_filter = c("GO:BP","GO:CC"))





up.set10.Timeseries_Sphere_001_gprofiler.result <- gprofiler(  query = up.set10.Timeseries_Sphere.probe,
                                                                    organism = "hsapiens",
                                                                    ordered_query = T,
                                                                    max_set_size = 1000,
                                                                    min_isect_size = 0,
                                                                    significant = T,
                                                                    underrep = T,
                                                                    exclude_iea = T,
                                                                    correction_method = "fdr",
                                                                    hier_filtering = "moderate",
                                                                    custom_bg = hgu133plus2.probe,
                                                                    src_filter = c("GO:BP","GO:CC"))

down.set10.STimeseries_Sphere_001_gprofiler.result <- gprofiler(  query = down.set10.Timeseries_Sphere.001$Probe ,
                                                                  organism = "hsapiens",
                                                                  ordered_query = T,
                                                                  max_set_size = 1000,
                                                                  min_isect_size = 0,
                                                                  significant = T,
                                                                  underrep = T,
                                                                  exclude_iea = T,
                                                                  correction_method = "fdr",
                                                                  hier_filtering = "moderate",
                                                                  custom_bg = hgu133plus2.probe,
                                                                  src_filter = c("GO:BP","GO:CC"))


# SPIA analysis -----------------------------------------------------------
library(affy)
library(SPIA)
load("~/Documents/Repository/code_in_lab/Set9_10_11_Exprs.Rdata")
load("~/Documents/Repository/code_in_lab/set10_adjust_sep.Rdata")
kml.spia.dir <- "/Users/Weitinglin/Documents/R_scripts/Lab/kegg/spiaXML/"
# 直接變換在spiaXML下的hsaSPIA，就可以在KEGG, BioGrid, Reactome間轉換
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

hash.table <- vector()
for ( i in 1:length(set10.Sphere_CLS$Probe)){
    hash.table[[set10.Sphere_CLS$Probe[i]]] <- set10.Sphere_CLS$ENTREZID[i]
}

# set 9
expression <- exprs(Set9.Exprs.data)
set9.norm <- quantile_normalization(expression, method="median")

# set 10
expression <- exprs(Set10.Exprs.data)
set10.norm <- quantile_normalization(expression, method = "median")



#Sphere_CLS
cutoff.p <- 0.05
Sphere_CLS.subset <- set10.norm[,c(4,5,6,1,2,3)]
Sphere_CLS.subset.log2change <- calculateLog2foldchange(Sphere_CLS.subset)
set10.Sphere_CLS.sig <- set10.Sphere_CLS %>% filter(BH.p.adj < cutoff.p) %>% with(Probe)

SPIA.vector <- Sphere_CLS.subset.log2change$logFC
names(SPIA.vector) <- hash.table[Sphere_CLS.subset.log2change$Probe]
SPIA.vector <- SPIA.vector[!is.na(names(SPIA.vector))]
mean.SPIA.vector <- meanDuplicateEntrezID(SPIA.vector)
set10.Sphere_CLS.sig.id <- hash.table[set10.Sphere_CLS.sig] %>% unname %>% unique
set10.Sphere_CLS.sig.id <- set10.Sphere_CLS.sig.id[!is.na(set10.Sphere_CLS.sig.id)]
DE.SPIA.vector <- mean.SPIA.vector[names(mean.SPIA.vector) %in% (set10.Sphere_CLS.sig.id)]
set10.Sphere_CLS.spia.result.nb.2000 <- spia(de=DE.SPIA.vector,
                                         all=names(mean.SPIA.vector),
                                        organism = "hsa",
                                        nB = 2000,
                                        plots = FALSE,
                                        combine = "fisher",
                                        data.dir = kml.spia.dir)
write.csv(set10.Sphere_CLS.spia.result.nb.2000, file="set10_005_spia_SphereCLS_2000_kegg.csv")



#Sphere_Timeseries
Sphere_Timeseries.subset <- set10.norm[,c(4,5,6,7,8,9)]
Sphere_Timeseries.subset.log2change <- calculateLog2foldchange(Sphere_Timeseries.subset) 
set10.Sphere_Timeseries.sig <- set10.Timeseries_Sphere %>% filter(BH.p.adj < cutoff.p) %>% with(Probe)


SPIA.vector <- Sphere_Timeseries.subset.log2change$logFC
names(SPIA.vector) <- hash.table[Sphere_Timeseries.subset.log2change$Probe]
SPIA.vector <- SPIA.vector[!is.na(names(SPIA.vector))]
mean.SPIA.vector <- meanDuplicateEntrezID(SPIA.vector)
set10.Sphere_Timeseries.sig.id <- hash.table[set10.Sphere_Timeseries.sig] %>% unname %>% unique
set10.Sphere_Timeseries.sig.id <- set10.Sphere_Timeseries.sig.id[!is.na(set10.Sphere_Timeseries.sig.id)]
DE.SPIA.vector <- mean.SPIA.vector[names(mean.SPIA.vector) %in% (set10.Sphere_Timeseries.sig.id)]
set10.Sphere_Timeseries.spia.nb2000.result <- spia(de=DE.SPIA.vector,
                                                all=names(mean.SPIA.vector),
                                                organism = "hsa",
                                                nB = 2000,
                                                plots = FALSE,
                                                combine = "fisher",
                                                data.dir = kml.spia.dir)
write.csv(set10.Sphere_Timeseries.spia.nb2000.result, file="set10_005_spia_SphereTimeseries_2000_kegg.csv")



# summarise the Signal pathway analysis with SPIA -------------------------


set10.Sphere_CLS.spia.result.kegg           <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereCLS_kegg.csv", delim = ",")           %>% filter(pGFdr < 0.05)  %>%  mutate(name_dir=paste0(Name, "-", Status))
set10.Sphere_Timeseries.spia.result.kegg    <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereTimeseries_KEGG.csv",delim=",")       %>% filter(pGFdr < 0.05)   %>%  mutate(name_dir=paste0(Name, "-", Status))   
 
set10.Sphere_CLS.spia.result.Reactom        <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereCLS_reactome.csv", delim=",")         %>% filter(pGFdr < 0.05)  %>%  mutate(name_dir=paste0(Name, "-", Status))    
set10.Sphere_Timeseries.spia.result.Reactom <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereTimeseries_Reactome.csv", delim=",")  %>% filter(pGFdr < 0.05)  %>%  mutate(name_dir=paste0(Name, "-", Status))          


# KEGG
set10.kegg.intersect.index <- intersect(set10.Sphere_CLS.spia.result.kegg$name_dir, set10.Sphere_Timeseries.spia.result.kegg$name_dir)
set10.Sphere_CLS.spia.result.kegg.filter <- set10.Sphere_CLS.spia.result.kegg %>%  filter(name_dir %in% set10.kegg.intersect.index)
set10.Sphere_Timeseries.spia.result.kegg.filter <- set10.Sphere_Timeseries.spia.result.kegg %>%  filter(name_dir %in% set10.kegg.intersect.index)

# Reactome

set10.reactome.intersect.index <- intersect(set10.Sphere_CLS.spia.result.Reactom$name_dir, set10.Sphere_Timeseries.spia.result.Reactom$name_dir)
set10.Sphere_CLS.spia.result.Reactom.filter <- set10.Sphere_CLS.spia.result.Reactom %>% filter(name_dir %in% set10.reactome.intersect.index)
set10.Sphere_Timeseries.spia.result.Reactom.filter <- set10.Sphere_Timeseries.spia.result.Reactom %>% filter(name_dir %in% set10.reactome.intersect.index)

TF.interaction.protein <- set10.005.total.symbol.interaction 
# 101578

# Network create ----------------------------------------------------------
CleanInteraction <- function(TF.interaction.protein){
    
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
    id.symbol.from.uniprot <- AnnotationDbi::select(org.Hs.eg.db,
                                                    keys = id.uniprot.TF,
                                                    columns = c("SYMBOL"),
                                                    keytype = "UNIPROT")
    
    id.symbol.from.entrez  <- AnnotationDbi::select(org.Hs.eg.db,
                                                    keys = id.entrez.TF,
                                                    columns = c("SYMBOL"),
                                                    keytype = "ENTREZID")
    
    
    id.symbol.from.uniprot <- id.symbol.from.uniprot %>% filter(!is.na(SYMBOL))
    id.symbol.from.entrez  <- id.symbol.from.entrez %>% filter(!is.na(SYMBOL)) 
    #remove prefix
    uniprot.TF.interaction.protein <- uniprot.TF.interaction.protein %>% mutate(A=str_replace_all(A,"uniprotkb:",""),B=str_replace_all(B,"uniprotkb:",""))
    entrez.TF.interaction.protein <- entrez.TF.interaction.protein %>% mutate(A=str_replace_all(A,"entrez gene\\/locuslink:",""),B=str_replace_all(B,"entrez gene\\/locuslink:",""))
    
    
    #filter without annotation id
    uniprot.symbol <- vector()
    
    for ( i in 1:nrow(id.symbol.from.uniprot)){
        uniprot.symbol[id.symbol.from.uniprot$UNIPROT[i]] <- id.symbol.from.uniprot$SYMBOL[i]
    }
    
    uniprot.TF.interaction.protein <- uniprot.TF.interaction.protein %>% mutate(A=unname(uniprot.symbol[A]), B=unname(uniprot.symbol[B]))
    uniprot.TF.interaction.protein <- uniprot.TF.interaction.protein %>% filter(!is.na(A),!is.na(B))
    
    # annotation the entrez with symbol
    entrez.symbol <- vector()
    for ( i in 1:nrow(id.symbol.from.entrez)){
        entrez.symbol[id.symbol.from.entrez$ENTREZID[i]] <- id.symbol.from.entrez$SYMBOL[i]
    }
    
    entrez.TF.interaction.protein <- entrez.TF.interaction.protein %>% mutate(A=unname(entrez.symbol[A]), B=unname(entrez.symbol[B]))
    entrez.TF.interaction.protein <- entrez.TF.interaction.protein %>% filter(!is.na(A),!is.na(B))
    
    
    
    total.clean.TF.interaction.protein <- bind_rows(entrez.TF.interaction.protein, uniprot.TF.interaction.protein)
    return(total.clean.TF.interaction.protein)
}

library(org.Hs.eg.db)

clean.set9.005.total.ppi <- CleanInteraction(set9.005.total.symbol.interaction)
clean.set10.005.total.ppi <- CleanInteraction(set10.005.total.symbol.interaction) %>%  unique



# prepare the output file for cytoscape


# Set 9 network ------------------------------------------------------------


# edge file
load("~/Documents/Repository/code_in_lab/set9_005_cortest.Rdata")

clean.set9.005.total.ppi.output <-  set9.005.interaction.cor.test %>% dplyr::select(A,B,estimate, p.value) %>% unique
# 89018
filter.clean.set9.005.total.ppi.output <- clean.set9.005.total.ppi.output %>%  filter(p.value < 0.05)
# 5665
write_delim(filter.clean.set9.005.total.ppi.output, path ="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set9/interaction_edge.txt", delim = "\t")


# node file
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

load("~/Documents/Repository/code_in_lab/set9_t_bh_final.Rdata")
load("~/Documents/Repository/code_in_lab/hgu133plus2_annotation_v1.RData")

hgu133plus2.annotate.set <- hgu133plus2.probe.annotate %>% dplyr::select(name, Probe, ENTREZID)
set9.t.bamjamin.final.result  <- set9.t.bamjamin.final.result %>% dplyr::rename(Probe=gene)
set9.t.bamjamin.final.result <- left_join(set9.t.bamjamin.final.result, hgu133plus2.annotate.set, by="Probe")

set9.final.ppi.symbol.list<- unique(c(filter.clean.set9.005.total.ppi.output$A, filter.clean.set9.005.total.ppi.output$B))

set9.node.annotation <- set9.t.bamjamin.final.result %>% filter(adjusted.p < 0.05) %>%  filter(name %in% set9.final.ppi.symbol.list)

# functional annotation: search from the panther database
set9.node.panther.term <- search_panther_database(set9.node.annotation$ENTREZID)
set9.node.panther.class.term <- set9.node.panther.term %>%  filter(TERM_TYPE == "CLASS_TERM")
set9.node.panther.TF.related <- set9.node.panther.class.term %>% filter(str_detect(PANTHER_TERM,"transcription"))
set9.node.panther.TF.related.entrezid <- set9.node.panther.TF.related$ENTREZ %>% unique()

set9.node.panther.class.term.merged <- set9.node.panther.class.term %>%  group_by(ENTREZ) %>% summarise(PATHER_TERM=paste(PANTHER_TERM, collapse=","))
set9.node.panther.class.term.merged <- set9.node.panther.class.term.merged %>% dplyr::rename(ENTREZID = ENTREZ)
# add annotation on
set9.node.annotation$TFrelated[set9.node.annotation$ENTREZID %in% set9.node.panther.TF.related.entrezid] <- "TFrelated"
set9.node.annotation$TFrelated[is.na(set9.node.annotation$TFrelated)] <- "nonTFrelated"
set9.node.annotation$type[set9.node.annotation$estimate < 0] <- "up"
set9.node.annotation$type[set9.node.annotation$estimate > 0] <- "down"

set9.node.annotation <- left_join(set9.node.annotation, set9.node.panther.class.term.merged, by="ENTREZID")

set9.node.annotation <- set9.node.annotation %>% dplyr::select(-estimate1, -estimate2, -parameter, -conf.low, -conf.high, -method, -alternative, -adjusted.method, -Probe)

tmp <- set9.node.annotation %>% filter(!is.na(ENTREZID) )%>%  split(.$ENTREZID) %>% map(function(x) arrange(x,desc(adjusted.p))) %>% map(~tail(.,n=1))
set9.node.merged <- tmp %>% reduce(bind_rows)

set9.node.merged <- set9.node.merged[!is.na(set9.node.merged$ENTREZID),]

set9.node.merged.GO <- AnnotationDbi::select(org.Hs.eg.db,
                                                             keys = as.character(set9.node.merged$ENTREZID),
                                                             columns = "GO",
                                                             keytype = "ENTREZID")

set9.node.merged.GO.CC <- set9.node.merged.GO %>%  filter(EVIDENCE != "IEA" & ONTOLOGY == "CC")

set9.node.merged.GO.CC <- set9.node.merged.GO.CC %>% dplyr::rename(GOID=GO)
set9.node.merged.GO.CC <- set9.node.merged.GO.CC %>%  filter(EVIDENCE %in% c("IDA", "EXP", "IPI", "IMP","IGI","IEP"))
set9.node.merged.GO.CC <- unique(set9.node.merged.GO.CC)

total.set9.GO.annotation <- AnnotationDbi::select(GO.db,
                                                  keys = as.character(set9.node.merged.GO.CC$GO),
                                                  columns = c("GOID","DEFINITION","ONTOLOGY","TERM"),
                                                  keytypes = "GOID")


set9.node.merged.GO.CC <- left_join(set9.node.merged.GO.CC, total.set9.GO.annotation[,c(1,4)], by="GOID")


set9.node.merged.GO.CC <- set9.node.merged.GO.CC %>% dplyr::select(-EVIDENCE) %>% unique()


output.set9.go.term <- left_join(set9.node.merged[,c(5,6)], set9.node.merged.GO.CC[,c(1,2,4)], by="ENTREZID")
output.set9.go.term <- output.set9.go.term %>%  dplyr::rename(Id=name)
output.set9.go.term <- output.set9.go.term %>% filter(!is.na(GOID))
write_delim(output.set9.go.term, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set9/set9_GO_term.txt", delim="\t")
write_delim(set9.node.annotation, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set9/set9_node.txt", delim="\t")



# topology analysis -------------------------------------------------------

set9.ppi.topology <- read_csv(file="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set9/cytohubba.csv")
set9.ppi.topology.transform <- set9.ppi.topology %>%  mutate(Degree=log2(Degree+1), Betweenness=log2(Betweenness+1), BottleNeck=log2(BottleNeck+1))

ggplot(data=set9.ppi.topology.transform) + geom_density(aes(x=Degree))
ggplot(data=set9.ppi.topology.transform) + geom_density(aes(x=Betweenness))
ggplot(data=set9.ppi.topology.transform) + geom_density(aes(x=BottleNeck))
ggplot(data=set9.ppi.topology) + geom_density(aes(x=Degree))
ggplot(data=set9.ppi.topology) + geom_density(aes(x=Betweenness))
ggplot(data=set9.ppi.topology) + geom_density(aes(x=BottleNeck))




# subnetwork --------------------------------------------------------------
top_100_betweenness <- read_csv(file="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set9/top_100_betweenness_set9_nodes.csv")
top_100_bottleneck  <- read_csv(file="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set9/top_100_bottleneck_set9_nodes.csv")
top_100_degree      <- read_csv(file="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set9/top_100_degree_set9_nodes.csv")

top_set9_filter <- bind_rows(top_100_betweenness, top_100_bottleneck, top_100_degree) %>%  unique
sub.top_set9_filter <- top_set9_filter %>%  dplyr::select(Betweenness, BottleNeck, Degree, ENTREZID) %>% mutate(ENTREZID=as.character(ENTREZID))
# pre-filter




topologic_filter_set9_symbol <- unique(c(top_100_betweenness$`shared name`, top_100_bottleneck$`shared name`, top_100_degree$`shared name`) )


# edge

topologic.filter.clean.set9.005.total.ppi.output <- filter.clean.set9.005.total.ppi.output %>% filter(A %in% topologic_filter_set9_symbol & B %in% topologic_filter_set9_symbol)
write_delim(topologic.filter.clean.set9.005.total.ppi.output, path ="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set9/topologic_filter_interaction_edge.txt", delim = "\t")

# node annotation
topologic.filter.output.set9.go.term <- output.set9.go.term %>% filter(Id %in% topologic_filter_set9_symbol)
write_delim(topologic.filter.output.set9.go.term, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set9/topologic_filter_set9_GO_term.txt", delim="\t")


topologic.filter.set9.node.annotation <- set9.node.annotation %>%  filter(name %in% topologic_filter_set9_symbol)
topologic.filter.set9.node.annotation <- left_join(topologic.filter.set9.node.annotation,sub.top_set9_filter, by="ENTREZID" )
write_delim(topologic.filter.set9.node.annotation, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set9/topologic_filter_set9_node.txt", delim="\t")


# set10 network-------------------------------------------------------------------

# raw interaction query data
load("~/Documents/Repository/code_in_lab/set10_005_total_symbol_interaction.Rdata")
load("~/Documents/Repository/code_in_lab/set10_005_cortest.Rdata")
clean.set10.005.total.ppi.reactome <- clean.set10.005.total.ppi %>%  filter(sourceDatabases == "psi-mi:MI:0467(reactome)")
clean.set10.005.total.ppi.reactome <- clean.set10.005.total.ppi.reactome %>%  mutate(A_B=paste0(A,"-",B))
tmp <- set10.005.interaction.cor.test %>% mutate(A_B=paste0(A,"-",B))
tmp <- tmp %>% filter(A_B %in% clean.set10.005.total.ppi.reactome$A_B)
# method 1
# 24119
clean.set10.005.total.ppi.output <-  set10.005.interaction.cor.test %>% dplyr::select(A,B,estimate, p.value) %>% unique
# method 2
clean.set10.005.total.ppi.output <- tmp %>% dplyr::select(A,B,estimate, p.value) %>% unique
# 1609
# method 2.1
filter.clean.set10.005.total.ppi.output <- clean.set10.005.total.ppi.output %>%  filter(p.value < 0.05)
# method 2.2
filter.clean.set10.005.total.ppi.output <- clean.set10.005.total.ppi.output
write_delim(filter.clean.set10.005.total.ppi.output, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/set10_interaction_edge_FI_005.txt")

load("~/Documents/Repository/code_in_lab/hgu133plus2_annotation_v1.RData")
hgu133plus2.annotate.set <- hgu133plus2.probe.annotate %>% dplyr::select(name, Probe, ENTREZID)


set10.final.ppi.symbol.list<- unique(c(filter.clean.set10.005.total.ppi.output$A, filter.clean.set10.005.total.ppi.output$B))

set10.node.annotation <- set10.Sphere_CLS %>%  filter(name %in% set10.final.ppi.symbol.list)
set10.node.annotation <- set10.node.annotation %>% dplyr::select(-case, -lwr, -upr,-p.adj,-annova.p.value, -Probe)
tmp <- set10.node.annotation %>% filter(!is.na(ENTREZID) )%>%  split(.$ENTREZID) %>% map(function(x) arrange(x,desc(BH.p.adj))) %>% map(~tail(.,n=1))
set10.node.annotation <- tmp %>% reduce(bind_rows)

# functional annotation: search from the panther database
set10.node.panther.term <- search_panther_database(set10.node.annotation$ENTREZID)
set10.node.panther.class.term <- set10.node.panther.term %>%  filter(TERM_TYPE == "CLASS_TERM")
set10.node.panther.TF.related <- set10.node.panther.class.term %>% filter(str_detect(PANTHER_TERM,"transcription"))
set10.node.panther.TF.related.entrezid <- set10.node.panther.TF.related$ENTREZ %>% unique()

set10.node.panther.class.term.merged <- set10.node.panther.class.term %>%  group_by(ENTREZ) %>% summarise(PATHER_TERM=paste(PANTHER_TERM, collapse=","))
set10.node.panther.class.term.merged <- set10.node.panther.class.term.merged %>% dplyr::rename(ENTREZID = ENTREZ)
# add annotation on
set10.node.annotation$TFrelated[set10.node.annotation$ENTREZID %in% set10.node.panther.TF.related.entrezid] <- "TFrelated"
set10.node.annotation$TFrelated[is.na(set10.node.annotation$TFrelated)] <- "nonTFrelated"
set10.node.annotation$type[set10.node.annotation$diff > 0] <- "up"
set10.node.annotation$type[set10.node.annotation$diff < 0] <- "down"

set10.node.annotation <- left_join(set10.node.annotation, set10.node.panther.class.term.merged, by="ENTREZID")


set10.node.merged <- set10.node.annotation[!is.na(set10.node.annotation$ENTREZID),]

set10.node.merged <- set10.node.merged[!is.na(set10.node.merged$ENTREZID),]

set10.node.merged.GO <- AnnotationDbi::select(org.Hs.eg.db,
                                              keys = as.character(set10.node.merged$ENTREZID),
                                              columns = "GO",
                                              keytype = "ENTREZID")

set10.node.merged.GO.CC <- set10.node.merged.GO %>%  filter(EVIDENCE != "IEA" & ONTOLOGY == "CC")

set10.node.merged.GO.CC <- set10.node.merged.GO.CC %>% dplyr::rename(GOID=GO)
set10.node.merged.GO.CC <- set10.node.merged.GO.CC %>%  filter(EVIDENCE %in% c("IDA", "EXP", "IPI", "IMP","IGI","IEP"))
set10.node.merged.GO.CC <- unique(set10.node.merged.GO.CC)
library(GO.db)
total.set10.GO.annotation <- AnnotationDbi::select(GO.db,
                                                   keys = as.character(set10.node.merged.GO.CC$GO),
                                                   columns = c("GOID","DEFINITION","ONTOLOGY","TERM"),
                                                   keytypes = "GOID")


set10.node.merged.GO.CC <- left_join(set10.node.merged.GO.CC, total.set10.GO.annotation[,c(1,4)], by="GOID")


set10.node.merged.GO.CC <- set10.node.merged.GO.CC %>% dplyr::select(-EVIDENCE) %>% unique()
output.set10.go.term <- left_join(set10.node.merged[,c(3,4)], set10.node.merged.GO.CC[,c(1,2,4)], by="ENTREZID")
output.set10.go.term <- output.set10.go.term %>%  dplyr::rename(Id=name)
output.set10.go.term <- output.set10.go.term %>% filter(!is.na(GOID))
#  Annotation the pathway information -------------------------------------
# Input GSEA --------------------------------------------------------------
library(GSEABase)
ReactomePathway.geneset <- getGmt("/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/ReactomePathways.gmt")

ReactomePathway.title  <- ReactomePathway.geneset %>% map(setName)
ReactomePathway.Id     <- ReactomePathway.geneset %>% map(description)
ReactomePathway.symbol <- ReactomePathway.geneset %>% map(geneIds)


ReactomePathway.symbol <- ReactomePathway.symbol %>% map(~discard(.,function(x) x == "Reactome Pathway"))

ReactomePathway.sym.df <- ReactomePathway.symbol %>%
    map(. %>% data.frame) %>%
    list(.,ReactomePathway.title) %>%
    pmap(~mutate(.x, Title=.y)) %>%
    list(.,ReactomePathway.Id) %>%
    pmap(~mutate(.x, ReactomeId=.y))
ReactomePathway.sym.df <- ReactomePathway.sym.df %>% reduce(bind_rows)
colnames(ReactomePathway.sym.df) <- c("Symbol","PathwayTitle","ReactomeId")

ReactomePathwaySymbolNum <- ReactomePathway.sym.df %>% group_by(PathwayTitle) %>% summarise(SymbolNum=n())
ReactomePathway.sym.summarise.df <- ReactomePathway.sym.df %>%
    group_by(Symbol) %>%
    summarise(Title=paste(PathwayTitle,collapse=","), ReactomeId=paste(ReactomeId, collapse=","))




write_delim(output.set10.go.term, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/set10_GO_term_fi.txt", delim="\t")
write_delim(set10.node.annotation, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/set10_node_fi.txt", delim="\t")


set10.node.annotation <- set10.node.annotation %>% dplyr::rename(Symbol=name)
set10.node.annotation.addpathway <- left_join(set10.node.annotation, ReactomePathway.sym.summarise.df, by="Symbol")


# Add the intermediate item of pathway ------------------------------------

filter.clean.set10.005.total.ppi.output <- clean.set10.005.total.ppi.output
# filter the pathway with significant
ReactomePathway.sym.df.filter <- ReactomePathway.sym.df %>% filter(Symbol %in% set10.final.ppi.symbol.list) 
ReactomePathway_sig <- c(set10.Sphere_CLS.spia.2000.pathway.id$X1,
                         "R-HSA-937061", "R-HSA-937061", "R-HSA-381426", "R-HSA-8978868", "R-HSA-8979227", "R-HSA-74182",
                         "R-HSA-556833", "R-HSA-109688", "R-HSA-6791226", "R-HSA-69202")
# 14771
ReactomePathway.sym.df.filter <- ReactomePathway.sym.df.filter %>% filter(ReactomeId %in% ReactomePathway_sig)
ReactomePathway.sym.df.filter <- left_join(ReactomePathway.sym.df.filter, ReactomePathwaySymbolNum, by="PathwayTitle")

# secondary filter
# filter with large gene set
ReactomePathway.sym.df.filter <- ReactomePathway.sym.df.filter %>% filter(SymbolNum < 400)

# filter not relevant pathway
ReactomePathway.sym.df.filter <- ReactomePathway.sym.df.filter %>% filter(!ReactomeId %in% c("R-HSA-5663202", "R-HSA-72203", "R-HSA-72766",
                                                                                             "R-HSA-72312","R-HSA-3247509", "R-HSA-4839726","R-HSA-72172",
                                                                                             "R-HSA-162909", "R-HSA-1500931", "R-HSA-202403", "R-HSA-2500257",
                                                                                             "R-HSA-381340", "R-HSA-380287", "R-HSA-72187", "R-HSA-429914"))

ReactomePathway.sym.df.filter.tmp <- ReactomePathway.sym.df.filter %>%
                                        dplyr::select(Symbol, PathwayTitle) %>%
                                        dplyr::mutate(estimate=0,p.value=1) %>% 
                                        dplyr::rename(A=Symbol, B=PathwayTitle)
# 5267
write_delim(ReactomePathway.sym.df.filter.tmp, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/set10_only_pathway_gene_005.txt", delim="\t")

set10_interaction_add_pathway <- bind_rows(ReactomePathway.sym.df.filter.tmp, filter.clean.set10.005.total.ppi.output)


write_delim(set10_interaction_add_pathway, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/set10_interaction_edge_FI_pathway_005.txt", delim="\t")



# Topological module
# Prepare module gene set
set10.module <- read.csv("/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/annotated_set10_node_with_module.csv")
set10.module <- as.data.frame(set10.module) %>% dplyr::select(AverageShortestPathLength, 
                                                                BetweennessCentrality, 
                                                                ClosenessCentrality,
                                                                ClusteringCoefficient,
                                                                Eccentricity,
                                                                Degree,
                                                                module,
                                                                name)

set10.module <- set10.module %>% dplyr::rename(Symbol=name)
set10.node.annotation.addpathway.module <- left_join(set10.node.annotation.addpathway,set10.module, by="Symbol")
set10.node.annotation.addpathway.module$module[is.na(set10.node.annotation.addpathway.module$module)] <- 24
# category the node
set10.node.annotation.addpathway.module$Actor[set10.node.annotation.addpathway.module$module != 24 & !is.na(set10.node.annotation.addpathway.module$Title)] <-  "player"
set10.node.annotation.addpathway.module$Actor[set10.node.annotation.addpathway.module$module != 24 &  is.na(set10.node.annotation.addpathway.module$Title)] <-  "interactor"
set10.node.annotation.addpathway.module$Actor[set10.node.annotation.addpathway.module$module == 24 & !is.na(set10.node.annotation.addpathway.module$Title)] <-  "unknown interactor"
set10.node.annotation.addpathway.module$Actor[set10.node.annotation.addpathway.module$module == 24 & is.na(set10.node.annotation.addpathway.module$Title)] <-   "unknown"

write_delim(set10.node.annotation.addpathway.module, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/set10_interaction_node_FI_pathway_module_005.txt",delim="\t")

# Use pathway item to adjust the functional interaction -------------------



# function for enrichment
getModuleGprofiler <- function(set10_pathwayadjusted_module, module_numb){
    
    set10.module0 <- set10_pathwayadjusted_module %>% filter(module == module_numb) 
    query <- set10.module0 %>% filter(!is.na(Actor)) %>% with(name)
    set10.module0.enrichment<-  gprofiler(query,
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
    return(set10.module0.enrichment)
}


# Module 0 ----------------------------------------------------------------

set10.module0 <- set10_pathwayadjusted_module %>% filter(module == 0) 
query <- set10.module0 %>% filter(!is.na(Actor)) %>% with(name)
set10.module0.enrichment<-  gprofiler(query,
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


# Module 1 ----------------------------------------------------------------

set10.module1 <- set10_pathwayadjusted_module %>% filter(module == 1) 
set10.module1.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 1)

# Module 2 ----------------------------------------------------------------
set10.module2 <- set10_pathwayadjusted_module %>% filter(module == 2) 
set10.module2.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 2)


# Module 3 ----------------------------------------------------------------
set10.module3 <- set10_pathwayadjusted_module %>% filter(module == 3) 
set10.module3.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 3)


# Module 4 ----------------------------------------------------------------
set10.module4 <- set10_pathwayadjusted_module %>% filter(module == 4) 
set10.module4.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 4)


# Module 5 ----------------------------------------------------------------
set10.module5 <- set10_pathwayadjusted_module %>% filter(module == 5) 
set10.module5.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 5)


# Module 6 ----------------------------------------------------------------
set10.module6 <- set10_pathwayadjusted_module %>% filter(module == 6) 
set10.module6.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 6)


# Module 7 ----------------------------------------------------------------
set10.module7 <- set10_pathwayadjusted_module %>% filter(module == 7) 
set10.module7.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 7)


# Module 8 ----------------------------------------------------------------
set10.module8 <- set10_pathwayadjusted_module %>% filter(module == 8) 
set10.module8.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 8)


# Module 9 ----------------------------------------------------------------
set10.module9 <- set10_pathwayadjusted_module %>% filter(module == 9) 
set10.module9.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 9)


# Module 10 ---------------------------------------------------------------
set10.module10 <- set10_pathwayadjusted_module %>% filter(module == 10) 
set10.module10.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 10)


# Module 11 ---------------------------------------------------------------
set10.module11 <- set10_pathwayadjusted_module %>% filter(module == 11) 
set10.module11.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 11)


# Module 12 ---------------------------------------------------------------
set10.module12 <- set10_pathwayadjusted_module %>% filter(module == 12) 
set10.module12.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 12)


# Module 13 ---------------------------------------------------------------
set10.module13 <- set10_pathwayadjusted_module %>% filter(module == 13) 
set10.module13.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 13)


# Module 14 ---------------------------------------------------------------
set10.module14 <- set10_pathwayadjusted_module %>% filter(module == 14) 
set10.module14.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 14)


# Module 15 ---------------------------------------------------------------
set10.module15 <- set10_pathwayadjusted_module %>% filter(module == 15) 
set10.module15.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 15)


# Module 16 ---------------------------------------------------------------
set10.module16 <- set10_pathwayadjusted_module %>% filter(module == 16) 
set10.module16.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 16)


# Module 17 ---------------------------------------------------------------
set10.module17 <- set10_pathwayadjusted_module %>% filter(module == 17) 
set10.module17.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 17)


# Module 18 ---------------------------------------------------------------
set10.module18 <- set10_pathwayadjusted_module %>% filter(module == 18) 
set10.module18.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 18)


# Module 19 ---------------------------------------------------------------
set10.module19 <- set10_pathwayadjusted_module %>% filter(module == 19) 
set10.module19.enrichment<-  getModuleGprofiler(set10_pathwayadjusted_module, module_numb = 19)


# Module 20 ---------------------------------------------------------------




# topologic calculate ------------------------------------------------------


set10.ppi.topology <- read_csv(file="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/set10_hubba_export.csv")
set10.ppi.topology.transform <- set10.ppi.topology %>%  mutate(Degree=log2(Degree+1), Betweenness=log2(Betweenness+1), BottleNeck=log2(BottleNeck+1))

ggplot(data=set10.ppi.topology.transform) + geom_density(aes(x=Degree))
ggplot(data=set10.ppi.topology.transform) + geom_density(aes(x=Betweenness))
ggplot(data=set10.ppi.topology.transform) + geom_density(aes(x=BottleNeck))

ggplot(data=set10.ppi.topology) + geom_density(aes(x=Degree))
ggplot(data=set10.ppi.topology) + geom_density(aes(x=Betweenness))
ggplot(data=set10.ppi.topology) + geom_density(aes(x=BottleNeck))

