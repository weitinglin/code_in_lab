# 07 survival analysis
# take the network analysis result modules to see whether can seperate the survival curve


# import library
library(survival)
library(OIsurv)
library(affy)
library(readr)
library(dplyr)

# Input validate sample ---------------------------------------------------

# GSE19188 ----------------------------------------------------------------




# GSE31210 ----------------------------------------------------------------
GSE31210.phenotype.information <- read_delim(file ="/Volumes/2016weiting/lab_data/validationdataset/Raw_data/GSE31210/GSE31210_series_matrix_cut.txt",
                                    delim = "\t",
                                    n_max = 33)
total.survival.information <- read_delim(file="/Volumes/2016weiting/lab_data/validationdataset/20160817_PD1_PDL1_Clinica1.txt", delim = "\t") 

GSE31210.survial.information <- total.survival.information %>% filter(Study == 4)

GSE31210.mutation.status <- GSE31210.phenotype.information[1:33,] 

data.path <- "/Volumes/2016weiting/lab_data/validationdataset/Raw_data/GSE31210/GSE31210_RAW"
experiment.set <- GSE31210.phenotype.information[2,] %>%  unlist %>%  unname
experiment.set <- experiment.set[-1]
phenodata.set   <- matrix ( rep ( experiment.set, 2) , ncol = 2 )
phenodata.set   <- as.data.frame ( phenodata.set )
colnames ( phenodata.set )   <- c ( "Name" , "FileName" )
phenodata.set$experiment.set <- experiment.set
write.table(phenodata.set, paste(data.path,"/phenodata_set.txt", sep = ""),quote=F,sep="\t",row.names=F)

celfile.set <- ReadAffy ( celfile.path = data.path ,
                          phenoData = paste(data.path,"/phenodata_set.txt", sep = ""),
                          compress = T)
save(celfile.set, file="GSE312310celfile.Rdata")
Exprs.data <- mas5 ( celfile.set , normalize = FALSE, analysis = "absolute", sc = 500)
save(Exprs.data, file="GSE31210exprsdata.Rdata")

ecelfile.set <- exprs(Exprs.data)
rm(Exprs.data)
norm.GSE31210 <- quantile_normalization(ecelfile.set, method = "median")
rm(ecelfile.set)
save(norm.GSE31210, file="QNnormGSE31210.Rdata")

load("~/Documents/Repository/code_in_lab/GSE312310norm.Rdata")
load("~/Documents/Repository/code_in_lab/hgu133plus2_annotation_v1.RData")

hgu133plus2.annotate.set <- hgu133plus2.probe.annotate %>% dplyr::select(name, Probe, ENTREZID)
# remove probe without annotation 
hgu133plus2.annotate.set <- hgu133plus2.annotate.set[!is.na(hgu133plus2.annotate.set$ENTREZID),]
probe.list <- rownames(norm.GSE312310) 
filter.index <- probe.list %in% hgu133plus2.annotate.set$Probe 

# first filter
norm.GSE312310 <- norm.GSE312310[filter.index,]

# get mean of the same gene
norm.probe.GSE312310 <- as.data.frame(norm.GSE312310)
norm.probe.GSE312310$Probe <- rownames(norm.probe.GSE312310)
norm.probe.GSE312310 <- left_join(norm.probe.GSE312310, hgu133plus2.annotate.set, by="Probe")
mean.norm.probe.GSE312310 <- norm.probe.GSE312310 %>%
                                    dplyr::select(-Probe) %>%
                                    group_by(name, ENTREZID) %>% 
                                    summarise_each(funs(mean))




# Create the normal lung tissue expression level --------------------------
reference.expression.GSE312310 <- mean.norm.probe.GSE312310[,229:248]

reference.expression.mean.GSE312310 <- reference.expression.GSE312310 %>% t() %>% 
                                            as.data.frame()  %>%
                                             summarise_all(funs(median)) 

reference.expression.mean.GSE312310 <- reference.expression.mean.GSE312310 %>% t
colnames(reference.expression.mean.GSE312310) <- "normal"

tumor.expression.GSE312310 <- mean.norm.probe.GSE312310[,1:228]
total.expression.GSE312310 <- bind_cols(tumor.expression.GSE312310,as.data.frame(reference.expression.mean.GSE312310))
# Prepare the data for GSEA analysis

gsea.prepare.data.path <- "/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/GSEAgct_GSE31210/"

norm.name        <- rownames(norm.tmp)
norm.description <- rep("na", nrow(norm.tmp))
tumor.index <- 1:60
normal.index <- 61:120
for ( i in 1:60){
    norm.tmp <- as.data.frame(norm.GSE19804[,c(tumor.index[i],normal.index[i])])
    norm.tmp$NAME <- norm.name
    norm.tmp$DESCRIPTION <- norm.description
    colnames(norm.tmp) <- c("Tumor","Normal", "Name", "Description")
    norm.tmp <- norm.tmp[,c(3,4,1,2)]
    write.table(norm.tmp, file=paste0(gsea.prepare.data.path,"tmp_norm_qn_",i,".gct"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
}
# GSE19804 ----------------------------------------------------------------


# biosample code GSE19804 120 .cel files

# file code GSM494556 - GSM494675 	

# Phenotype information
sample.code <- list(as.list("GSM"),as.list(494556:494675)) %>%  pmap(paste0) %>% unlist

phenotype.information <- read_delim(file = "/Volumes/2016weiting/lab_data/validationdataset/Raw_data/GSE19804/GSE19804_series_matrix_cut.txt",
                                    delim = "\t",
                                    n_max = 33)

total.survival.information <- read_delim(file="/Volumes/2016weiting/lab_data/validationdataset/20160817_PD1_PDL1_Clinica1.txt", delim = "\t") 
GSE19804.survial.information <- total.survival.information %>% filter(Study == 10)

# Preparation of the expression


data.path  <- file.path("")
experiment.set <- phenotype.information[1,] %>%  unlist %>%  unname
experiment.set <- experiment.set[-1]
phenodata.set   <- matrix ( rep ( set, 2) , ncol = 2 )
phenodata.set   <- as.data.frame ( phenodata.set )
colnames ( phenodata.set )   <- c ( "Name" , "FileName" )
phenodata.set$experiment.set <- experiment.set
write.table(phenodata.set, paste(data.path,"/phenodata_set.txt", sep = ""),quote=F,sep="\t",row.names=F)

# read file
celfile.set <- ReadAffy ( celfile.path = data.path ,
                          phenoData = paste(data.path,"/phenodata_set.txt", sep = ""),
                          compress = T)

Exprs.data <- mas5 ( celfile.set , normalize = FALSE, analysis = "absolute", sc = 500)
save(Exprs.data, file="GSE19804exprsdata.Rdata")

ecelfile.set <- exprs(Exprs.data)

norm.GSE19804 <- quantile_normalization(ecelfile.set, method = "median")
save(norm.GSE19804, file="QNnormGSE19803.Rdata")

# Prepare the data for GSEA analysis

gsea.prepare.data.path <- "/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/GSEAgct/"

norm.name        <- rownames(norm.tmp)
norm.description <- rep("na", nrow(norm.tmp))
tumor.index <- 1:60
normal.index <- 61:120
for ( i in 1:60){
norm.tmp <- as.data.frame(norm.GSE19804[,c(tumor.index[i],normal.index[i])])
norm.tmp$NAME <- norm.name
norm.tmp$DESCRIPTION <- norm.description
colnames(norm.tmp) <- c("Tumor","Normal", "Name", "Description")
norm.tmp <- norm.tmp[,c(3,4,1,2)]
write.table(norm.tmp, file=paste0(gsea.prepare.data.path,"tmp_norm_qn_",i,".gct"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
}



# Prepare module gene set
# this gene set was the first version with the reactome FI network clustering result
set10.module <- read.csv("/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/annotated_set10_node_with_module.csv")
set10.module <- as.data.frame(set10.module)

set10.module.geneset <- set10.module %>% group_by(module) %>% summarise(geneset=paste(shared.name, collapse=",")) %>% mutate(description="na")
set10.module.geneset <- set10.module.geneset %>% mutate(all=paste(module,description,geneset,sep=","))
set10.module.geneset <- set10.module.geneset %>% dplyr::select(all)
write.table(set10.module.geneset, file=paste0(gsea.prepare.data.path,"module.gmt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

# this gene set was use the FI network with enrichment pathway for adjustment then use FI network clustering
gsea.prepare.data.path <- "/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/GSEAgct_two/"

set10_pathwayadjusted_module <- read_csv(file="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/set10_interaction_pathwayadjusted_module.csv") %>% dplyr::select(-SUID, -selected, -`shared name`)
set10_pathwayadjusted_module.symbol <- set10_pathwayadjusted_module %>% filter(!is.na(Actor))
set10_pathwayadjusted_module.symbol.geneset <- set10_pathwayadjusted_module.symbol %>% group_by(module) %>%
                                                    summarise(geneset=paste(name, collapse=",")) %>% 
                                                    mutate(description="na")
set10_pathwayadjusted_module.symbol.geneset <- set10_pathwayadjusted_module.symbol.geneset %>% mutate(all=paste(module,description,geneset,sep=",")) %>% dplyr::select(all)
write.table(set10_pathwayadjusted_module.symbol.geneset, file=paste0(gsea.prepare.data.path,"set10_pathadjust_module.gmt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

set10_pathwayadjusted_module$Actor[is.na(set10_pathwayadjusted_module$Actor) ] <- "Pathway"
set10_pathwayadjusted_module$type[is.na(set10_pathwayadjusted_module$type)] <- "none"
write_delim(set10_pathwayadjusted_module, path="/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/set10_FI_pathway_module_fixednoeannotation_005.txt", delim="\t")
# input the GSE19804 GSEA result
library(readxl)

for (i in 1:60){
file_path <- paste0("/Users/Weitinglin/Documents/R_scripts/Lab/GSEA/set10module/sample_",i,"_weighted_s2n.Gsea")
t.file.list <- list.files(path=file_path, full.names = TRUE, pattern = "gsea_report_for_Tumor_")
n.file.list <- list.files(path=file_path, full.names = TRUE, pattern = "gsea_report_for_Normal_")
t.xlsfile <- t.file.list[2]
n.xlsfile <- n.file.list[2]
n.sample.name <- paste0("tmp_normal_sample_",i)
t.sample.name <- paste0("tmp_tumor_sample_",i)
sample.name <- paste0("sample_",i)
assign(n.sample.name, read_delim(file=n.xlsfile, delim="\t"))
assign(t.sample.name, read_delim(file=t.xlsfile, delim="\t"))
assign(sample.name, bind_rows(get(n.sample.name),get(t.sample.name)))
}




# Combine GSEA and phenotype ----------------------------------------------
sampleGSE19804.gsea.list <- list()
for (i in 1:60){
    sampleGSE19804.gsea.list[[i]] <- get(paste0("sample_",i))
}

sampleGSE19804.gsea.list <- sampleGSE19804.gsea.list %>% map(. %>% arrange(NAME))

sampleGSE19804.es <- sampleGSE19804.gsea.list %>% map(. %>% dplyr::select(NES)) %>% reduce(bind_cols)
col.name.vectors <- list(as.list("sample_"), as.list(1:60)) %>% pmap(paste0) %>% unlist
colnames(sampleGSE19804.es) <- c(col.name.vectors)



# Make heatmap

library(RColorBrewer)
library(plotly)
library(heatmaply)
RdYlGn <- colorRampPalette(brewer.pal(11, "RdYlGn"))
heatmaply(sampleGSE19804.es, k_col = 4, k_row = 4,
          xaxis_font_size = "10pt", yaxis_font_size = "10pt",
          row_text_angle = -40, column_text_angle =  90,
          col_side_colors = c(as.numeric(GSE19804.survial.information$Stage[1:60])),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "green", high = "red", midpoint = 0) ) %>%
    layout(margin = list(l = 180, b = 100))


index.stage2.3 <- GSE19804.survial.information$Stage == 1 | GSE19804.survial.information$Stage == 3
index.stage2.3 <- index.stage2.3[1:60]
heatmaply(sampleGSE19804.es[,index.stage2.3], k_col = 4, k_row = 4,
          xaxis_font_size = "10pt", yaxis_font_size = "10pt",
          row_text_angle = -40, column_text_angle =  90,
          col_side_colors = c(as.numeric(GSE19804.survial.information$Stage[1:60][index.stage2.3])),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "green", high = "red", midpoint = 0) ) %>%
    layout(margin = list(l = 180, b = 100))


# use annova to select significant modules --------------------------------
library(foreach)
library(doParallel)
ex.design <- factor(GSE19804.survial.information$Stage[1:60])
test.annova.result  <- foreach(i = 1:nrow(sampleGSE19804.es)) %dopar% {
    aov(as.matrix(sampleGSE19804.es)[i,] ~ ex.design) 
}

test.tukey.result  <- foreach(i = 1:nrow(sampleGSE19804.es)) %dopar% {
    TukeyHSD(test.annova.result[[i]]) %>% broom::tidy()
}




sampleGSE19804.fdr <- sampleGSE19804.gsea.list %>% map(. %>% dplyr::select(`FDR q-val`)) %>% reduce(bind_cols)
col.name.vectors <- list(as.list("sample_"), as.list(1:60)) %>% pmap(paste0) %>% unlist
colnames(sampleGSE19804.fdr) <- c(col.name.vectors)





# Use GAGE for custom gene set enrichment ---------------------------------
library(gage)
module.geneset <- readList("/Users/Weitinglin/Documents/R_scripts/Lab/gephi/set10/GSEAgct_two/set10_trans_pathadjust_module.gmt")
data("egSymb")
modele.geneset.entrezid <- lapply(module.geneset, sym2eg)
#228 , 1
total.expression.GSE312310.matrix <- total.expression.GSE312310 %>% ungroup %>% dplyr::select(-name, -ENTREZID) %>% as.matrix
rownames(total.expression.GSE312310.matrix) <- total.expression.GSE312310$name
test.result <- list()

library(foreach)
library(doParallel)
registerDoParallel(cores=4)

test.result <- foreach(i=1:266) %dopar% {
    gage(total.expression.GSE312310.matrix, gsets = module.geneset, ref = c(227), samp = c(i))
}
test.result <- list()
for (i in 1:266){
    test.result[[i]] <- gage(total.expression.GSE312310.matrix, gsets = module.geneset, ref = c(227), samp = c(i))
}



