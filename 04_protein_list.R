############################################################################
#     The Analysis for the  proteomic data result
# 
#
############################################################################

library(readxl)
raw <- read_excel("/Users/Weitinglin/Documents/R_scripts/Lab/Proteomic/20110310_Lung_cancer_stem_cell_normalized_single.xlsx")
success_id_gene <- raw$`Gene Symbol`
raw <- raw[,c(1:4,6:9,11:32)]


# result gene list --------------------------------------------------------
load("~/Documents/Repository/code_in_lab/updown_tukey_aov_result.Rdata")
## Upregulation
cutoff.p <- 0.0001
up.CLF_CLS    <- up.tukey.aov.CLF_CLS %>% filter(BH.p.adj < cutoff.p)
up.Sphere_CLS <- up.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)


# Downregulation
cutoff.p <- 0.0001
down.CLF_CLS    <- down.tukey.aov.CLF_CLS  %>% filter(BH.p.adj < cutoff.p)
down.Sphere_CLS <- down.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)

# _down.stream analysis ---------------------------------------------------
# GET THE Entrezid
union.up       <- union(up.CLF_CLS$EntrezID, up.Sphere_CLS$EntrezID)
intersect.up   <- intersect(up.CLF_CLS$EntrezID, up.Sphere_CLS$EntrezID)

union.down     <- union(down.CLF_CLS$EntrezID, down.Sphere_CLS$EntrezID)
intersect.down <- intersect(down.CLF_CLS$EntrezID, down.Sphere_CLS$EntrezID)

# Get the Probe
union.up       <- union(up.CLF_CLS$Probe, up.Sphere_CLS$Probe)
intersect.up   <- intersect(up.CLF_CLS$Probe, up.Sphere_CLS$Probe)

union.down     <- union(down.CLF_CLS$Probe, down.Sphere_CLS$Probe)
intersect.down <- intersect(down.CLF_CLS$Probe, down.Sphere_CLS$Probe)


# _down.stream analysis ---------------------------------------------------
# Hierachial cluster use heatmaply
union.up.symbol       <- union(up.CLF_CLS$Symbol, up.Sphere_CLS$Symbol)
intersect.up.symbol   <- intersect(up.CLF_CLS$Symbol, up.Sphere_CLS$Symbol)

union.down.symbol     <- union(down.CLF_CLS$Symbol, down.Sphere_CLS$Symbol)
intersect.down.symbol <- intersect(down.CLF_CLS$Symbol, down.Sphere_CLS$Symbol)



hit_protein.union.up   <- union.up.symbol[!is.na(union.up.symbol)][(union.up.symbol[!is.na(union.up.symbol)] %in% success_id_gene %>% which)]
hit_protein.union.down   <- union.down.symbol[!is.na(union.down.symbol)][(union.down.symbol[!is.na(union.down.symbol)] %in% success_id_gene %>% which)]


hit_protein.intersect.up   <- intersect.up.symbol[!is.na(intersect.up.symbol)][(intersect.up.symbol[!is.na(intersect.up.symbol)] %in% success_id_gene %>% which)]
hit_protein.intersect.down   <- intersect.down.symbol[!is.na(intersect.down.symbol)][(intersect.down.symbol[!is.na(intersect.down.symbol)] %in% success_id_gene %>% which)]

hit.data.frame <- raw %>% filter(`Gene Symbol` %in% c(hit_protein.intersect.up,hit_protein.intersect.down))
hit.data.frame$Regulation[hit.data.frame$`Gene Symbol` %in% hit_protein.intersect.up]   <- "Up_Regulation"
hit.data.frame$Regulation[hit.data.frame$`Gene Symbol` %in% hit_protein.intersect.down] <- "Down_Regulation"
hit.data.frame <- hit.data.frame[,c(31,1:30)]
dim(hit.data.frame)
write.csv(hit.data.frame, file="hit_ID_memebrane_protein_in_microarray_gene_list.csv")

