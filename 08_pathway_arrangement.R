# arrange the SPIA result
library(readr)
library(dplyr)





# 01_SPIAresult_summary ---------------------------------------------------



# Load reactome/kegg result -----------------------------------------------



set10.Sphere_CLS.spia.result.kegg           <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereCLS_kegg.csv", delim = ",")           %>% filter(pGFdr < 0.05)  %>%  mutate(name_dir=paste0(Name, "-", Status))
set10.Sphere_Timeseries.spia.result.kegg    <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereTimeseries_KEGG.csv",delim=",")       %>% filter(pGFdr < 0.05)   %>%  mutate(name_dir=paste0(Name, "-", Status))   

set10.Sphere_CLS.spia.result.Reactom        <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereCLS_reactome.csv", delim=",")         %>% filter(pGFdr < 0.05)  %>%  mutate(name_dir=paste0(Name, "-", Status))    
set10.Sphere_Timeseries.spia.result.Reactom <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereTimeseries_Reactome.csv", delim=",")  %>% filter(pGFdr < 0.05)  %>%  mutate(name_dir=paste0(Name, "-", Status))          



# _Kegg -------------------------------------------------------------------


set10.kegg.intersect.index <- intersect(set10.Sphere_CLS.spia.result.kegg$name_dir, set10.Sphere_Timeseries.spia.result.kegg$name_dir)
set10.Sphere_CLS.spia.result.kegg.filter <- set10.Sphere_CLS.spia.result.kegg %>%  filter(name_dir %in% set10.kegg.intersect.index)
set10.Sphere_Timeseries.spia.result.kegg.filter <- set10.Sphere_Timeseries.spia.result.kegg %>%  filter(name_dir %in% set10.kegg.intersect.index)


# _Reactome ---------------------------------------------------------------


set10.reactome.intersect.index <- intersect(set10.Sphere_CLS.spia.result.Reactom$name_dir, set10.Sphere_Timeseries.spia.result.Reactom$name_dir)
set10.Sphere_CLS.spia.result.Reactom.filter <- set10.Sphere_CLS.spia.result.Reactom %>% filter(name_dir %in% set10.reactome.intersect.index)
set10.Sphere_Timeseries.spia.result.Reactom.filter <- set10.Sphere_Timeseries.spia.result.Reactom %>% filter(name_dir %in% set10.reactome.intersect.index)



# read the reactome database schema in order to know the hierarchy --------

esemble2reactome_all_levels <- read_delim(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/Ensembl2Reactome_All_Levels.txt", col_names = FALSE, delim="\t")
esembl2reactome <- read_delim(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/Ensembl2Reactome.txt", col_names = FALSE, delim="\t")
ensemble2reactomereactoms <- read_delim(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/Ensembl2ReactomeReactions.txt", col_names = FALSE, delim="\t")


reactome.pathway.map          <- read_delim(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/ReactomePathways.txt", col_names = FALSE, delim="\t") %>% filter(X3 %in% "Homo sapiens")
reactome.pathway.relationship <- read_delim(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/ReactomePathwaysRelation.txt", col_names = FALSE, delim="\t")



# Gather the pathway identifier in our result
reactome.pathway.result <-    reactome.pathway.map %>% filter(X2 %in% set10.Sphere_CLS.spia.result.Reactom.filter$Name)

# origin 177 pathway, 13 cannot diectly mapped, add on 14 
manual.add.pathway.result <- reactome.pathway.map %>% filter(X1 %in% c("R-HSA-74160", "R-HSA-937061", "R-HSA-933541","R-HSA-975110","R-HSA-933542", "R-HSA-975110",
                                                                       "R-HSA-933542", "R-HSA-975138", "R-HSA-937072", "R-HSA-381426","R-HSA-74182", "R-HSA-8979227",
                                                                       "R-HSA-8978868", "R-HSA-6791226", "R-HSA-69202", "R-HSA-166166", "R-HSA-109688"))

reactome.pathway.result <- bind_rows(reactome.pathway.result, manual.add.pathway.result) %>% unique


reactome.pathway.result.relationship <- reactome.pathway.relationship %>%  filter(X1 %in% reactome.pathway.result$X1 & X2 %in% reactome.pathway.result$X1)

reactome.identifiers <- unique(c(reactome.pathway.result$X1))

reactome.pathway.result.relationship.nodes <- reactome.pathway.map %>% filter(X1 %in% reactome.identifiers)

reactome.pathway.result.relationship.nodes <- reactome.pathway.result.relationship.nodes %>%  dplyr::rename(ReactomeID=X1, Name=X2, Species=X3)

reactome.pathway.result.nodes.pvalue <- set10.Sphere_CLS.spia.result.Reactom.filter %>% dplyr::select(Name, pGFdr, pSize, NDE,Status)

reactome.pathway.result.relationship.nodes <- left_join(reactome.pathway.result.relationship.nodes, reactome.pathway.result.nodes.pvalue)

reactome.pathway.result.manual.nodes.pvalue <- reactome.pathway.result.relationship.nodes %>% filter(is.na(NDE))
reactome.pathway.result.relationship.nodes <- reactome.pathway.result.relationship.nodes %>% filter(!is.na(NDE))

reactome.pathway.result.manual.nodes.pvalue <- reactome.pathway.result.manual.nodes.pvalue %>% dplyr::select(ReactomeID, Name, Species)
manual.add.nodes.info <- data.frame(ReactomeID=c("R-HSA-74160","R-HSA-166166", "R-HSA-937061", "R-HSA-69202", "R-HSA-109688", "R-HSA-975110", "R-HSA-937072"),
                                                          pGFdr=c(8.761657e-23, 6.277138e-05, 6.277138e-05, 0.036525189,0.0108200970, 4.555800e-05, 4.555800e-05),
                                                          pSize=c(1505, 97, 97, 65, 61, 72, 72),
                                                          NDE=c(513, 45, 45, 25, 26, 36, 36))
reactome.pathway.result.manual.nodes.pvalue <- left_join(reactome.pathway.result.manual.nodes.pvalue, manual.add.nodes.info)
reactome.pathway.result.relationship.nodes <- bind_rows(reactome.pathway.result.relationship.nodes, reactome.pathway.result.manual.nodes.pvalue)



write_delim(reactome.pathway.result.relationship, path="/Users/Weitinglin/Documents/R_scripts/Lab/Cytoscape/reactome_spia_result_visualization/reactome_spia_set10_relation.txt", delim = "\t")
write_delim(reactome.pathway.result.relationship.nodes, path="/Users/Weitinglin/Documents/R_scripts/Lab/Cytoscape/reactome_spia_result_visualization/reactome_spia_set10_node.txt", delim = "\t")


origin.pathway.name <- c("Gene Expression","TRIF-mediated TLR3/TLR4 signaling")
new.pathway.name    <- c("Gene expression (Transcription)","TRIF-mediated TLR3/TLR4 signaling") 
#R-HSA-74160   Gene expression (Transcription)  
#R-HSA-937061  TRIF-mediated TLR3/TLR4 signaling
#TRAF6 Mediated Induction of proinflammatory cytokines ?
#R-HSA-933541                                               TRAF6 mediated IRF7 activation Homo sapiens
#R-HSA-975110                      TRAF6 mediated IRF7 activation in TLR7/8 or 9 signaling Homo sapiens
#R-HSA-933542                                              TRAF6 mediated NF-kB activation Homo sapiens
#R-HSA-975138 TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation Homo sapiens
#R-HSA-937072                                     TRAF6 mediated induction of TAK1 complex Homo sapiens
#R-HSA-381426  Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) 
#R-HSA-74182       Ketone body metabolism
#R-HSA-8979227   Triglyceride metabolism Homo sapiens
#R-HSA-8978868   Fatty acid metabolism Homo sapiens
#R-HSA-6791226 Major pathway of rRNA processing in the nucleolus and cytosol Homo sapiens
#R-HSA-69202                                        Cyclin E associated events during G1/S transition
#R-HSA-166166           MyD88-independent TLR3/TLR4 cascade
#R-HSA-109688 Cleavage of Growing Transcript in the Termination Region

#Gene Expression ok
#TRAF6 Mediated Induction of proinflammatory cytokines
#TRIF-mediated TLR3/TLR4 signaling ok
#MyD88-independent TLR3/TLR4 cascade ok
#Regulation of IGF Activity by IGFBP
#Fatty acid, triacylglycerol, and ketone body metabolism
#Metabolism of lipids and lipoproteins
#RNA Polymerase I, RNA Polymerase III, and Mitochondrial Transcription
#Post-Elongation Processing of Intron-Containing pre-mRNA
#Post-Elongation Processing of the Transcript
#Major pathway of rRNA processing in the nucleolus
#Cyclin E associated events during G1/S transition


setdiff(set10.Sphere_CLS.spia.result.Reactom.filter$Name,reactome.pathway.result$X2) 



reactome.pathway.map %>% filter(str_detect(X2,"Growing Transcript"))


# 02_SPIAresult_summary ---------------------------------------------------



# SPIA with permutation 2000 times ----------------------------------------
set10.Sphere_CLS.spia.2000.result.kegg           <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereCLS_2000_kegg.csv", delim = ",")           %>% filter(pGFdr < 0.05)  %>%  mutate(name_dir=paste0(Name, "-", Status))
set10.Sphere_Timeseries.spia.2000.result.kegg    <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereTimeseries_2000_KEGG.csv",delim=",")       %>% filter(pGFdr < 0.05)   %>%  mutate(name_dir=paste0(Name, "-", Status))   

set10.Sphere_CLS.spia.2000.result.Reactom        <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereCLS_2000_Reactome.csv", delim=",")         %>% filter(pGFdr < 0.05)  %>%  mutate(name_dir=paste0(Name, "-", Status))    
set10.Sphere_Timeseries.spia.2000.result.Reactom <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereTimeseries_2000_Reactome.csv", delim=",")  %>% filter(pGFdr < 0.05)  %>%  mutate(name_dir=paste0(Name, "-", Status))          


# _Kegg -------------------------------------------------------------------


set10.kegg.2000.intersect.index <- intersect(set10.Sphere_CLS.spia.2000.result.kegg$name_dir, set10.Sphere_Timeseries.spia.2000.result.kegg$name_dir)
set10.Sphere_CLS.spia.2000.result.kegg.filter <- set10.Sphere_CLS.spia.2000.result.kegg %>%  filter(name_dir %in% set10.kegg.2000.intersect.index)
set10.Sphere_Timeseries.spia.2000.result.kegg.filter <- set10.Sphere_Timeseries.spia.2000.result.kegg %>%  filter(name_dir %in% set10.kegg.2000.intersect.index)


# _Reactome ---------------------------------------------------------------


set10.reactome.2000.intersect.index <- intersect(set10.Sphere_CLS.spia.2000.result.Reactom$name_dir, set10.Sphere_Timeseries.spia.2000.result.Reactom$name_dir)
set10.Sphere_CLS.spia.2000.result.Reactom.filter <- set10.Sphere_CLS.spia.2000.result.Reactom %>% filter(name_dir %in% set10.reactome.2000.intersect.index)
set10.Sphere_Timeseries.spia.2000.result.Reactom.filter <- set10.Sphere_Timeseries.spia.2000.result.Reactom %>% filter(name_dir %in% set10.reactome.2000.intersect.index)



# Gather the pathway identifier in our result
reactome.pathway.result <- reactome.pathway.map %>% filter(X2 %in% set10.Sphere_CLS.spia.2000.result.Reactom.filter$Name)

# origin 177 pathway, 13 cannot diectly mapped, add on 14 
manual.add.pathway.result <- reactome.pathway.map %>% filter(X1 %in% c("R-HSA-74160", "R-HSA-937061", "R-HSA-933541","R-HSA-975110","R-HSA-933542", "R-HSA-975110",
                                                                       "R-HSA-933542", "R-HSA-975138", "R-HSA-937072", "R-HSA-381426","R-HSA-74182", "R-HSA-8979227",
                                                                       "R-HSA-8978868", "R-HSA-6791226", "R-HSA-69202", "R-HSA-166166", "R-HSA-109688"))

reactome.pathway.result <- bind_rows(reactome.pathway.result, manual.add.pathway.result) %>% unique


reactome.pathway.result.relationship <- reactome.pathway.relationship %>%  filter(X1 %in% reactome.pathway.result$X1 & X2 %in% reactome.pathway.result$X1)

reactome.identifiers <- unique(c(reactome.pathway.result$X1))

reactome.pathway.result.relationship.nodes <- reactome.pathway.map %>% filter(X1 %in% reactome.identifiers)

reactome.pathway.result.relationship.nodes <- reactome.pathway.result.relationship.nodes %>%  dplyr::rename(ReactomeID=X1, Name=X2, Species=X3)

reactome.pathway.result.nodes.pvalue <- set10.Sphere_CLS.spia.result.Reactom.filter %>% dplyr::select(Name, pGFdr, pSize, NDE,Status)

reactome.pathway.result.relationship.nodes <- left_join(reactome.pathway.result.relationship.nodes, reactome.pathway.result.nodes.pvalue)

reactome.pathway.result.manual.nodes.pvalue <- reactome.pathway.result.relationship.nodes %>% filter(is.na(NDE))
reactome.pathway.result.relationship.nodes <- reactome.pathway.result.relationship.nodes %>% filter(!is.na(NDE))

reactome.pathway.result.manual.nodes.pvalue <- reactome.pathway.result.manual.nodes.pvalue %>% dplyr::select(ReactomeID, Name, Species)
manual.add.nodes.info <- data.frame(ReactomeID=c("R-HSA-74160","R-HSA-166166", "R-HSA-937061", "R-HSA-69202", "R-HSA-109688", "R-HSA-975110", "R-HSA-937072"),
                                    pGFdr=c(8.761657e-23, 6.277138e-05, 6.277138e-05, 0.03652519,0.0108201, 4.555800e-05, 4.555800e-05),
                                    pSize=c(1505, 97, 97, 65, 61, 72, 72),
                                    NDE=c(513, 45, 45, 25, 26, 36, 36))
reactome.pathway.result.manual.nodes.pvalue <- left_join(reactome.pathway.result.manual.nodes.pvalue, manual.add.nodes.info)
reactome.pathway.result.relationship.nodes <- bind_rows(reactome.pathway.result.relationship.nodes, reactome.pathway.result.manual.nodes.pvalue)



write_delim(reactome.pathway.result.relationship, path="/Users/Weitinglin/Documents/R_scripts/Lab/Cytoscape/reactome_spia_result_visualization/reactome_spia_set10_2000_relation.txt", delim = "\t")
write_delim(reactome.pathway.result.relationship.nodes, path="/Users/Weitinglin/Documents/R_scripts/Lab/Cytoscape/reactome_spia_result_visualization/reactome_spia_set10_2000_node.txt", delim = "\t")



# 03_SPIAresult_summary ---------------------------------------------------

# SPIA with permutation 2000 times ----------------------------------------

set10.Sphere_CLS.spia.2000.result.Reactom        <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereCLS_2000_Reactome.csv", delim=",")         %>% filter(pGFdr < 0.05)  
set10.Sphere_Timeseries.spia.2000.result.Reactom <- read_delim(file="/Users/Weitinglin/Documents/Repository/code_in_lab/set10_005_spia_SphereTimeseries_2000_Reactome.csv", delim=",")  %>% filter(pGFdr < 0.05)  





# read the reactome database schema in order to know the hierarchy --------

esemble2reactome_all_levels <- read_delim(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/Ensembl2Reactome_All_Levels.txt", col_names = FALSE, delim="\t")
esembl2reactome <- read_delim(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/Ensembl2Reactome.txt", col_names = FALSE, delim="\t")
ensemble2reactomereactoms <- read_delim(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/Ensembl2ReactomeReactions.txt", col_names = FALSE, delim="\t")


reactome.pathway.map          <- read_delim(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/ReactomePathways.txt", col_names = FALSE, delim="\t") %>% filter(X3 %in% "Homo sapiens")
reactome.pathway.relationship <- read_delim(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/ReactomePathwaysRelation.txt", col_names = FALSE, delim="\t")


# find the reactome id mapping
set10.Sphere_CLS.spia.2000.pathway.id        <-  reactome.pathway.map %>%  filter(X2 %in% set10.Sphere_CLS.spia.2000.result.Reactom$Name )

set10.Sphere_Timeseries.spia.2000.pathway.id <-   reactome.pathway.map %>%  filter(X2 %in% set10.Sphere_Timeseries.spia.2000.result.Reactom$Name)




# _Sphere_CLS_result ------------------------------------------------------



set10.Sphere_CLS.spia.2000.pathway.relationship <- reactome.pathway.relationship %>% filter(X1 %in% set10.Sphere_CLS.spia.2000.pathway.id$X1 | X2 %in% set10.Sphere_CLS.spia.2000.pathway.id$X1)

#output interaction
write_delim(set10.Sphere_CLS.spia.2000.pathway.relationship, path="/Users/Weitinglin/Documents/R_scripts/Lab/Cytoscape/reactome_spia_result_visualization/reactome_spia_set10_2000_SphereCLS_relation.txt", delim = "\t")


set10.Sphere_CLS.spia.2000.pathway.related.node <- reactome.pathway.map %>% filter(X1 %in% unique(c(set10.Sphere_CLS.spia.2000.pathway.relationship$X1, set10.Sphere_CLS.spia.2000.pathway.relationship$X2)))
set10.Sphere_CLS.spia.2000.pathway.related.node <- set10.Sphere_CLS.spia.2000.pathway.related.node %>% dplyr::rename(ReactomeID=X1, PathwayName=X2, Species=X3)
set10.Sphere_CLS.spia.2000.result.Reactom.tmp <- set10.Sphere_CLS.spia.2000.result.Reactom %>% dplyr::select(Name, pSize, NDE, pNDE,pGFdr,Status) %>% dplyr::rename(PathwayName=Name)
set10.Sphere_CLS.spia.2000.pathway.related.node <- left_join(set10.Sphere_CLS.spia.2000.pathway.related.node, set10.Sphere_CLS.spia.2000.result.Reactom.tmp, by="PathwayName")
set10.Sphere_CLS.spia.2000.pathway.related.node$pGFdr[is.na(set10.Sphere_CLS.spia.2000.pathway.related.node$pGFdr)] <- 1

write_delim(set10.Sphere_CLS.spia.2000.pathway.related.node, path="/Users/Weitinglin/Documents/R_scripts/Lab/Cytoscape/reactome_spia_result_visualization/reactome_spia_set10_2000_SphereCLS_node.txt", delim = "\t")





# Input GSEA --------------------------------------------------------------
library(GSEABase)
ReactomePathway.geneset <- getGmt("/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/ReactomePathways.gmt")

ReactomePathway.title  <- ReactomePathway.geneset %>% map(setName)
ReactomePathway.Id     <- ReactomePathway.geneset %>% map(description)
ReactomePathway.symbol <- ReactomePathway.geneset %>% map(geneIds)
