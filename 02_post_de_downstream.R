# management the probe in different parameter setting statistics


library(readr)

# =*********Load the data*************= -------------------------------------------
# Wilcox_adj_nofilter ---------
QNlogwilcox.CLF_CLS    <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/CLF_CLS_wilcox_adj_nofilter.csv")
QNlogwilcox.CLF_Sphere <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/CLF_Sphere_wilcox_adj_nofilter.csv")
QNlogwilcox.CLS_Sphere <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/CLS_Sphere_wilcox_adj_nofilter.csv")
# Wilcox_no log transformation -------
QNnologwilcox.CLF_CLS    <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/nologQNwilcoxCLF_CLS.txt", 
                                       skip = 1,
                                       col_names = c("Num","gene","statistic","p.value","adjusted.p","adjusted.method"),
                                       delim = "\t")
QNnologwilcox.CLF_Sphere <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/nologQNwilcoxCLF_Sphere.txt", 
                                       skip = 1,
                                       col_names = c("Num","gene","statistic","p.value","adjusted.p","adjusted.method"),
                                       , delim = "\t")
QNnologwilcox.CLS_Sphere <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/nologQNwilcoxCLS_Sphere.txt",, 
                                       skip = 1,
                                       col_names = c("Num","gene","statistic","p.value","adjusted.p","adjusted.method"),
                                       delim = "\t")

# QN, Log2, t. test , adjust --------------------------------------------------------------------

tQNlogCLF_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/QNttestCLF_CLS.txt", delim = "\t")
tQNlogCLF_Sphere <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/QNttestCLF_Sphere.txt", delim = "\t")
tQNlogCLS_Sphere <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/QNttestCLS_Sphere.txt", delim = "\t")

# QN, Log2, filter variance 50% , t.test , adjust ------------------------------------------------
tQNfilter50logCLF_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/filter50QNttestCLF_CLS.txt", delim = "\t")
tQNfilter50logCLF_Sphere <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/filter50QNttestCLF_Sphere.txt", delim = "\t")
tQNfilter50logCLS_Sphere <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/filter50QNttestCLS_Sphere.txt", delim = "\t")



library(hgu133plus2.db)
library(annotate)

QNlogwilcox.CLF_CLS.probe<-QNlogwilcox.CLF_CLS$gene
QNlogwilcox.CLF_Sphere.probe<-QNlogwilcox.CLF_Sphere$gene
QNlogwilcox.CLS_Sphere.probe<-QNlogwilcox.CLS_Sphere$gene


QNnologwilcox.CLF_CLS.probe<-QNnologwilcox.CLF_CLS$gene
QNnologwilcox.CLF_Sphere.probe<-QNnologwilcox.CLF_Sphere$gene
QNnologwilcox.CLS_Sphere.probe<-QNnologwilcox.CLS_Sphere$gene


tQNlogCLF_CLS.probe<-tQNlogCLF_CLS$gene
tQNlogCLF_Sphere.probe<-tQNlogCLF_Sphere$gene
tQNlogCLS_Sphere.probe<-tQNlogCLS_Sphere$gene


tQNfilter50logCLF_CLS.probe<-tQNfilter50logCLF_CLS$gene
tQNfilter50logCLF_Sphere.probe<-tQNfilter50logCLF_Sphere$gene
tQNfilter50logCLS_Sphere.probe<-tQNfilter50logCLS_Sphere$gene

#debug
CLF_CLS <- tQNfilter50logCLF_CLS    
CLF_Sphere <- tQNfilter50logCLF_Sphere 
CLS_Sphere <- tQNfilter50logCLS_Sphere


take_intersect <- function(CLF_CLS, CLF_Sphere, CLS_Sphere, method="adjusted.p" ,p=0.05){
  if ( method == "adjusted.p"){
  CLF_CLS.probe    <- (CLF_CLS %>% filter(adjusted.p < p))$gene
  CLF_Sphere.probe <- (CLF_Sphere %>% filter(adjusted.p < p))$gene
  CLS_Sphere.probe <- (CLS_Sphere %>% filter(adjusted.p < p))$gene}
  else if ( method == "p.value"){
  CLF_CLS.probe    <- (CLF_CLS %>% filter(p.value < p))$gene
  CLF_Sphere.probe <- (CLF_Sphere %>% filter(p.value < p))$gene
  CLS_Sphere.probe <- (CLS_Sphere %>% filter(p.value < p))$gene}
  
  target.list <- setdiff(intersect(CLF_CLS.probe, CLS_Sphere.probe), CLF_Sphere.probe )
  return(target.list)
}

QNlogwilcox.target.probe       <- take_intersect(QNlogwilcox.CLF_CLS,
                                                 QNlogwilcox.CLF_Sphere,
                                                 QNlogwilcox.CLS_Sphere,
                                                 method = "p.value",
                                                 p = 0.05)
QNnologlogwilcox.target.probe  <- take_intersect(QNnologwilcox.CLF_CLS ,
                                                 QNnologwilcox.CLF_Sphere,
                                                 QNnologwilcox.CLS_Sphere,
                                                 method = "p.value",
                                                 p = 0.05)
QNtlog.target.probe_05      <- take_intersect(tQNlogCLF_CLS,
                                                 tQNlogCLF_Sphere,
                                                 tQNlogCLS_Sphere,
                                                 method = "adjusted.p",
                                                 p = 0.05)
QNtlog.target.probe_01      <- take_intersect(tQNlogCLF_CLS,
                                                    tQNlogCLF_Sphere,
                                                    tQNlogCLS_Sphere,
                                                    method = "adjusted.p",
                                                    p = 0.01)
QNtlog.target.probe_001      <- take_intersect(tQNlogCLF_CLS,
                                                    tQNlogCLF_Sphere,
                                                    tQNlogCLS_Sphere,
                                                    method = "adjusted.p",
                                                    p = 0.001)
QNtfilterlog.target.probe_05 <- take_intersect(tQNfilter50logCLF_CLS,
                                                 tQNfilter50logCLF_Sphere,
                                                 tQNfilter50logCLS_Sphere,
                                                 method = "adjusted.p",
                                                 p = 0.05)
QNtfilterlog.target.probe_01<- take_intersect(tQNfilter50logCLF_CLS,
                                                 tQNfilter50logCLF_Sphere,
                                                 tQNfilter50logCLS_Sphere,
                                                 method = "adjusted.p",
                                                 p = 0.01)
QNtfilterlog.target.probe_001<- take_intersect(tQNfilter50logCLF_CLS,
                                                 tQNfilter50logCLF_Sphere,
                                                 tQNfilter50logCLS_Sphere,
                                                 method = "adjusted.p",
                                                 p = 0.001)

# =***********************************= -------------------------------------------


# Find out two ------------------------------------------------------------


#為何會有這個奇怪的現象？
setdiff(QNtfilterlogwilcox.target.probe_01,QNtfilterlogwilcox.target.probe_001)
setdiff(QNtfilterlogwilcox.target.probe_001,QNtfilterlogwilcox.target.probe_01) ?
#[1] "1555788_a_at" "202912_at"   

which(tQNfilter50logCLF_CLS$gene == "1555788_a_at")
tQNfilter50logCLF_CLS[1351,]$adjusted.p 
#[1] 0.0005417983
which(tQNfilter50logCLF_Sphere$gene == "1555788_a_at")
tQNfilter50logCLF_Sphere[1293,]$adjusted.p
#0.001282036
which(tQNfilter50logCLS_Sphere$gene == "1555788_a_at")
#1508
tQNfilter50logCLS_Sphere[1508,]$adjusted.p
#[1] 0.0005491053


# make the upset plot -----------------------------------------------------


library(UpSetR)
listInput <- list(wilcox = QNlogwilcox.target.probe,
                  wilcoxnolog = QNnologlogwilcox.target.probe,
                  t.test_05 = QNtlogwilcox.target.probe_05,
                  t.test_01 =QNtlogwilcox.target.probe_01,
                  t.test_001 = QNtlogwilcox.target.probe_001,
                  t.test_filterby50_05 = QNtfilterlogwilcox.target.probe_05,
                  t.test_filterby50_01 = QNtfilterlogwilcox.target.probe_01,
                  t.test_filterby50_001 = QNtfilterlogwilcox.target.probe_001)

upset(fromList(listInput), nsets = 8, order.by = "freq")


# annotation information from the internet --------------------------------
install.packages("readxl")
library(readxl)
cancer_stem_marker <- read_excel("/Users/Weitinglin/Documents/R_scripts/Lab/CSCdbdata/marker_list.xls")
#57
cancer_stem_related_marker <- read_excel("/Users/Weitinglin/Documents/R_scripts/Lab/CSCdbdata/annotation.xls")
#13,308
study.interest <- c("NANOG","IGF1R","SDF1","POU5F2","POU5F1P3","POU5F1P4","POU5F1B","HGF","POU5F1","IGF2","LIF","TGF","SOX2","THY1","JAK1",
                                 "CXCR4","IGF2","THY-1","ABCG2","ABCG","ABCG1","ALDH1","CD133","SDF2",
                                 "CD44","PROM1","MYC","BMI1","EZH2","KIT","PAR","IGFBP","IGFBP1","IGFBP2","GCSF",
                                 "CHBD","ANGPT","NT","IL3","IL4","IL5","IL7","CCL18","CCL11","FGF2","EPHA1","TGFBR1","TGFBR2","LIFR","SMAD2","YAP1","TCF21")
#52

# Specify the Gene Synbol -------------------------------------------------
probe.relatedTF.data.path <- file.path("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/proberelateTF.txt")
table.probe.relate.TF <- read_delim(probe.relatedTF.data.path, delim = " ", col_types = c("ciiii"))
unname(getSYMBOL(table.probe.relate.TF$Probe, "hgu133plus2"))
table.probe.relate.TF <- table.probe.relate.TF %>% mutate(Symbol = unname(getSYMBOL(table.probe.relate.TF$Probe, "hgu133plus2")))
#3949

sum(table.probe.relate.TF$Symbol   %in% cancer_stem_marker$symbol)
#label the cancer stem marker
table.probe.relate.TF$Cancer_Stem_Marker[table.probe.relate.TF$Symbol   %in% cancer_stem_marker$symbol] <- 1
table.probe.relate.TF$Cancer_Stem_Marker[!table.probe.relate.TF$Symbol   %in% cancer_stem_marker$symbol] <- 0
#50

sum(table.probe.relate.TF$Symbol   %in% cancer_stem_related_marker$symbol)
#885
#label the cancer related 
table.probe.relate.TF$Cancer_Stem_related[table.probe.relate.TF$Symbol   %in% cancer_stem_related_marker$symbol] <- 1
table.probe.relate.TF$Cancer_Stem_related[!table.probe.relate.TF$Symbol   %in% cancer_stem_related_marker$symbol] <- 0


# t.test_filter_05 =====================================================

t.test_filter_05 <- getSYMBOL(QNtfilterlog.target.probe_05,"hgu133plus2")
t.test_filter_05 <- data.frame(probe = names(t.test_filter_05), symbol = unname(t.test_filter_05))



#tagged the TF
t.test_filter_05$TF[t.test_filter_05$symbol %in% table.probe.relate.TF$Symbol] <- 1
t.test_filter_05$TF[!t.test_filter_05$symbol %in% table.probe.relate.TF$Symbol] <- 0
#cancer stem marker
t.test_filter_05$CS[t.test_filter_05$symbol %in% cancer_stem_marker$symbol] <- 1
t.test_filter_05$CS[!t.test_filter_05$symbol %in% cancer_stem_marker$symbol] <- 0
#cancer stem related
t.test_filter_05$CSRelated[t.test_filter_05$symbol  %in% cancer_stem_related_marker$symbol] <- 1
t.test_filter_05$CSRelated[!t.test_filter_05$symbol  %in% cancer_stem_related_marker$symbol] <- 0
# t.test_filter_01 =====================================================
t.test_filter_01 <- getSYMBOL(QNtfilterlog.target.probe_01,"hgu133plus2")
t.test_filter_01 <- data.frame(probe = names(t.test_filter_01), symbol = unname(t.test_filter_01))
#tagged the TF
t.test_filter_01$TF[t.test_filter_01$symbol %in% table.probe.relate.TF$Symbol] <- 1
t.test_filter_01$TF[!t.test_filter_01$symbol %in% table.probe.relate.TF$Symbol] <- 0
#cancer stem marker
t.test_filter_01$CS[t.test_filter_01$symbol %in% cancer_stem_marker$symbol] <- 1
t.test_filter_01$CS[!t.test_filter_01$symbol %in% cancer_stem_marker$symbol] <- 0
#cancer stem related
t.test_filter_01$CSRelated[t.test_filter_01$symbol  %in% cancer_stem_related_marker$symbol] <- 1
t.test_filter_01$CSRelated[!t.test_filter_01$symbol  %in% cancer_stem_related_marker$symbol] <- 0
# t.test_filter_001 =====================================================
t.test_filter_001 <- getSYMBOL(QNtfilterlog.target.probe_001, "hgu133plus2")
t.test_filter_001 <- data.frame(probe = names(t.test_filter_001), symbol = unname(t.test_filter_001))
#tagged the TF
t.test_filter_001$TF[t.test_filter_001$symbol %in% table.probe.relate.TF$Symbol] <- 1
t.test_filter_001$TF[!t.test_filter_001$symbol %in% table.probe.relate.TF$Symbol] <- 0
#cancer stem marker
t.test_filter_001$CS[t.test_filter_001$symbol %in% cancer_stem_marker$symbol] <- 1
t.test_filter_001$CS[!t.test_filter_001$symbol %in% cancer_stem_marker$symbol] <- 0
#cancer stem related
t.test_filter_001$CSRelated[t.test_filter_001$symbol  %in% cancer_stem_related_marker$symbol] <- 1
t.test_filter_001$CSRelated[!t.test_filter_001$symbol  %in% cancer_stem_related_marker$symbol] <- 0
# t.test_filter three group venndiagram -----------------------------------


CLF_CLS.probe    <- t.test_filter_05$symbol
CLF_sphere.probe <- t.test_filter_01$symbol
CLS_sphere.probe <- t.test_filter_001$symbol
area.1 <- length(CLF_CLS.probe)
area.2 <- length(CLF_sphere.probe)
area.3 <- length(CLS_sphere.probe)
area.12 <- length(intersect(CLF_CLS.probe, CLF_sphere.probe))
area.23 <- length(intersect(CLF_sphere.probe,CLS_sphere.probe))
area.13 <- length(intersect(CLF_CLS.probe, CLS_sphere.probe))
area.123 <- length(intersect(intersect(CLF_CLS.probe, CLF_sphere.probe),CLS_sphere.probe))
area.list <- list(area.1, area.2, area.3, area.12, area.23, area.13, area.123)

library(VennDiagram)
draw.triple.venn(area1 = area.list[[1]], area2 = area.list[[2]], area3 = area.list[[3]], n12 = area.list[[4]], n23 = area.list[[5]], n13 = area.list[[6]],
                 n123 = area.list[[7]], category = c("p =0.05", "p = 0.01", "p = 0.001"), lty = "blank",
                 fill = c("skyblue", "pink1", "mediumorchid"))


# t.test__05 =====================================================

t.test__05 <- getSYMBOL(QNtlog.target.probe_05,"hgu133plus2")
t.test__05 <- data.frame(probe = names(t.test__05), symbol = unname(t.test__05))



#tagged the TF
t.test__05$TF[t.test__05$symbol %in% table.probe.relate.TF$Symbol] <- 1
t.test__05$TF[!t.test__05$symbol %in% table.probe.relate.TF$Symbol] <- 0
#cancer stem marker
t.test__05$CS[t.test__05$symbol %in% cancer_stem_marker$symbol] <- 1
t.test__05$CS[!t.test__05$symbol %in% cancer_stem_marker$symbol] <- 0
#cancer stem related
t.test__05$CSRelated[t.test__05$symbol  %in% cancer_stem_related_marker$symbol] <- 1
t.test__05$CSRelated[!t.test__05$symbol  %in% cancer_stem_related_marker$symbol] <- 0
# t.test__01 =====================================================
t.test__01 <- getSYMBOL(QNtlog.target.probe_01,"hgu133plus2")
t.test__01 <- data.frame(probe = names(t.test__01), symbol = unname(t.test__01))
#tagged the TF
t.test__01$TF[t.test__01$symbol %in% table.probe.relate.TF$Symbol] <- 1
t.test__01$TF[!t.test__01$symbol %in% table.probe.relate.TF$Symbol] <- 0
#cancer stem marker
t.test__01$CS[t.test__01$symbol %in% cancer_stem_marker$symbol] <- 1
t.test__01$CS[!t.test__01$symbol %in% cancer_stem_marker$symbol] <- 0
#cancer stem related
t.test__01$CSRelated[t.test__01$symbol  %in% cancer_stem_related_marker$symbol] <- 1
t.test__01$CSRelated[!t.test__01$symbol  %in% cancer_stem_related_marker$symbol] <- 0


# t.test__001 =====================================================
t.test__001 <- getSYMBOL(QNtlog.target.probe_001, "hgu133plus2")
t.test__001 <- data.frame(probe = names(t.test__001), symbol = unname(t.test__001))
#tagged the TF
t.test__001$TF[t.test__001$symbol %in% table.probe.relate.TF$Symbol] <- 1
t.test__001$TF[!t.test__001$symbol %in% table.probe.relate.TF$Symbol] <- 0
#cancer stem marker
t.test__001$CS[t.test__001$symbol %in% cancer_stem_marker$symbol] <- 1
t.test__001$CS[!t.test__001$symbol %in% cancer_stem_marker$symbol] <- 0
#cancer stem related
t.test__001$CSRelated[t.test__001$symbol  %in% cancer_stem_related_marker$symbol] <- 1
t.test__001$CSRelated[!t.test__001$symbol  %in% cancer_stem_related_marker$symbol] <- 0


# t.test_ three group venndiagram -----------------------------------


CLF_CLS.probe    <- t.test__05$symbol
CLF_sphere.probe <- t.test__01$symbol
CLS_sphere.probe <- t.test__001$symbol
area.1 <- length(CLF_CLS.probe)
area.2 <- length(CLF_sphere.probe)
area.3 <- length(CLS_sphere.probe)
area.12 <- length(intersect(CLF_CLS.probe, CLF_sphere.probe))
area.23 <- length(intersect(CLF_sphere.probe,CLS_sphere.probe))
area.13 <- length(intersect(CLF_CLS.probe, CLS_sphere.probe))
area.123 <- length(intersect(intersect(CLF_CLS.probe, CLF_sphere.probe),CLS_sphere.probe))
area.list <- list(area.1, area.2, area.3, area.12, area.23, area.13, area.123)

library(VennDiagram)
draw.triple.venn(area1 = area.list[[1]], area2 = area.list[[2]], area3 = area.list[[3]], n12 = area.list[[4]], n23 = area.list[[5]], n13 = area.list[[6]],
                 n123 = area.list[[7]], category = c("p =0.05", "p = 0.01", "p = 0.001"), lty = "blank",
                 fill = c("skyblue", "pink1", "mediumorchid"))



# save the result ---------------------------------------------------------
path <- "/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/result"
setwd(path)

t.test_filter_01 <- t.test_filter_01 %>%
  mutate(total = TF + CS + CSRelated) %>%
  arrange(desc(total), desc(TF), desc(CS)) %>% mutate(p_01 = 1)
  
t.test_filter_05 <- t.test_filter_05 %>%
  mutate(total = TF + CS + CSRelated) %>%
  arrange(desc(total), desc(TF), desc(CS)) %>% mutate(p_05 = 1)
  
t.test_filter_001 <- t.test_filter_001 %>%
  mutate(total = TF + CS + CSRelated) %>%
  arrange(desc(total), desc(TF), desc(CS)) %>% mutate(p_001 = 1)

library(dplyr)
tmp <- full_join(t.test_filter_01, t.test_filter_05, by = c("probe", "symbol", "TF", "CS", "CSRelated","total")) 
tmp <- full_join(tmp,t.test_filter_001 , by =  c("probe", "symbol", "TF", "CS", "CSRelated","total"))
tmp$p_01[is.na(tmp$p_01)] <- 0
tmp$p_05[is.na(tmp$p_05)] <- 0
tmp$p_001[is.na(tmp$p_001)] <- 0
tmp <- tmp %>% mutate(all = p_01 + p_05 + p_001) %>% arrange(desc(total), desc(all))
tmp$interst[tmp$symbol %in% study.interest] <- 1
tmp$interst[!tmp$symbol %in% study.interest] <- 0
tmp <- tmp %>% arrange(desc(interst)) 
write.csv(tmp, "t.test_filter_total.csv")

write.csv(t.test_filter_05,"t.test_filter_05.csv")
write.csv(t.test_filter_01, "t.test_filter_01.csv")
write.csv(t.test_filter_001,"t.test_filter_001.csv" )

t.test__01 <- t.test__01 %>%
  mutate(total = TF + CS + CSRelated) %>%
  arrange(desc(total), desc(TF), desc(CS)) %>% mutate(p_01 = 1)
t.test__05 <- t.test__05 %>%
  mutate(total = TF + CS + CSRelated) %>%
  arrange(desc(total), desc(TF), desc(CS))  %>% mutate(p_05 = 1)
t.test__001 <-t.test__001 %>%
  mutate(total = TF + CS + CSRelated) %>%
  arrange(desc(total), desc(TF), desc(CS))  %>% mutate(p_001 = 1)

tmp <- full_join(t.test__01, t.test__05, by = c("probe", "symbol", "TF", "CS", "CSRelated","total")) 
tmp <- full_join(tmp,t.test__001 , by =  c("probe", "symbol", "TF", "CS", "CSRelated","total"))
tmp$p_01[is.na(tmp$p_01)] <- 0
tmp$p_05[is.na(tmp$p_05)] <- 0
tmp$p_001[is.na(tmp$p_001)] <- 0
tmp <- tmp %>% mutate(all = p_01 + p_05 + p_001) %>% arrange(desc(total), desc(all))
tmp$interst[tmp$symbol %in% study.interest] <- 1
tmp$interst[!tmp$symbol %in% study.interest] <- 0
tmp <- tmp %>% arrange(desc(interst)) 
write.csv(tmp, "t.test__total.csv")

write.csv(t.test__05, "t.test__05.csv")
write.csv(t.test__01, "t.test__01.csv")
write.csv(t.test__001, "t.test__001.csv")


tmp %>% filter(p_05 ==1) %>% dim()
tmp %>% filter(p_05 ==1 & CS == 1) %>% dim()
tmp %>% filter(p_05 ==1 & TF == 1) %>% dim()
tmp %>% filter(p_05 ==1 & CSRelated == 1) %>% dim()
tmp %>% filter(p_01 ==1 & CS == 1) %>% dim()
tmp %>% filter(p_01 ==1 & TF == 1) %>% dim()
tmp %>% filter(p_01 ==1 & CSRelated == 1) %>% dim()
tmp %>% filter(p_01 ==1 ) %>% dim()
tmp %>% filter(p_001 ==1 ) %>% dim()
tmp %>% filter(p_001 ==1 & CS == 1) %>% dim()
tmp %>% filter(p_001 ==1 & TF == 1) %>% dim()
tmp %>% filter(p_001 ==1 & CSRelated == 1) %>% dim()


# *********************************** -------------------------------------


# Heatmap -----------------------------------------------------------------


# _load the raw data ------------------------------------------------------
library(readr)
library(dplyr)
norm.exprs.path <- file.path("/Users/Weitinglin/Documents/R_scripts/Lab/microarry20161107andshiny/norm.txt")
norm <- read_delim(norm.exprs.path, delim = "\t", col_names = TRUE)




# _t.test_05 --------------------------------------------------------------
t.test_05 <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/result/t.test__05.csv")
t.test_05.probe <- t.test_05$probe


# _t.test_01 --------------------------------------------------------------
t.test_01 <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/result/t.test__01.csv")

# _t.test_001 -------------------------------------------------------------
t.test_001 <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/result/t.test__001.csv")

# _t.test_filter_05 -------------------------------------------------------
t.test_filter_05 <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/result/t.test_filter_05.csv")

# _t.test_filter_01 -------------------------------------------------------
t.test_filter_01 <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/result/t.test_filter_01.csv")

# _t.test_filter_001 ------------------------------------------------------
t.test_filter_001 <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/result/t.test_filter_001.csv")


# 3DscatterPlot -----------------------------------------------------------

# ＿arrange the data -------------------------------------------------------
#combine the set into mean



# _t.test_05 --------------------------------------------------------------

DE.annotate <- read_csv("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/result/t.test__total.csv")

TF.list <- (t.test_05 %>% filter(TF==1))$probe
nonTF.list <- (t.test_05 %>% filter(!TF==1))$probe


scatter <- norm %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                           CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                           Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) %>%
  dplyr::select(Probe, CLF, CLS, Sphere)

scatter$tag[scatter$Probe %in% TF.list] <- "TF"
scatter$tag[scatter$Probe %in% nonTF.list] <- "nonTF"
scatter$tag[is.na(scatter$tag)] <- "notDE"
library(hgu133plus2.db)
library(annotate)
scatter$Symbol <- getSYMBOL(scatter$Probe, "hgu133plus2")


library(RColorBrewer)
usr.col <- brewer.pal(3, "Set1")
library(plotly)
p <- plot_ly(scatter, x = ~CLF, y = ~CLS, z = ~Sphere,
             color = ~tag, colors = usr.col, text= ~paste('Symbol:',Symbol,'<br>Probe:',Probe)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'CLF'),
                      yaxis = list(title = 'CLS'),
                      zaxis = list(title = 'Sphere')))




# _t.test_01 --------------------------------------------------------------

TF.list <- (t.test_01 %>% filter(TF==1))$probe
nonTF.list <- (t.test_01 %>% filter(!TF==1))$probe

norm <- norm %>% rename(Probe = `“Probe”`)
scatter <- norm %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                           CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                           Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) %>%
  dplyr::select(Probe, CLF, CLS, Sphere)
scatter$tag <- NA
scatter$tag[scatter$Probe %in% TF.list] <- "TF"
scatter$tag[scatter$Probe %in% nonTF.list] <- "nonTF"
scatter$tag[is.na(scatter$tag)] <- "notDE"
library(hgu133plus2.db)
library(annotate)
scatter$Symbol <- getSYMBOL(scatter$Probe, "hgu133plus2")


library(RColorBrewer)
usr.col <- brewer.pal(3, "Set1")
library(plotly)
p <- plot_ly(scatter, x = ~CLF, y = ~CLS, z = ~Sphere,
             color = ~tag, colors = usr.col, text= ~paste('Symbol:',Symbol,'<br>Probe:',Probe)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'CLF'),
                      yaxis = list(title = 'CLS'),
                      zaxis = list(title = 'Sphere')))


# _t.test_001 -------------------------------------------------------------



TF.list <- (t.test_001 %>% filter(TF==1))$probe
nonTF.list <- (t.test_001 %>% filter(!TF==1))$probe

norm <- norm %>% rename(Probe = `“Probe”`)
scatter <- norm %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                           CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                           Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) %>%
  dplyr::select(Probe, CLF, CLS, Sphere)
scatter$tag <- NA
scatter$tag[scatter$Probe %in% TF.list] <- "TF"
scatter$tag[scatter$Probe %in% nonTF.list] <- "nonTF"
scatter$tag[is.na(scatter$tag)] <- "notDE"
library(hgu133plus2.db)
library(annotate)
scatter$Symbol <- getSYMBOL(scatter$Probe, "hgu133plus2")


library(RColorBrewer)
usr.col <- brewer.pal(3, "Set1")
library(plotly)
p <- plot_ly(scatter, x = ~CLF, y = ~CLS, z = ~Sphere,
             color = ~tag, colors = usr.col, text= ~paste('Symbol:',Symbol,'<br>Probe:',Probe)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'CLF'),
                      yaxis = list(title = 'CLS'),
                      zaxis = list(title = 'Sphere')))



# _t.test_filter_05 -------------------------------------------------------
TF.list <- (t.test_filter_05 %>% filter(TF==1))$probe
nonTF.list <- (t.test_filter_05 %>% filter(!TF==1))$probe

norm <- norm %>% rename(Probe = `“Probe”`)
scatter <- norm %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                           CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                           Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) %>%
  dplyr::select(Probe, CLF, CLS, Sphere)
scatter$tag <- NA
scatter$tag[scatter$Probe %in% TF.list] <- "TF"
scatter$tag[scatter$Probe %in% nonTF.list] <- "nonTF"
scatter$tag[is.na(scatter$tag)] <- "notDE"
library(hgu133plus2.db)
library(annotate)
scatter$Symbol <- getSYMBOL(scatter$Probe, "hgu133plus2")


library(RColorBrewer)
usr.col <- brewer.pal(3, "Set1")
library(plotly)
p <- plot_ly(scatter, x = ~CLF, y = ~CLS, z = ~Sphere,
             color = ~tag, colors = usr.col, text= ~paste('Symbol:',Symbol,'<br>Probe:',Probe)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'CLF'),
                      yaxis = list(title = 'CLS'),
                      zaxis = list(title = 'Sphere')))

# _t.test_filter_01 -------------------------------------------------------
TF.list <- (t.test_filter_01 %>% filter(TF==1))$probe
nonTF.list <- (t.test_filter_01 %>% filter(!TF==1))$probe

norm <- norm %>% rename(Probe = `“Probe”`)
scatter <- norm %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                           CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                           Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) %>%
  dplyr::select(Probe, CLF, CLS, Sphere)
scatter$tag <- NA
scatter$tag[scatter$Probe %in% TF.list] <- "TF"
scatter$tag[scatter$Probe %in% nonTF.list] <- "nonTF"
scatter$tag[is.na(scatter$tag)] <- "notDE"
library(hgu133plus2.db)
library(annotate)
scatter$Symbol <- getSYMBOL(scatter$Probe, "hgu133plus2")


library(RColorBrewer)
usr.col <- brewer.pal(3, "Set1")
library(plotly)
p <- plot_ly(scatter, x = ~CLF, y = ~CLS, z = ~Sphere,
             color = ~tag, colors = usr.col, text= ~paste('Symbol:',Symbol,'<br>Probe:',Probe)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'CLF'),
                      yaxis = list(title = 'CLS'),
                      zaxis = list(title = 'Sphere')))

# _t.test_filter_001 ------------------------------------------------------

TF.list <- (t.test_filter_001 %>% filter(TF==1))$probe
nonTF.list <- (t.test_filter_001 %>% filter(!TF==1))$probe

norm <- norm %>% rename(Probe = `“Probe”`)
scatter <- norm %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                           CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                           Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) %>%
  dplyr::select(Probe, CLF, CLS, Sphere)
scatter$tag <- NA
scatter$tag[scatter$Probe %in% TF.list] <- "TF"
scatter$tag[scatter$Probe %in% nonTF.list] <- "nonTF"
scatter$tag[is.na(scatter$tag)] <- "notDE"
library(hgu133plus2.db)
library(annotate)
scatter$Symbol <- getSYMBOL(scatter$Probe, "hgu133plus2")


library(RColorBrewer)
usr.col <- brewer.pal(3, "Set1")
library(plotly)
p <- plot_ly(scatter, x = ~CLF, y = ~CLS, z = ~Sphere,
             color = ~tag, colors = usr.col, text= ~paste('Symbol:',Symbol,'<br>Probe:',Probe)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'CLF'),
                      yaxis = list(title = 'CLS'),
                      zaxis = list(title = 'Sphere')))


# *********************************** -------------------------------------

# check Cancer Stem Marker ------------------------------------------------

#use the scatter from the previous step
SYMBOL <- getSYMBOL(scatter$Probe,"hgu133plus2")

scatter <- scatter %>% mutate(Symbol = SYMBOL)
scatter <- scatter[,c(1,5,3,2,4)]
csc.markder <- c("SOX2","OCT4","NANOG")

#CAF crosstalk with Lung CSCs


#search sting pattern with stringr
index <- which(str_detect(SYMBOL,"POU5"))
index <- which(str_detect(SYMBOL,"KLF4"))
index <- which(str_detect(SYMBOL,"NOTCH1"))
index <- which(str_detect(SYMBOL,"ABCG2"))
index <- which(str_detect(SYMBOL,"IGF1R"))
index <- c(index, which(str_detect(SYMBOL,"SOX2")))
index <- c(index, which(str_detect(SYMBOL,"NANOG")))
#basic standardization method

#test
temp <- scatter[index,c(3,2,4,1,5)]
temp
filter <- which(!(str_detect(temp$Probe,c("_s_") | str_detect(temp$Probe, c("_x_"))))) 
temp <- temp[filter,]
temp <- temp %>% mutate(Name=paste(Probe, Symbol, sep=":")) %>% tidyr::gather("Type","Value",1:3)
temp$Type <- factor(temp$Type, levels=c("CLS", "CLF", "Sphere"))
ggplot(data= temp, aes(x=Type, y= Name))+geom_tile(aes(fill=Value), colour="black", size=2)+scale_fill_gradient( low = "#FF3300" , high = "#006600")





# *********************************** -------------------------------------

# QN, Log2, t. test , adjust --------------------------------------------------------------------

tQNlogCLF_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/QNttestCLF_CLS.txt", delim = "\t")
tQNlogCLF_Sphere <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/QNttestCLF_Sphere.txt", delim = "\t")
tQNlogCLS_Sphere <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/QNttestCLS_Sphere.txt", delim = "\t")

# QN, Log2, filter variance 50% , t.test , adjust ------------------------------------------------
tQNfilter50logCLF_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/filter50QNttestCLF_CLS.txt", delim = "\t")
tQNfilter50logCLF_Sphere <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/filter50QNttestCLF_Sphere.txt", delim = "\t")
tQNfilter50logCLS_Sphere <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray20161111/filter50QNttestCLS_Sphere.txt", delim = "\t")




# use gprofiler to query the NA value
#total.probe.data.frame <- convertNA_probe(total.probe.data.frame)




gconvert("228498_at",
         organism = "hsapiens",
         target = "AFFY_HG_U133_PLUS_2",
         df = T,
         numeric_ns = "AFFY_HG_U133_PLUS_2",
         mthreshold = 3)

# function depend on gprofiler for transforming the tag
convertNA_probe <- function(x){
    n <- dim(x)[[1]]
    
    for ( i in 1:n){
        if (is.na(x[[i,6]])){
            tmp <-   gconvert(x[[i,1]],
                              organism = "hsapiens",
                              target = "AFFY_HG_U133_PLUS_2",
                              df = T,
                              numeric_ns = "AFFY_HG_U133_PLUS_2",
                              mthreshold = 1)
            
            if(is.null(tmp$name)){
                x[[i,6]] <- NA
            }else{
                x[[i,6]] <- tmp$name[[1]]
            }
        }
    }
    return(x)
}

# function to annotated the t.test result
annotate_the_ttest <- function(x){
    x <- x %>% mutate(Symbol = getSYMBOL(Probe,"hgu133plus2"))
    
    
}



# ******************************** ----------------------------------------


# Data input --------------------------------------------------------------


# No filter ---------------------------------------------------------------
#for CLS-CLF: No filter
# _greater

t.greater.CLF_CLS      <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/greaterQNttestCLF_CLS.txt", delim = "\t")
t.greater.Sphere_CLF   <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/greaterQNttestCLF_Sphere.txt", delim = "\t")
t.greater.Sphere_CLS   <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/greaterQNttestCLS_Sphere.txt", delim = "\t")

# _less 
t.less.CLF_CLS      <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/lessQNttestCLF_CLS.txt", delim="\t")
t.less.Sphere_CLF   <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/lessQNttestCLF_Sphere.txt", delim="\t")
t.less.Sphere_CLS   <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/lessQNttestCLS_Sphere.txt", delim="\t")
# UP-regulation 


t.greater.CLF_CLS <- t.greater.CLF_CLS %>% mutate(Case = "t.greater.CLF_CLS", Method = "Nofilter")
t.greater.Sphere_CLF <- t.greater.Sphere_CLF %>% mutate(Case = "t.greater.Sphere_CLF", Method = "Nofilter")
t.greater.Sphere_CLS <- t.greater.Sphere_CLS %>% mutate(Case = "t.greater.Sphere_CLS", Method = "Nofilter")
total.t.greater <- bind_rows(t.greater.CLF_CLS, t.greater.Sphere_CLF, t.greater.Sphere_CLS)


# DOWN-regulation 
t.less.CLF_CLS <- t.less.CLF_CLS %>% mutate(Case = "t.less.CLF_CLS", Method = "Nofilter")
t.less.Sphere_CLF <- t.less.Sphere_CLF %>% mutate(Case = "t.less.Sphere_CLF", Method = "Nofilter")
t.less.Sphere_CLS <- t.less.Sphere_CLS %>% mutate(Case = "t.less.Sphere_CLS", Method = "Nofilter")
total.t.less <- bind_rows(t.less.CLF_CLS, t.less.Sphere_CLF, t.less.Sphere_CLS)

# group_without_filter 
total_without_filter <- bind_rows(total.t.greater, total.t.less)


# Filter ------------------------------------------------------------------

#for CLS-CLF: filter by 0.5
# _greater 
f.t.greater.CLF_CLS      <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/greaterfilter50QNttestCLF_CLS.txt", delim="\t")
f.t.greater.Sphere_CLF   <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/greaterfilter50QNttestCLF_Sphere.txt", delim="\t")
f.t.greater.Sphere_CLS   <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/greaterfilter50QNttestCLS_Sphere.txt", delim="\t")

# _less 
f.t.less.CLF_CLS      <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/lessfilter50QNttestCLF_CLS.txt", delim="\t")
f.t.less.Sphere_CLF   <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/lessfilter50QNttestCLF_Sphere.txt", delim="\t")
f.t.less.Sphere_CLS   <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/lessfilter50QNttestCLS_Sphere.txt", delim="\t")




# t.test with filter 


# UP-regulation 
f.t.greater.CLF_CLS <- f.t.greater.CLF_CLS %>% mutate(Case = "f.t.greater.CLF_CLS", Method = "Filter")
f.t.greater.Sphere_CLF <- f.t.greater.Sphere_CLF %>% mutate(Case = "f.t.greater.Sphere_CLF", Method = "Filter")
f.t.greater.Sphere_CLS <- f.t.greater.Sphere_CLS %>% mutate(Case = "f.t.greater.Sphere_CLS", Method = "Filter")
total.f.t.greater <- bind_rows(f.t.greater.CLF_CLS, f.t.greater.Sphere_CLF, f.t.greater.Sphere_CLS)


# DOWN-regulation 
f.t.less.CLF_CLS <- f.t.less.CLF_CLS %>% mutate(Case = "f.t.less.CLF_CLS", Method = "Filter")
f.t.less.Sphere_CLF <- f.t.less.Sphere_CLF %>% mutate(Case = "f.t.less.Sphere_CLF", Method = "Filter") 
f.t.less.Sphere_CLS <- f.t.less.Sphere_CLS %>% mutate(Case = "f.t.less.Sphere_CLS", Method = "Filter") 
total.f.t.less <- bind_rows(f.t.less.CLF_CLS, f.t.less.Sphere_CLF, f.t.less.Sphere_CLS)



# group_with_filter 
total_with_filter    <- bind_rows(total.f.t.greater, total.f.t.less)



# combine all 
total_ttest_result <- bind_rows(total_without_filter, total_with_filter)
save(total_ttest_result, file = "total_ttest_result.Rdata")

total_ttest_result %>% mutate(A = 0.5*(estimate1 + estimate2)) %>% group_by(Case) %>% 
    summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
              A_0.05 = mean(A[adjusted.p < 0.05], na.rm = TRUE),
              min_0.05 = min(A[adjusted.p < 0.05], na.rm = TRUE),
              max_0.05 = max(A[adjusted.p < 0.05], na.rm = TRUE),
              n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
              A_0.01 = mean(A[adjusted.p < 0.01]),
              min_0.01 = min(A[adjusted.p < 0.05], na.rm = TRUE),
              max_0.01 = max(A[adjusted.p < 0.05], na.rm = TRUE),
              n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
              A_0.001 = mean(A[adjusted.p < 0.001]),
              min_0.001 = min(A[adjusted.p < 0.05], na.rm = TRUE),
              max_0.001 = max(A[adjusted.p < 0.05], na.rm = TRUE),
              n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE),
              A_0.0001 = mean(A[adjusted.p < 0.0001], na.rm = TRUE),
              min_0.0001 = min(A[adjusted.p < 0.05], na.rm = TRUE),
              max_0.0001 = max(A[adjusted.p < 0.05], na.rm = TRUE))

total_ttest_result %>%mutate(A = 0.5*(estimate1 + estimate2)) %>% ggplot() + geom_boxplot(aes(x = Case, y = A)) + facet_grid(Method ~ .)
total_ttest_result %>%mutate(A = 0.5*(estimate1 + estimate2)) %>% ggplot() + geom_violin(aes(x = Case, y = A)) + facet_grid(Method ~ .) +
    geom_hline(yintercept = )

# QN, Log2, withoutfilter upregulation and downregulation ----------------




# t-test without filter
# up-regulation
total.t.greater.table <-  total.t.greater %>% group_by(Case) %>% summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                                                                           n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                                                                           n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                                                                           n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE))
total.t.greater <- total.t.greater %>% dplyr::rename(Probe = gene)
total.t.greater.annotate <- left_join(total.t.greater, total.probe.data.frame, by = "Probe")
total.t.greater.annotate %>% group_by(Case) %>% summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                                                          n_0.05_TF = sum(adjusted.p < 0.05 & !is.na(Total), na.rm = TRUE),
                                                          n_0.05_CS = sum(adjusted.p < 0.05 & CSmarker == 1, na.rm = TRUE),
                                                          n_0.05_CSrelated = sum(adjusted.p < 0.05 & CSrelated == 1, na.rm = TRUE),
                                                          n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                                                          n_0.01_TF = sum(adjusted.p < 0.01 & !is.na(Total), na.rm = TRUE),
                                                          n_0.01_CS = sum(adjusted.p < 0.01 & CSmarker == 1, na.rm = TRUE),
                                                          n_0.01_CSrelated = sum(adjusted.p < 0.01 & CSrelated == 1, na.rm = TRUE),
                                                          n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                                                          n_0.001_TF = sum(adjusted.p < 0.001 & !is.na(Total),na.rm = TRUE),
                                                          n_0.001_CS = sum(adjusted.p < 0.001 & CSmarker == 1, na.rm = TRUE),
                                                          n_0.001_CSrelated = sum(adjusted.p < 0.001 & CSrelated == 1, na.rm = TRUE),
                                                          n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE),
                                                          n_0.0001_TF = sum(adjusted.p < 0.0001 & !is.na(Total), na.rm = TRUE),
                                                          n_0.0001_CS = sum(adjusted.p < 0.0001 & CSmarker ==1, na.rm = TRUE),
                                                          n_0.0001_CSrelated = sum(adjusted.p < 0.0001 & CSrelated ==1, na.rm = TRUE)) %>%
    tidyr::unite("n_0.0001",c(14:17),sep = "/") %>%
    tidyr::unite("n_0.001", c(10:13), sep = "/") %>%
    tidyr::unite("n_0.01", c(6:9), sep = "/") %>%
    tidyr::unite("n_0.05", c(2:5), sep = "/") 




# down-regulation
total.t.less.table <- total.t.less %>% group_by(Case) %>% summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                                                                    n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                                                                    n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                                                                    n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE))

total.t.less <- total.t.less %>% dplyr::rename(Probe = gene)
total.t.less.annotate <- left_join(total.t.less, total.probe.data.frame, by = "Probe")
total.t.less.annotate %>% group_by(Case) %>% summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                                                       n_0.05_TF = sum(adjusted.p < 0.05 & !is.na(Total), na.rm = TRUE),
                                                       n_0.05_CS = sum(adjusted.p < 0.05 & CSmarker == 1, na.rm = TRUE),
                                                       n_0.05_CSrelated = sum(adjusted.p < 0.05 & CSrelated == 1, na.rm = TRUE),
                                                       n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                                                       n_0.01_TF = sum(adjusted.p < 0.01 & !is.na(Total), na.rm = TRUE),
                                                       n_0.01_CS = sum(adjusted.p < 0.01 & CSmarker == 1, na.rm = TRUE),
                                                       n_0.01_CSrelated = sum(adjusted.p < 0.01 & CSrelated == 1, na.rm = TRUE),
                                                       n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                                                       n_0.001_TF = sum(adjusted.p < 0.001 & !is.na(Total),na.rm = TRUE),
                                                       n_0.001_CS = sum(adjusted.p < 0.001 & CSmarker == 1, na.rm = TRUE),
                                                       n_0.001_CSrelated = sum(adjusted.p < 0.001 & CSrelated == 1, na.rm = TRUE),
                                                       n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE),
                                                       n_0.0001_TF = sum(adjusted.p < 0.0001 & !is.na(Total), na.rm = TRUE),
                                                       n_0.0001_CS = sum(adjusted.p < 0.0001 & CSmarker ==1, na.rm = TRUE),
                                                       n_0.0001_CSrelated = sum(adjusted.p < 0.0001 & CSrelated ==1, na.rm = TRUE)) %>%
    tidyr::unite("n_0.0001",c(14:17),sep = "/") %>%
    tidyr::unite("n_0.001", c(10:13), sep = "/") %>%
    tidyr::unite("n_0.01", c(6:9), sep = "/") %>%
    tidyr::unite("n_0.05", c(2:5), sep = "/") 

# gprofiler result load
# output the genesybol for the gprofiler 
subset <- total.t.less %>% filter(Case == "t.less.CLF_CLS", adjusted.p < 0.001 )   %>% mutate(symbol = getSYMBOL(gene,"hgu133plus2.db"))
test <- subset %>% filter(is.na(symbol)) %>% mutate(new = convertNA_probe(gene)) 
test$new[test$new == "1"] <- NA
subset$symbol[is.na(subset$symbol)] <- test$new
sum(is.na(subset$symbol))

subset <- total.t.less %>% filter(Case == "t.less.Sphere_CLF", adjusted.p < 0.001)  %>% mutate(symbol = getSYMBOL(gene,"hgu133plus2.db"))
test <- subset %>% filter(is.na(symbol)) %>% mutate(new = convertNA_probe(gene))
test$new[test$new == "1"] <- NA
subset$symbol[is.na(subset$symbol)] <- test$new
sum(is.na(subset$symbol))
subset$symbol[!is.na(subset$symbol)] %>% cat

subset <- total.t.less %>% filter(Case == "t.less.Sphere_CLS", adjusted.p < 0.003)  %>% mutate(symbol = getSYMBOL(gene,"hgu133plus2.db"))
test <- subset %>% filter(is.na(symbol)) %>% mutate(new = convertNA_probe(gene))
test$new[test$new == "1"] <- NA
subset$symbol[is.na(subset$symbol)] <- test$new
sum(is.na(subset$symbol))
subset$symbol[!is.na(subset$symbol)] %>% cat
system("say finished")




# Without Filter 
gprofilter.greater.CLF_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/gprofiler/ttestnofilter10_3hierchyfilter/gprofiler_results_upregulation_CLF_CLS.txt", delim = "\t")
gprofilter.greater.Sphere_CLF <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/gprofiler/ttestnofilter10_3hierchyfilter/gprofiler_results_upregulation_Sphere_CLF.txt", delim = "\t")
gprofilter.greater.Sphere_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/gprofiler/ttestnofilter10_3hierchyfilter/gprofiler_results_upregulation_Sphere_CLS.txt", delim = "\t")

gprofilter.less.CLF_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/gprofiler/ttestWITHfilter10_3hierchyfilter /gprofiler_results_filter_correct_downregulation_CLF_CLS.txt", delim = "\t")
gprofilter.less.Sphere_CLF <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/gprofiler/ttestnofilter10_3hierchyfilter/gprofiler_results_downregulation_Sphere_CLF.txt", delim = "\t")
gprofilter.less.Sphere_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/gprofiler/ttestWITHfilter10_3hierchyfilter /gprofiler_results_filter_downregulation_correct_Sphere_CLS.txt", delim = "\t")


gprofilter.greater.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "keg") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "BP")
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "CC")
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "MF")
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "keg")
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "hp")
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "omi")
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "rea")
gprofilter.greater.CLF_CLS$`t type` %>% table()


gprofilter.greater.Sphere_CLF    %>% filter(`t type` == "BP")
gprofilter.greater.Sphere_CLF  %>% filter(`t type` == "CC")
gprofilter.greater.Sphere_CLF  %>% filter(`t type` == "MF")
gprofilter.greater.Sphere_CLF     %>% filter(`t type` == "keg")
gprofilter.greater.Sphere_CLF     %>% filter(`t type` == "hp")
gprofilter.greater.Sphere_CLF     %>% filter(`t type` == "omi")
gprofilter.greater.Sphere_CLF     %>% filter(`t type` == "rea")
gprofilter.greater.Sphere_CLF$`t type` %>% table()

gprofilter.greater.Sphere_CLS   %>% filter(`t type` == "BP")
gprofilter.greater.Sphere_CLS %>% filter(`t type` == "CC")
gprofilter.greater.Sphere_CLS %>% filter(`t type` == "MF")
gprofilter.greater.Sphere_CLS     %>% filter(`t type` == "keg")
gprofilter.greater.Sphere_CLS     %>% filter(`t type` == "hp")
gprofilter.greater.Sphere_CLS     %>% filter(`t type` == "omi")
gprofilter.greater.Sphere_CLS     %>% filter(`t type` == "rea")
gprofilter.greater.Sphere_CLS$`t type` %>% table()

gprofilter.less.CLF_CLS      %>% filter(`t type` == "BP")
gprofilter.less.CLF_CLS       %>% filter(`t type` == "CC")
gprofilter.less.CLF_CLS       %>% filter(`t type` == "MF")
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "keg")
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "hp")
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "omi")
gprofilter.greater.CLF_CLS     %>% filter(`t type` == "rea")
gprofilter.greater.CLF_CLS$`t type` %>% table()

gprofilter.less.Sphere_CLF    %>% filter(`t type` == "BP")
gprofilter.less.Sphere_CLF    %>% filter(`t type` == "CC")
gprofilter.less.Sphere_CLF    %>% filter(`t type` == "MF")
gprofilter.greater.Sphere_CLF     %>% filter(`t type` == "keg")
gprofilter.greater.Sphere_CLF     %>% filter(`t type` == "hp")
gprofilter.greater.Sphere_CLF     %>% filter(`t type` == "omi")
gprofilter.greater.Sphere_CLF     %>% filter(`t type` == "rea")
gprofilter.greater.Sphere_CLF$`t type` %>% table()

gprofilter.less.Sphere_CLS    %>% filter(`t type` == "BP")
gprofilter.less.Sphere_CLS    %>% filter(`t type` == "CC")
gprofilter.less.Sphere_CLS     %>% filter(`t type` == "MF")
gprofilter.greater.Sphere_CLS     %>% filter(`t type` == "keg")
gprofilter.greater.Sphere_CLS     %>% filter(`t type` == "hp")
gprofilter.greater.Sphere_CLS     %>% filter(`t type` == "omi")
gprofilter.greater.Sphere_CLS     %>% filter(`t type` == "rea")
gprofilter.greater.Sphere_CLS$`t type` %>% table()



down.CLF_CLS.nofilter    <- (total.t.less.annotate %>% filter(Case == "t.less.CLF_CLS" & adjusted.p < 0.001  ))$Probe %>% unique
down.Sphere_CLS.nofilter <- (total.t.less.annotate %>% filter(Case == "t.less.Sphere_CLS" & adjusted.p < 0.001 ))$Probe %>% unique

area.1 <- length(down.CLF_CLS )
area.2 <- length(down.Sphere_CLS)
area.12 <- length(intersect(down.CLF_CLS, down.Sphere_CLS))
area.list <- list(area.1, area.2, area.12)


#  *************************************************------








# Summary the test result 

# up-regulation
total.f.t.greater.table <-  total.f.t.greater %>% group_by(Case) %>% summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                                                                               n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                                                                               n_0.005 = sum(adjusted.p < 0.005, na.rm = TRUE),
                                                                               n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                                                                               n_0.0005 = sum(adjusted.p < 0.0005, na.rm = TRUE),
                                                                               n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE))
total.f.t.greater <- total.f.t.greater %>% dplyr::rename(Probe = gene)
total.f.t.greater.annotate <- left_join(total.f.t.greater, total.probe.data.frame, by = "Probe")
total.f.t.greater.annotate %>% group_by(Case) %>% summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                                                       n_0.05_TF = sum(adjusted.p < 0.05 & !is.na(Total), na.rm = TRUE),
                                                       n_0.05_CS = sum(adjusted.p < 0.05 & CSmarker == 1, na.rm = TRUE),
                                                       n_0.05_CSrelated = sum(adjusted.p < 0.05 & CSrelated == 1, na.rm = TRUE),
                                                       n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                                                       n_0.01_TF = sum(adjusted.p < 0.01 & !is.na(Total), na.rm = TRUE),
                                                       n_0.01_CS = sum(adjusted.p < 0.01 & CSmarker == 1, na.rm = TRUE),
                                                       n_0.01_CSrelated = sum(adjusted.p < 0.01 & CSrelated == 1, na.rm = TRUE),
                                                       n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                                                       n_0.001_TF = sum(adjusted.p < 0.001 & !is.na(Total),na.rm = TRUE),
                                                       n_0.001_CS = sum(adjusted.p < 0.001 & CSmarker == 1, na.rm = TRUE),
                                                       n_0.001_CSrelated = sum(adjusted.p < 0.001 & CSrelated == 1, na.rm = TRUE),
                                                       n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE),
                                                       n_0.0001_TF = sum(adjusted.p < 0.0001 & !is.na(Total), na.rm = TRUE),
                                                       n_0.0001_CS = sum(adjusted.p < 0.0001 & CSmarker ==1, na.rm = TRUE),
                                                       n_0.0001_CSrelated = sum(adjusted.p < 0.0001 & CSrelated ==1, na.rm = TRUE)) %>%
  tidyr::unite("n_0.0001",c(14:17),sep = "/") %>%
  tidyr::unite("n_0.001", c(10:13), sep = "/") %>%
  tidyr::unite("n_0.01", c(6:9), sep = "/") %>%
  tidyr::unite("n_0.05", c(2:5), sep = "/") 


# down-regulation
total.f.t.less.table <- total.f.t.less %>% group_by(Case) %>% summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                                                                        n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                                                                        n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                                                                        n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE))
total.f.t.less <- total.f.t.less %>% dplyr::rename(Probe = gene)
total.f.t.less.annotate <- left_join(total.f.t.less, total.probe.data.frame, by = "Probe")
total.f.t.less.annotate %>% group_by(Case) %>% summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                                                       n_0.05_TF = sum(adjusted.p < 0.05 & !is.na(Total), na.rm = TRUE),
                                                       n_0.05_CS = sum(adjusted.p < 0.05 & CSmarker == 1, na.rm = TRUE),
                                                       n_0.05_CSrelated = sum(adjusted.p < 0.05 & CSrelated == 1, na.rm = TRUE),
                                                       n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                                                       n_0.01_TF = sum(adjusted.p < 0.01 & !is.na(Total), na.rm = TRUE),
                                                       n_0.01_CS = sum(adjusted.p < 0.01 & CSmarker == 1, na.rm = TRUE),
                                                       n_0.01_CSrelated = sum(adjusted.p < 0.01 & CSrelated == 1, na.rm = TRUE),
                                                       n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                                                       n_0.001_TF = sum(adjusted.p < 0.001 & !is.na(Total),na.rm = TRUE),
                                                       n_0.001_CS = sum(adjusted.p < 0.001 & CSmarker == 1, na.rm = TRUE),
                                                       n_0.001_CSrelated = sum(adjusted.p < 0.001 & CSrelated == 1, na.rm = TRUE),
                                                       n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE),
                                                       n_0.0001_TF = sum(adjusted.p < 0.0001 & !is.na(Total), na.rm = TRUE),
                                                       n_0.0001_CS = sum(adjusted.p < 0.0001 & CSmarker ==1, na.rm = TRUE),
                                                       n_0.0001_CSrelated = sum(adjusted.p < 0.0001 & CSrelated ==1, na.rm = TRUE)) %>%
  tidyr::unite("n_0.0001",c(14:17),sep = "/") %>%
  tidyr::unite("n_0.001", c(10:13), sep = "/") %>%
  tidyr::unite("n_0.01", c(6:9), sep = "/") %>%
  tidyr::unite("n_0.05", c(2:5), sep = "/") 







# t-test with filter
# up-regulation


# output the genesybol for the gprofiler 
subset <- total.f.t.greater %>% filter(Case == "f.t.greater.CLF_CLS", adjusted.p < 0.001 )   %>% mutate(symbol = getSYMBOL(Probe,"hgu133plus2.db"))
test <- subset %>% filter(is.na(symbol)) %>% mutate(new = convertNA_probe(gene))
test$new[test$new == "1"] <- NA
subset$symbol[is.na(subset$symbol)] <- test$new
sum(is.na(subset$symbol))
subset$symbol[!is.na(subset$symbol)] %>% cat
system("say finished")


subset <- total.f.t.greater %>% filter(Case == "f.t.greater.Sphere_CLF", adjusted.p < 0.001 )   %>% mutate(symbol = getSYMBOL(gene,"hgu133plus2.db"))
test <- subset %>% filter(is.na(symbol)) %>% mutate(new = convertNA_probe(gene))
test$new[test$new == "1"] <- NA
subset$symbol[is.na(subset$symbol)] <- test$new
sum(is.na(subset$symbol))
subset$symbol[!is.na(subset$symbol)] %>% cat
system("say finished")


subset <- total.f.t.greater %>% filter(Case == "f.t.greater.Sphere_CLS", adjusted.p < 0.001 )   %>% mutate(symbol = getSYMBOL(gene,"hgu133plus2.db"))
test <- subset %>% filter(is.na(symbol)) %>% mutate(new = convertNA_probe(gene))
test$new[test$new == "1"] <- NA
subset$symbol[is.na(subset$symbol)] <- test$new
sum(is.na(subset$symbol))
subset$symbol[!is.na(subset$symbol)] %>% cat
system("say finished")

subset <- total.f.t.greater %>% filter(Case == "f.t.greater.CLF_CLS", adjusted.p < 0.001 ) 



# down-regulation
library(gProfileR)
hgu133plus2.probe.list <- rownames(ecelfile.set)
probe <- subset$Probe
result<- gprofiler(query = probe,
          organism = "hsapiens",
          ordered_query = T,
          max_set_size = 1000,
          min_isect_size = 0,
          significant = T,
          exclude_iea = T,
          underrep = T,
          correction_method = "fdr",
          hier_filtering = "moderate",
          custom_bg = hgu133plus2.probe)


# output the genesybol for the gprofiler 
subset <- total.f.t.less %>% filter(Case == "f.t.less.CLF_CLS", adjusted.p < 0.001 )   %>% mutate(symbol = getSYMBOL(Probe,"hgu133plus2.db"))
test <- subset %>% filter(is.na(symbol)) %>% mutate(new = convertNA_probe(gene))
test$new[test$new == "1"] <- NA
subset$symbol[is.na(subset$symbol)] <- test$new
sum(is.na(subset$symbol))
subset$symbol[!is.na(subset$symbol)] %>% cat
subset$Probe %>% cat
system("say finished")


subset <- total.f.t.less %>% filter(Case == "f.t.less.Sphere_CLF", adjusted.p < 0.001 )   %>% mutate(symbol = getSYMBOL(gene,"hgu133plus2.db"))
test <- subset %>% filter(is.na(symbol)) %>% mutate(new = convertNA_probe(gene))
test$new[test$new == "1"] <- NA
subset$symbol[is.na(subset$symbol)] <- test$new
sum(is.na(subset$symbol))
subset$symbol[!is.na(subset$symbol)] %>% cat
system("say finished")


subset <- total.f.t.less %>% filter(Case == "f.t.less.Sphere_CLS", adjusted.p < 0.001 )   %>% mutate(symbol = getSYMBOL(Probe,"hgu133plus2.db"))

test <- subset %>% filter(is.na(symbol)) %>% mutate(new = convertNA_probe(gene))
test$new[test$new == "1"] <- NA
subset$symbol[is.na(subset$symbol)] <- test$new
sum(is.na(subset$symbol))
subset$symbol[!is.na(subset$symbol)] %>% unique%>% cat
system("say finished")




# Filter
gprofilterfilter.greater.CLF_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/gprofiler/ttestWITHfilter10_3hierchyfilter /gprofiler_results_filter_upregulation_CLF_CLS.txt", delim = "\t")
gprofilterfilter.greater.Sphere_CLF <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/gprofiler/ttestWITHfilter10_3hierchyfilter /gprofiler_results_filter_upregulation_Sphere_CLF.txt", delim = "\t")
gprofilterfilter.greater.Sphere_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/gprofiler/ttestWITHfilter10_3hierchyfilter /gprofiler_results_filter_003_greater_Sphere_CLS.txt", delim = "\t")
gprofilterfilter.less.CLF_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/lessfilter50QNttestCLF_CLS.txt", delim = "\t")
gprofilterfilter.less.Sphere_CLF <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/lessfilter50QNttestCLF_Sphere.txt", delim = "\t")
gprofilterfilter.less.Sphere_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/microarray/result/Differential_Expression/upanddownseperatemeasure/lessfilter50QNttestCLS_Sphere.txt", delim = "\t")

gprofilterfilter.less.CLF_CLS <- read_delim("/Users/Weitinglin/Documents/R_scripts/Lab/gprofiler/ttestfilter/gprofiler_results_201612filter_CLF_CLS_greater.txt", delim="\t")


gprofilterfilter.greater.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "BP") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "MF") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "CC") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "rea") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "keg") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)

gprofilterfilter.greater.Sphere_CLF[,-c(1,2)]   %>% filter(`t type` == "BP") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.Sphere_CLF[,-c(1,2)]   %>% filter(`t type` == "MF") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.Sphere_CLF[,-c(1,2)]   %>% filter(`t type` == "CC") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.Sphere_CLF[,-c(1,2)]   %>% filter(`t type` == "rea") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.Sphere_CLF[,-c(1,2)]   %>% filter(`t type` == "keg") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)

gprofilterfilter.greater.Sphere_CLS[,-c(1,2)]   %>% filter(`t type` == "BP") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.Sphere_CLS[,-c(1,2)]   %>% filter(`t type` == "MF") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.Sphere_CLS[,-c(1,2)]   %>% filter(`t type` == "CC") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.Sphere_CLS[,-c(1,2)]   %>% filter(`t type` == "rea") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.greater.Sphere_CLS[,-c(1,2)]   %>% filter(`t type` == "keg") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)

gprofilterfilter.less.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "BP") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "MF") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "CC") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "rea") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.CLF_CLS[,-c(1,2)]   %>% filter(`t type` == "keg") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)

gprofilterfilter.less.Sphere_CLF[,-c(1,2)]   %>% filter(`t type` == "BP") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.Sphere_CLF[,-c(1,2)]   %>% filter(`t type` == "MF") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.Sphere_CLF[,-c(1,2)]   %>% filter(`t type` == "CC") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.Sphere_CLF[,-c(1,2)]   %>% filter(`t type` == "rea") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.Sphere_CLF[,-c(1,2)]   %>% filter(`t type` == "keg") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)

gprofilterfilter.less.Sphere_CLS[,-c(1,2)]   %>% filter(`t type` == "BP") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.Sphere_CLS[,-c(1,2)]   %>% filter(`t type` == "MF") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.Sphere_CLS[,-c(1,2)]   %>% filter(`t type` == "CC") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.Sphere_CLS[,-c(1,2)]   %>% filter(`t type` == "rea") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)
gprofilterfilter.less.Sphere_CLS[,-c(1,2)]   %>% filter(`t type` == "keg") %>% mutate_each(funs(as.numeric),c(1:3)) %>% arrange(`p-value`)












# total merge dataframe
total_without_filter.table <- total_without_filter %>% group_by(Case) %>% summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                                                      n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                                                      n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                                                      n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE))


total_with_filter.table <- total_with_filter %>% group_by(Case) %>% summarise(n_0.05 = sum(adjusted.p < 0.05, na.rm = TRUE),
                                                   n_0.01 = sum(adjusted.p < 0.01, na.rm = TRUE),
                                                   n_0.001 = sum(adjusted.p < 0.001, na.rm = TRUE),
                                                   n_0.0001 = sum(adjusted.p < 0.0001, na.rm = TRUE))


total.f.t.greater.annotate %>% filter(Case == "f.t.greater.CLF_CLS" & !is.na(Total) & adjusted.p < 0.001 )
total.f.t.greater.annotate %>% filter(Case == "f.t.greater.CLF_CLS" & CSmarker == 1 & adjusted.p < 0.001 )
(total.f.t.greater.annotate %>% filter(Case == "f.t.greater.CLF_CLS" & CSrelated == 1 & adjusted.p < 0.001 ))$Symbol %>% table()
(total.f.t.less.annotate %>% filter(Case == "f.t.less.CLF_CLS" & !is.na(Total) &  adjusted.p < 0.001 ))$Symbol %>% table()
(total.f.t.less.annotate %>% filter(Case == "f.t.less.CLF_CLS" & CSmarker == 1 &  adjusted.p < 0.001 ))$Symbol %>% table()
(total.f.t.less.annotate %>% filter(Case == "f.t.less.CLF_CLS" & CSrelated == 1 &  adjusted.p < 0.001 ))$Symbol %>% table()

(total.f.t.greater.annotate %>% filter(Case == "f.t.greater.Sphere_CLS" & !is.na(Total) & adjusted.p < 0.001 ))$Symbol %>% table
(total.f.t.greater.annotate %>% filter(Case == "f.t.greater.Sphere_CLS" & CSmarker == 1 & adjusted.p < 0.001 ))$Symbol
(total.f.t.greater.annotate %>% filter(Case == "f.t.greater.CLF_CLS" & CSrelated == 1 & adjusted.p < 0.001 ))$Symbol %>% table()
(total.f.t.less.annotate %>% filter(Case == "f.t.less.Sphere_CLS" & !is.na(Total) &  adjusted.p < 0.001 ))$Symbol %>% table()
(total.f.t.less.annotate %>% filter(Case == "f.t.less.Sphere_CLS" & CSmarker == 1 &  adjusted.p < 0.001 ))$Symbol %>% table()
(total.f.t.less.annotate %>% filter(Case == "f.t.less.Sphere_CLS" & CSrelated == 1 &  adjusted.p < 0.001 ))$Symbol %>% table()


up.CLF_CLS    <- (total.f.t.greater.annotate %>% filter(Case == "f.t.greater.CLF_CLS" & adjusted.p < 0.001 ))$Probe %>% unique
up.Sphere_CLS <- (total.f.t.greater.annotate %>% filter(Case == "f.t.greater.Sphere_CLS" & adjusted.p < 0.001 ))$Probe %>% unique
area.1 <- length(up.CLF_CLS )
area.2 <- length(up.Sphere_CLS)
area.12 <- length(intersect(up.CLF_CLS, up.Sphere_CLS))
area.list <- list(area.1, area.2, area.12)

#target.probe.up.set
#SEC61A1
#CLF_CLS:"217716_s_at"
#Sphere_CLS:"222385_x_at"
#PHLDA3
#CLF_CLS:"218634_at","235534_at"
#Sphere_CLS:"218634_at"
#SMIM7
#CLF_CLS:"221988_at"
#Sphere_CLS:"229650_s_at"
target.probe.up.set    <- c("217716_s_at","222385_x_at","218634_at","235534_at","218634_at","221988_at","229650_s_at")
target.probe.down.set  <- c("1554510_s_at", "1554747_a_at", "1565717_s_at", "201931_at","208622_s_at","209774_x_at",
                            "212514_x_at" ,"213988_s_at","217234_s_at","218172_s_at","228498_at","222987_s_at","222988_s_at",
                            "226065_at","226069_at","226065_at","226069_at","230708_at")

time.expression <- norm
time.expression <- ecelfile.set
time.expression <- as.data.frame(time.expression)
time.expression <- time.expression[ row.names(time.expression) %in% target.probe.up.set,] 
time.expression <- t(time.expression) %>% scale %>% t()
attributes(time.expression)$`scaled:center` <- NULL
attributes(time.expression)$`scaled:scale` <- NULL
time.expression <- as.data.frame(time.expression)
time.expression$probe <- row.names(time.expression)
time.expression <- time.expression %>% tidyr::gather("time","value",-probe)
time.expression$time <- factor(time.expression$time, levels = c("p3", "p6", "p8", "p14")) 
ggplot(data=time.expression ) + geom_point(aes(x=time, y=value)) + facet_grid(~probe) 
ggplot(data=time.expression) + geom_point(aes(x=time, y=value)) + facet_grid(~probe) 
# ggplot(data=time.expression %>% filter(probe %in% target.probe.up.set)) + geom_point(aes(x=time, y=value)) + facet_grid(~probe) 
# ggplot(data=time.expression %>% filter(probe %in% target.probe.down.set)) + geom_point(aes(x=time, y=value)) + facet_grid(~probe) 
#target.probe.down.set 
#same probe
#1554510_s_at" "1554747_a_at" "1565717_s_at" "201931_at"    "208622_s_at"  "209774_x_at"  
#"212514_x_at"  "213988_s_at"  "217234_s_at"  "218172_s_at" "228498_at" **


#same gene
#"TMEM9"
#CLF_CLS:"222987_s_at"
#Sphere_CLS:"222988_s_at"
#"PRICKLE1"
#CLF_CLS:"226065_at","226069_at"
#Sphere_CLS:"230708_at"


down.CLF_CLS    <- (total.f.t.less.annotate %>% filter(Case == "f.t.less.CLF_CLS" & adjusted.p < 0.001  ))$Probe %>% unique
down.Sphere_CLS <- (total.f.t.less.annotate %>% filter(Case == "f.t.less.Sphere_CLS" & adjusted.p < 0.001 ))$Probe %>% unique

area.1 <- length(down.CLF_CLS )
area.2 <- length(down.Sphere_CLS)
area.12 <- length(intersect(down.CLF_CLS, down.Sphere_CLS))
area.list <- list(area.1, area.2, area.12)

library(VennDiagram)
draw.pairwise.venn(area1 = area.list[[1]], area2 = area.list[[2]], n12 = area.list[[3]], c("CLF_CLS", "Sphere_CLS"), lty = "blank")


# QN, Log2, Annova, Tukey -------------------------------------------------

aov.Sphere_CLS <- test.tukey.final.result %>% filter(case ==  "Sphere-CLS")
aov.CLS_CLF    <- test.tukey.final.result %>% filter(case ==   "CLS-CLF")

# filter with annova p value < 0.05
aov.p.Sphere_CLS  <- aov.CLS_CLF %>% filter(annova.p.value < 0.05)
aov.p.CLS_CLF     <- aov.CLS_CLF %>% filter(annova.p.value < 0.05)

# adjusted the p.value with the BH

test.tukey.final.result %>% filter(case != "Sphere-CLF") %>%
                            group_by(case) %>%
                            mutate(BH.p.adj=p.adjust(p.adj, method="BH")) %>%
                            summarise(
                                      N_0.05=sum(BH.p.adj < 0.05),
                                      N_0.01=sum(BH.p.adj < 0.01),
                                      N_0.001=sum(BH.p.adj < 0.001),
                                      N_0.0001=sum(BH.p.adj < 0.0001),
                                      N_0.00001=sum(BH.p.adj < 0.00001))

#分出up,down的case
# Up-regulation(只得是CLF相對於CLS,或是Sphere相對CLS)
tukey.aov.adjust.method <- "BH"
up.tukey.aov.CLF_CLS        <- test.tukey.final.result %>%
                                filter(case == "CLS-CLF", diff < 0) %>% 
                                mutate(BH.p.adj=p.adjust(p.adj, method=tukey.aov.adjust.method))
up.tukey.aov.CLF_CLS %>% summarise(n_005=sum(BH.p.adj < 0.05), n_001=sum(BH.p.adj < 0.01), n_0001=sum(BH.p.adj < 0.001), n_00001=sum(BH.p.adj < 0.0001))                          
  
up.tukey.aov.Sphere_CLS     <- test.tukey.final.result %>%
                                filter(case == "Sphere-CLS", diff > 0) %>% 
                                mutate(BH.p.adj=p.adjust(p.adj, method=tukey.aov.adjust.method))
up.tukey.aov.Sphere_CLS  %>% summarise(n_005=sum(BH.p.adj < 0.05), n_001=sum(BH.p.adj < 0.01), n_0001=sum(BH.p.adj < 0.001), n_00001=sum(BH.p.adj < 0.0001))




# Down-regulation
down.tukey.aov.CLF_CLS      <- test.tukey.final.result %>%
                                filter(case == "CLS-CLF", diff > 0) %>% 
                                mutate(BH.p.adj=p.adjust(p.adj, method=tukey.aov.adjust.method))
down.tukey.aov.CLF_CLS %>% summarise(n_005=sum(BH.p.adj < 0.05), n_001=sum(BH.p.adj < 0.01), n_0001=sum(BH.p.adj < 0.001), n_00001=sum(BH.p.adj < 0.0001))

down.tukey.aov.Sphere_CLS   <- test.tukey.final.result %>%
                                filter(case == "Sphere-CLS", diff < 0) %>% 
                                mutate(BH.p.adj=p.adjust(p.adj, method=tukey.aov.adjust.method))
up.tukey.aov.Sphere_CLS %>% summarise(n_005=sum(BH.p.adj < 0.05), n_001=sum(BH.p.adj < 0.01), n_0001=sum(BH.p.adj < 0.001), n_00001=sum(BH.p.adj < 0.0001))

# Total
save(up.tukey.aov.CLF_CLS,up.tukey.aov.Sphere_CLS,down.tukey.aov.CLF_CLS,down.tukey.aov.Sphere_CLS,file = "updown_tukey_aov_result.Rdata")

tmp.up   <- bind_rows(up.tukey.aov.CLF_CLS, up.tukey.aov.Sphere_CLS) %>% mutate(type = "up")

tmp.down <- bind_rows(down.tukey.aov.CLF_CLS, down.tukey.aov.Sphere_CLS) %>% mutate(type = "down")

total.aov.tukey <- bind_rows(tmp.up, tmp.down) 


total.aov.tukey %>%
    group_by(case, type) %>%
    summarise(N_0.05=sum(BH.p.adj < 0.05),
              N_0.01=sum(BH.p.adj < 0.01),
              N_0.001=sum(BH.p.adj < 0.001),
              N_0.0001=sum(BH.p.adj < 0.0001),
              N_0.00001=sum(BH.p.adj < 0.00001))



MA.aov.data <- norm %>% as.data.frame %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                                CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                                Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) 
MA.aov.data$Probe <- rownames(norm)
MA.aov.data <- MA.aov.data %>% dplyr::select(CLF,CLS,Sphere,Probe)


#比較兩者交集的
#畫交集圖


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

# Upregulation
cutoff.p <- 0.00001
up.CLF_CLS    <- up.tukey.aov.CLF_CLS %>% filter(BH.p.adj < cutoff.p)
up.Sphere_CLS <- up.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)

area.1 <- length(up.CLF_CLS$Probe )
area.2 <- length(up.Sphere_CLS$Probe)
area.12 <- length(intersect(up.CLF_CLS$Probe, up.Sphere_CLS$Probe))
area.list <- list(area.1, area.2, area.12)
draw.pairwise.venn(area1 = area.list[[1]], area2 = area.list[[2]], cross.area = area.list[[3]],
                   category = c("P6+fibroblast > P6", "P6-Sphere > p6"),
                   fill = c("Blue","Red"))
# Downregulation
cutoff.p <- 0.00001
down.CLF_CLS    <- down.tukey.aov.CLF_CLS  %>% filter(BH.p.adj < cutoff.p)
down.Sphere_CLS <- down.tukey.aov.Sphere_CLS %>% filter(BH.p.adj < cutoff.p)
# The DE gene related to the 
area.1 <- length(down.CLF_CLS$Probe )
area.2 <- length(down.Sphere_CLS$Probe )
area.12 <- length(intersect(down.CLF_CLS$Probe, down.Sphere_CLS$Probe ))
area.list <- list(area.1, area.2, area.12)
draw.pairwise.venn(area1 = area.list[[1]], area2 = area.list[[2]], cross.area = area.list[[3]],
                   category = c("P6+fibroblast < P6", "P6-Sphere < p6"),
                   fill = c("Blue", "Red"))





# _down.stream analysis ---------------------------------------------------
# Hierachial cluster use heatmaply
union.up       <- union(up.CLF_CLS$Probe, up.Sphere_CLS$Probe)
intersect.up   <- intersect(up.CLF_CLS$Probe, up.Sphere_CLS$Probe)

union.down     <- union(down.CLF_CLS$Probe, down.Sphere_CLS$Probe)
intersect.down <- intersect(down.CLF_CLS$Probe, down.Sphere_CLS$Probe)

subset.norm <- norm[rownames(norm) %in% c(union.up),]
subset.norm <- norm[rownames(norm) %in% c(intersect.up),]
subset.norm <- norm[rownames(norm) %in% c(intersect.down),]
subset.norm <- norm[rownames(norm) %in% c(intersect.up, intersect.down),]
subset.norm <- norm[rownames(norm) %in% c(union.down),]
subset.norm <- norm[rownames(norm) %in% c(union.up,union.down),]
subset.norm <- t(scale(t(subset.norm)))
colnames(subset.norm) <- c("P6+Fib-1","P6+Fib-2","P6+Fib-3","P6-1","P6-2","P6-3","Sphere_p6-1","Sphere_p6-2","Sphere_p6-3")
RdYlGn <- colorRampPalette(brewer.pal(11, "RdYlGn"))
heatmaply(subset.norm, k_col = 2, k_row = 4,
          xaxis_font_size = "10pt", yaxis_font_size = "10pt",
          row_text_angle = -40, column_text_angle =  90,
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "green", high = "red", midpoint = 0) ) %>%
          layout(margin = list(l = 180, b = 100))

heatmaply(subset.norm, k_col = 2, k_row = 4,
          xaxis_font_size = "10pt", yaxis_font_size = "10pt",
          labRow = rep("",nrow(subset.norm)), column_text_angle =  90,
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "green", high = "red", midpoint = 0) ) %>%
    layout(margin = list(l = 180, b = 100))
library(heatmaply)
heatmaply(subset.norm, k_col = 2, k_row = 4,
          xaxis_font_size = "10pt", yaxis_font_size = "10pt",
          row_text_angle = -40, column_text_angle =  90,
          colors = rev(RdYlGn(256)),
          seriate = "GW" ) %>%
    layout(margin = list(l = 180, b = 100))


library(gplots)
heatmap.2(subset.norm, srtCol=30,labRow = FALSE, keysize = 1, offsetCol = 0, col=rev(RdYlGn(256)),
          trace = "none")




# scatter plot ------------------------------------------------------------
library(RColorBrewer)
usr.col <- brewer.pal(3, "Set1")

scatter.norm <- norm %>% as.data.frame %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                                              CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                                              Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) %>% dplyr::select(CLS, Sphere, CLF)
scatter.norm <- scatter.norm %>% mutate(Probe = probe.name)
scatter.norm <- scatter.norm %>% dplyr::select(Probe,CLF,CLS,Sphere)
scatter.norm <- left_join(scatter.norm, (annotation.probe.dataframe %>% dplyr::select(Probe,name,CSmarker,CSrelated,Total)), by ="Probe")

scatter.norm$tag[scatter.norm$Probe %in% union.up] <- "total upregulation"
scatter.norm$tag[scatter.norm$Probe %in% union.down] <- "total downregulation"
scatter.norm$tag[scatter.norm$Probe %in% intersect.up] <- "intersect of upregulation"
scatter.norm$tag[scatter.norm$Probe %in% intersect.down] <- "intersect of downregulation"
scatter.norm$tag[is.na(scatter.norm$tag)] <- "not differential expression"


plot_ly(scatter.norm, x = ~CLF, y = ~CLS, z = ~Sphere,
        color = ~tag, colors = usr.col, text= ~paste('Symbol:',name,'<br>Probe:',Probe)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'CLF'),
                        yaxis = list(title = 'CLS'),
                        zaxis = list(title = 'Sphere')))





#use the package seriation to reordering the cluster
#OLO GW mean non
library(AnnotationHub)
ah <- AnnotationHub()
human <- ah[ah$species == "Homo sapiens"]
human.orgs <- subset(human, human$rdataclass == "OrgDb")
data.human.orgs <- human.orgs[["AH49582"]]
genesym <- total.probe.data.frame$name[!is.na(total.probe.data.frame$name)]
geneid <- select(data.human.orgs, keys=genesym, keytype="SYMBOL",
                 columns="ENTREZID")
geneid <- geneid %>% dplyr::rename(name=SYMBOL)


# Use panther for further exploration 
geneid <- unique(geneid)
geneid.deduplicate <- geneid[!duplicated(geneid$name),]
annotation.probe.dataframe <- left_join(total.probe.data.frame, geneid.deduplicate, by="name")


# up intersect
up.annotation.tukey <- annotation.probe.dataframe %>% filter(Probe %in% intersect.up)
library(PANTHER.db)
pthOrganisms(PANTHER.db) <- "HUMAN"
key     <- up.annotation.tukey$ENTREZID
choices <- c("CLASS_TERM")
#protein function
up.result.protein.function <- PANTHER.db::select(PANTHER.db,
                             keys = key,
                             columns = choices,
                             keytype = "ENTREZ")
up.result.protein.function <- up.result.protein.function %>% filter(!is.na(CLASS_TERM))
#pathway
choices <- c("PATHWAY_TERM")
up.result.protein.pathway <- PANTHER.db::select(PANTHER.db,
                                              keys = key,
                                              columns = choices,
                                              keytype = "ENTREZ")
up.result.protein.pathway <- up.result.protein.pathway %>% filter(!is.na(PATHWAY_TERM))


# down intersect
down.annotation.tukey <- annotation.probe.dataframe %>% filter(Probe %in% intersect.down)




key     <- down.annotation.tukey$ENTREZID
choices <- c("CLASS_TERM")
#protein function
down.result.protein.function <- PANTHER.db::select(PANTHER.db,
                                              keys = key,
                                              columns = choices,
                                              keytype = "ENTREZ")
down.result.protein.function <- down.result.protein.function %>% filter(!is.na(CLASS_TERM))
#pathway
choices <- c("PATHWAY_TERM")
down.result.protein.pathway <- PANTHER.db::select(PANTHER.db,
                                             keys = key,
                                             columns = choices,
                                             keytype = "ENTREZ")
down.result.protein.pathway <- down.result.protein.pathway %>% filter(!is.na(PATHWAY_TERM))

save(up.result.protein.function, up.result.protein.pathway, down.result.protein.function, down.result.protein.pathway, file = "tukey_result_panther.Rdata")




#final annotation big table
up.annotation.tukey <- up.annotation.tukey %>% rename(ENTREZ = ENTREZID)
up.result.protein.function 
up.result.protein.pathway

up.result.aggregate.protein.function <- aggregate(CLASS_TERM ~ ENTREZ,data=up.result.protein.function,toString) 
up.result.aggregate.protein.pathway <- aggregate(PATHWAY_TERM ~ ENTREZ,data=up.result.protein.pathway,toString) 

up.total.annotation.tukey <- left_join(up.annotation.tukey,up.result.aggregate.protein.function, by="ENTREZ")
up.total.annotation.tukey <- left_join(up.total.annotation.tukey, up.result.aggregate.protein.pathway, by ="ENTREZ")


down.annotation.tukey <- down.annotation.tukey %>% rename(ENTREZ = ENTREZID)
down.result.aggregate.protein.function <- aggregate(CLASS_TERM ~ ENTREZ, data=down.result.protein.function,toString) 
down.result.aggregate.protein.pathway  <- aggregate(PATHWAY_TERM ~ ENTREZ, data=down.result.protein.pathway,toString) 

down.total.annotation.tukey <- left_join(down.annotation.tukey,down.result.aggregate.protein.function, by="ENTREZ")
down.total.annotation.tukey <- left_join(down.total.annotation.tukey, down.result.aggregate.protein.pathway, by ="ENTREZ")


up.p.value.CLF_CLS <- up.CLF_CLS %>% filter(Probe %in% intersect.up) %>%
    mutate(CLF_CLS.BH.p.adj=BH.p.adj) %>%
    dplyr::select(Probe,CLF_CLS.BH.p.adj)
up.p.value.Sphere_CLS <- up.Sphere_CLS %>% filter(Probe %in% intersect.up) %>%
    mutate(Sphere_CLS.BH.p.adj=BH.p.adj) %>%
    dplyr::select(Probe,Sphere_CLS.BH.p.adj)
up.total.annotation.tukey <- left_join(up.total.annotation.tukey, up.p.value.CLF_CLS, by="Probe")
up.total.annotation.tukey <- left_join(up.total.annotation.tukey, up.p.value.Sphere_CLS, by="Probe")
up.total.annotation.tukey <- up.total.annotation.tukey %>% dplyr::select(-c(BP,CC,GO))


down.p.value.CLF_CLS <- down.CLF_CLS %>% filter(Probe %in% intersect.down) %>%
    mutate(CLF_CLS.BH.p.adj=BH.p.adj) %>%
    dplyr::select(Probe,CLF_CLS.BH.p.adj)
down.p.value.Sphere_CLS <- down.Sphere_CLS %>% filter(Probe %in% intersect.down) %>%
    mutate(Sphere_CLS.BH.p.adj=BH.p.adj) %>%
    dplyr::select(Probe,Sphere_CLS.BH.p.adj)
down.total.annotation.tukey <- left_join(down.total.annotation.tukey, down.p.value.CLF_CLS, by="Probe")
down.total.annotation.tukey <- left_join(down.total.annotation.tukey, down.p.value.Sphere_CLS, by="Probe")
down.total.annotation.tukey <- down.total.annotation.tukey %>% dplyr::select(-c(BP,CC,GO))

write_csv(up.total.annotation.tukey, path = "upregulation_annotation_tukey_final.csv", col_names = T)
write_csv(down.total.annotation.tukey, path = "downregulation_annotation_tukey_final.csv", col_names = T)

write_csv(up.total.annotation.tukey %>% filter(Probe %in% intersect.up), path = "upregulation_00001_annotation_tukey_final.csv", col_names = T)
write_csv(down.total.annotation.tukey %>% filter(Probe %in% intersect.down), path = "downregulation_00001_annotation_tukey_final.csv", col_names = T)

up.total.annotation.tukey %>% filter(Probe %in% intersect.up)
down.total.annotation.tukey %>% filter(Probe %in% intersect.down)

#merge with the origin p.value


up.total.annotation.tukey %>% head






#experiment
a <- data.frame(x1=c("1293","2341","2345"),x2=c("A","B","C"))
b <- data.frame(x1=c("1293", "1293", "1293"),x3=c("function1","function2","function3"))




# ***************************************** -------------------------------
#find the direction
aov.norm <- norm %>% as.data.frame %>% mutate(CLF = (`CHW_CLF 13-6-1.CEL`+`CHW_CLF 13-6-2.CEL` + `CHW_CLF 13-6-3.CEL`)/3,
                                  CLS = (`CHW_CLS 1-2 p.6-1.CEL` + `CHW_CLS 1-2 p.6-2.CEL`+`CHW_CLS 1-2 p.6-3.CEL`)/3,
                                  Sphere = (`CHW_Sphere-1.CEL`+`CHW_Sphere-2.CEL`+ `CHW_Sphere-3.CEL`)/3) %>% dplyr::select(CLS, Sphere, CLF)
aov.norm <- aov.norm %>% mutate(Probe = probe.name)
aov.norm <- aov.norm %>% mutate(CLF_CLS = CLF - CLS, Sphere_CLS = Sphere - CLS)
aov.norm <- aov.norm %>% dplyr::select(-c(CLF,CLS,Sphere))
# QN, Log2, Annova, BF ----------------------------------------------------

aov.bf.final <- left_join(aov.bonferroni.final, aov.norm, by = "Probe")


up.aov.bf.CLF_CLS       <- aov.bf.final %>% filter(case == 'CLS-CLF', CLF_CLS > 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "up") %>% dplyr::select(-c(group1,group2))
down.aov.bf.CLF_CLS     <- aov.bf.final %>% filter(case == 'CLS-CLF', CLF_CLS < 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "down") %>% dplyr::select(-c(group1,group2))


up.aov.bf.Sphere_CLS    <- aov.bf.final %>% filter(case == 'Sphere-CLS', Sphere_CLS > 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "up") %>% dplyr::select(-c(group1,group2))
down.aov.bf.Sphere_CLS  <- aov.bf.final %>% filter(case == 'Sphere-CLS', Sphere_CLS < 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "down") %>% dplyr::select(-c(group1,group2))



total.aov.bf.final <- bind_rows(up.aov.bf.CLF_CLS, down.aov.bf.CLF_CLS, up.aov.bf.Sphere_CLS, down.aov.bf.Sphere_CLS)

total.aov.bf.final %>% group_by(case, type) %>% summarise(n_0.05 = sum(BH.mul.adj.p < 0.05),
                                                          n_0.01 = sum(BH.mul.adj.p < 0.01),
                                                          n_0.001 = sum(BH.mul.adj.p < 0.001),
                                                          n_0.0001 = sum(BH.mul.adj.p < 0.0001),
                                                          n_0.00001 = sum(BH.mul.adj.p < 0.00001))


draw.venndiagram <- function(total.aov.final, p.cutoff = 0.0001, direction = "up"){
    
    total.aov.bf.final <- total.aov.final
    BH.mul.adj.p.cutoff <- p.cutoff
    up.CLF_CLS    <- total.aov.bf.final %>% filter(BH.mul.adj.p < BH.mul.adj.p.cutoff, case == 'CLS-CLF', type == direction) %>% with(Probe)
    up.Sphere_CLS <- total.aov.bf.final %>% filter(BH.mul.adj.p < BH.mul.adj.p.cutoff, case == 'Sphere-CLS', type == direction) %>% with(Probe)
    
    title.area1 <- paste("P6+fibroblast", switch(direction, up = ">", down = "<"),"P6")
    title.area2 <- paste("P6-Sphere", switch(direction, up = ">", down = "<"), "P6")
    
    area.1 <- length(up.CLF_CLS  )
    area.2 <- length(up.Sphere_CLS)
    area.12 <- length(intersect(up.CLF_CLS, up.Sphere_CLS))
    area.list <- list(area.1, area.2, area.12)
    draw.pairwise.venn(area1 = area.list[[1]], area2 = area.list[[2]], cross.area = area.list[[3]],
                       category = c(title.area1, title.area2),
                       fill = c("Blue","Red"))}

# Upregulation
draw.venndiagram(total.aov.final = total.aov.bf.final, p.cutoff = 0.0001, direction = "up" )
# Downregulation
draw.venndiagram(total.aov.final = total.aov.bf.final, p.cutoff = 0.0001, direction = "down" )





# QN, Log2, Annova, Holm --------------------------------------------------

aov.hm.final <- left_join(aov.holm.final, aov.norm, by = "Probe")


up.aov.hm.CLF_CLS       <- aov.hm.final %>% filter(case == 'CLS-CLF', CLF_CLS > 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "up") %>% dplyr::select(-c(group1,group2))
down.aov.hm.CLF_CLS     <- aov.hm.final %>% filter(case == 'CLS-CLF', CLF_CLS < 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "down") %>% dplyr::select(-c(group1,group2))


up.aov.hm.Sphere_CLS    <- aov.hm.final %>% filter(case == 'Sphere-CLS', Sphere_CLS > 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "up") %>% dplyr::select(-c(group1,group2))
down.aov.hm.Sphere_CLS  <- aov.hm.final %>% filter(case == 'Sphere-CLS', Sphere_CLS < 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "down") %>% dplyr::select(-c(group1,group2))



total.aov.hm.final <- bind_rows(up.aov.hm.CLF_CLS, down.aov.hm.CLF_CLS, up.aov.hm.Sphere_CLS, down.aov.hm.Sphere_CLS)

total.aov.hm.final %>% group_by(case, type) %>% summarise(n_0.05 = sum(BH.mul.adj.p < 0.05),
                                                          n_0.01 = sum(BH.mul.adj.p < 0.01),
                                                          n_0.001 = sum(BH.mul.adj.p < 0.001),
                                                          n_0.0001 = sum(BH.mul.adj.p < 0.0001),
                                                          n_0.00001 = sum(BH.mul.adj.p < 0.00001))
# Upregulation
draw.venndiagram(total.aov.final = total.aov.hm.final, p.cutoff = 0.0001, direction = "up" )
# Downregulation
draw.venndiagram(total.aov.final = total.aov.hm.final, p.cutoff = 0.0001, direction = "down" )

# QN, Log2, Annova, BH ----------------------------------------------------

aov.bjh.final <- left_join(aov.bamjamin.final, aov.norm, by = "Probe")


up.aov.bjh.CLF_CLS       <- aov.bjh.final %>% filter(case == 'CLS-CLF', CLF_CLS > 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "up") %>% dplyr::select(-c(group1,group2))
down.aov.bjh.CLF_CLS     <- aov.bjh.final %>% filter(case == 'CLS-CLF', CLF_CLS < 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "down") %>% dplyr::select(-c(group1,group2))


up.aov.bjh.Sphere_CLS    <- aov.bjh.final %>% filter(case == 'Sphere-CLS', Sphere_CLS > 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "up") %>% dplyr::select(-c(group1,group2))
down.aov.bjh.Sphere_CLS  <- aov.bjh.final %>% filter(case == 'Sphere-CLS', Sphere_CLS < 0) %>% mutate(BH.mul.adj.p = p.adjust(p.value, method ="BH"), type = "down") %>% dplyr::select(-c(group1,group2))



total.aov.bjh.final <- bind_rows(up.aov.bjh.CLF_CLS, down.aov.bjh.CLF_CLS, up.aov.bjh.Sphere_CLS, down.aov.bjh.Sphere_CLS)

total.aov.bjh.final %>% group_by(case, type) %>% summarise(n_0.05 = sum(BH.mul.adj.p < 0.05),
                                                           n_0.01 = sum(BH.mul.adj.p < 0.01),
                                                           n_0.001 = sum(BH.mul.adj.p < 0.001),
                                                           n_0.0001 = sum(BH.mul.adj.p < 0.0001),
                                                           n_0.00001 = sum(BH.mul.adj.p < 0.00001))

save(total.aov.bf.final, total.aov.bjh.final, total.aov.hm.final, file = "total_aov_updown_result.Rdata")
# Upregulation
draw.venndiagram(total.aov.final = total.aov.bjh.final, p.cutoff = 0.0001, direction = "up" )
# Downregulation
draw.venndiagram(total.aov.final = total.aov.bjh.final, p.cutoff = 0.0001, direction = "down" )

# SOM -------------------------------------------------------------------

# Clustering & Heatmap ----------------------------------------------------




