# The Reciprocal Metabolic Reprogramming Caused by The  Cancer-Associated Fibroblast in Primary Non-Small Cell Lung Cancer Stem Cell   
Lin, Ting-Wei

## Biological Question
1. What's the impact of Cancer-associated fibroblast on the cancer stem cell?
2. Which part of these impact determined the drug-resistant and migration characteristics in the lung cancer stem cell?
3. What hypothesis can be generated for these phenomenon from our data (microarray, RNAseq, Proteomics)


## Thinking on methology approach
1. How you deal with your microarray data?
  - Preprocess method:MAS5
  - Normalization: Qunatile Normalization
  - Differential Expression: one way t-test
  - Downstreat analysis: 
    - Annotation with panther database, GO database
    - Enrichment analysis with gProfiler
    - Clustering and SOM method to find the related gene with the DE gene
2. How you deal with RNAseq data?
 - STAR, StringTie, Ballgown
 - Trinity for unmapped read
3. How you deal with the Proteomic data?
4. How you manage whole the project?
