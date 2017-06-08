# curated pathway
library(graphite)
edges <- data.frame(src=c("672","1"), dest=c("7157","11"), direction="undirected", type="binding")
pathway.1 <- buildPathway("#1", "example_1", edges, "hsapiens", "database", "ENTREZID")
pathway.2 <- buildPathway("#2", "example_2", edges, "hsapiens", "database", "ENTREZID")
prepareSPIA(list(pathway.1), pathwaySetName="test", print.names=FALSE)
pathwayDatabases()
reactome.pathway<- pathways("hsapiens", "reactome")




library(rBiopaxParser)
reactome.biopax <- readBiopax(file="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/biopax_level_3/Homo_sapiens.owl")

library(org.Hs.eg.db)
IGF1R.signaling.cascade <- path.info$'IGF1R signaling cascade'
Gene.products <- IGF1R.signaling.cascade$nodes

IGF1R.signaling.cascade.related <- AnnotationDbi::select(org.Hs.eg.db,
                      keys = Gene.products,
                      columns = c("SYMBOL","UNIPROT"),
                      keytype = "ENTREZID")

test<- path.info[1]$`2-LTR circle formation`

IGF.regulation.IGFBP <- path.info[['Regulation of IGF Activity by IGFBP']]   #R-HSA-381426
IGF1R.signaling.cascade <- path.info$'IGF1R signaling cascade'   #R-HSA-2428924

library(org.Hs.eg.db)


protein.in.164843 <- AnnotationDbi::select(org.Hs.eg.db,
                      keys = path.info[1]$`2-LTR circle formation`$nodes %>% as.character,
                      columns = "UNIPROT",
                      keytype = "ENTREZID")

result.entrzid <- AnnotationDbi::select(org.Hs.eg.db,
                      keys = result$Identifier %>% as.character,
                      columns = "ENTREZID",
                      keytype = "UNIPROT")
result.entrzid <- result.entrzid %>% filter(!is.na(ENTREZID)) %>% with(ENTREZID)
# input json file from neo4j query result ---------------------------------
json.path<-"/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/neo4j_query"
library(rjson)

pathway.id <- "R-HSA-2428924"
input.reactome     <- fromJSON(file=paste0(json.path,"/",pathway.id,"_input.json"))
output.reactome    <- fromJSON(file=paste0(json.path,"/",pathway.id,"_output.json"))
regulator.reactome <- fromJSON(file=paste0(json.path,"/",pathway.id,"_reguation.json"))
regulated.reactome <- fromJSON(file=paste0(json.path,"/",pathway.id,"_eventbyreguation.json"))


getInputData <- function(input.reactome){
colname.name <- input.reactome$results %>% map(`[[`,1) %>% unlist
data.list    <- input.reactome$results %>% map(`[[`,2) %>% transpose() %>% at_depth(2,"row") 

Pathwayname          <-  data.list %>% at_depth(2,1) %>% unlist
ReactionLikeEvent    <-  data.list %>% at_depth(2,2) %>% unlist
ReactionType         <-  data.list %>% at_depth(2,3) %>% unlist
ReactionLikeEvent_ID <-  data.list %>% at_depth(2,4) %>% unlist
PhysicalEntity       <-  data.list %>% at_depth(2,5) %>% unlist
Identifier           <-  data.list %>% at_depth(2,6) %>% unlist

result <- data.frame(Pathwayname, ReactionLikeEvent, ReactionType, ReactionLikeEvent_ID, PhysicalEntity, Identifier) %>% mutate(DataType="Input")


return(result)
}


getOutputData <- function(output.reactome){
    colname.name <- output.reactome$results %>% map(`[[`,1) %>% unlist
    data.list    <- output.reactome$results %>% map(`[[`,2) %>% transpose() %>% at_depth(2,"row") 
    
    Pathwayname          <-  data.list %>% at_depth(2,1) %>% unlist
    ReactionLikeEvent    <-  data.list %>% at_depth(2,2) %>% unlist
    ReactionType         <-  data.list %>% at_depth(2,3) %>% unlist
    ReactionLikeEvent_ID <-  data.list %>% at_depth(2,4) %>% unlist
    PhysicalEntity       <-  data.list %>% at_depth(2,5) %>% unlist
    Identifier           <-  data.list %>% at_depth(2,6) %>% unlist
    
    result <- data.frame(Pathwayname, ReactionLikeEvent, ReactionType, ReactionLikeEvent_ID, PhysicalEntity, Identifier) %>% mutate(DataType="Output")
    
    
    return(result)
}


# Browse the Graphite curated pathway -------------------------------------



IGF1R.signaling.cascade.symbol <-     AnnotationDbi::select(org.Hs.eg.db,
                                                             keys = IGF1R.signaling.cascade$nodes %>% as.character,
                                                             columns = "SYMBOL",
                                                             keytype = "ENTREZID")
IGF.regulation.IGFBP.symbol <-         AnnotationDbi::select(org.Hs.eg.db,
                                                             keys = IGF.regulation.IGFBP$nodes %>% as.character,
                                                             columns = "SYMBOL",
                                                             keytype = "ENTREZID")
# translated the symbol to uniprot
translate.vector <-vector()
for (i in 1:length(IGF1R.signaling.cascade.symbol$SYMBOL)){
    translate.vector[IGF1R.signaling.cascade.symbol$ENTREZID[i]] <- IGF1R.signaling.cascade.symbol$SYMBOL[i]
}


getTranslateEntrezToSymbol <- function(symbol.list){
    translate.vector <-vector()
    for (i in 1:length(symbol.list$SYMBOL)){
        translate.vector[symbol.list$ENTREZID[i]] <- symbol.list$SYMBOL[i]
    }    
    return(translate.vector)
}
getGraphTranslate <- function(graph,translate.vector){
    rownames(graph) <- unname(translate.vector[rownames(graph)])
    colnames(graph) <- unname(translate.vector[colnames(graph)])
    return(graph)
}

translate.vector <- getTranslateEntrezToSymbol(symbol.list = IGF.regulation.IGFBP.symbol)

IGF1R.signaling.cascade.graph <- getGraphTranslate(graph=IGF1R.signaling.cascade.graph, translate.vector = translate.vector)

IGF.regulation.IGFBP.graph <- IGF.regulation.IGFBP$`binding/association` 
IGF.regulation.IGFBP.graph <- getGraphTranslate(graph=IGF.regulation.IGFBP.graph, translate.vector = translate.vector)

# Observation the data ----------------------------------------------------------------
#load the library
library(statnet)

IGF1R.signaling.cascade.graph <- IGF1R.signaling.cascade$`binding/association`

#Creat a network object with sociometrix data format
net1 <- network(IGF1R.signaling.cascade.graph, matrix.type="adjacency")
net1 <- network(IGF.regulation.IGFBP.graph, matrix.type="adjacency")
##Visualization
gplot(net1)
gplot(net1, vertex.col = 5, displaylabels =TRUE)

net1.edge.list <- as.matrix(net1,matrix.type = "edgelist")
node.number <- attr(net1.edge.list,"n") 
node.list <- attr(net1.edge.list,"vnames")

attr(net1.edge.list,"n") <- NULL
attr(net1.edge.list,"vnames") <- NULL
net1.edge.dataframe <- net1.edge.list %>% as.data.frame()
net1.node.dataframe <- data.frame(1:node.number, node.list)
colnames(net1.node.dataframe) <- c("Id", "Symbol")
write_delim(net1.edge.dataframe, path="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/neo4j_cytoscape/net1_interaction.txt", delim="\t")
write_delim(net1.node.dataframe, path="/Users/Weitinglin/Documents/R_scripts/Lab/Reactome_database/neo4j_cytoscape/net1_node.txt", delim="\t")
