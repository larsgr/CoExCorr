source("R/phyloXML.R")
library(tidyverse)
library(ape)


### Load tree and taxonomy data ########
#
#

treeData <- readRDS("data/treeData/treeData.RDS")

# give the trees names..
names(treeData) <- sprintf("tree%05i",seq_along(treeData))

# load ncbi taxon ID mappings
taxidTbl <- read_tsv(file = "indata/taxid2taxnames.txt")

# tax ids for Os and At
OsAtTaxid <- c(At=3702,Os=39947)

# Find out which trees that contain both Os and At
sapply(treeData, function(td){
  all(OsAtTaxid %in% td$tip.data$taxonID)
}) -> hasOsAt

# make geneID to treeID mapping
lapply(OsAtTaxid, function(taxid){
  
  lapply(treeData[hasOsAt], function(td){
    td$tip.label[td$tip.data$taxonID == taxid]
  }) -> tree_geneID_list
  
  setNames( rep(names(tree_geneID_list),sapply(tree_geneID_list,length)),
            unlist(tree_geneID_list))
  
}) -> geneID2treeID



### Get geneIDs for treesWithAtOs ##### 
# 
# get the gene IDs for genes in trees with both At and Os that also have 
# expression data.
#
dir.create("data/subsets/treesWithAtOs")

At_geneIDs <- names(geneID2treeID$At)
At_geneIDs <- At_geneIDs[At_geneIDs %in% rownames(expAt)]
saveRDS(At_geneIDs, file = "data/subsets/treesWithAtOs/At_geneIDs.RDS")

Os_geneIDs <- names(geneID2treeID$Os)
Os_geneIDs <- Os_geneIDs[Os_geneIDs %in% rownames(expOs)]
saveRDS(Os_geneIDs, file = "data/subsets/treesWithAtOs/Os_geneIDs.RDS")



### load CCS #######

CCS <- readRDS("data/subsets/treesWithAtOs/AtOs_CCS.RDS")
CCSAt <- readRDS("data/subsets/treesWithAtOs/At_selfCCS.RDS")
CCSOs <- readRDS("data/subsets/treesWithAtOs/Os_selfCCS.RDS")
source("R/calcRanks.R")

AtOs11orthos <- lapply(
  c(At="data/subsets/treesWithAtOs/At_AtOs11_geneIDs.RDS",
    Os="data/subsets/treesWithAtOs/Os_AtOs11_geneIDs.RDS"),
  readRDS)



### Candidate tree ###################
# find some interresting candidate trees
# Find a tree with several 1:1 genes

AtOs11Orthos <- readRDS("data/subsets/PODC_AtOs/At_Os_11orthos.RDS")
as_tibble(table(geneID2treeID$At[AtOs11Orthos$At])) %>% 
  filter(n>1) %>% .$Var1 -> treesWithSeveral11


td <- treeData[[treesWithSeveral11[2]]]
tree <- as(treeDataToPhylo4(td),"phylo")


# add taxon on internal nodes
tree$node.label <- taxidTbl$taxname[match(td$node.data$taxonID,taxidTbl$taxid)]

# remove all other species
tr <- drop.tip(phy = tree, tip = which(!(td$tip.data$taxonID %in% OsAtTaxid)))



### helper functions ####

# change tip.label from geneIDs to geneNames
nameTree <- function(phylo){
  phylo$tip.label <- getGeneNames(phylo$tip.label)
  return(phylo)
}

# change column and rownames from geneIDs to geneNames
nameMat <- function(m){
  dimnames(m) <- lapply(dimnames(m), getGeneNames)
  return(m)
}


orthoPairRnks <- function(geneIDs, CCS, rename=T){
  if(!rename)
    getGeneNames <- function(x){x}
  
  # get all combinations of ortholog pairs
  allOrthoPairs <- 
    expand.grid(rowGenes=geneIDs[geneIDs %in% rownames(CCS)],
                colGenes=geneIDs[geneIDs %in% colnames(CCS)],
                stringsAsFactors = F)

  # get ranks for each pair  
  allOrthoPairs$rnks <- quickRanks(CCS,allOrthoPairs$rowGenes,allOrthoPairs$colGenes)
  
  # spread out like a matrix
  as_tibble(allOrthoPairs) %>% 
    mutate( rowGenes = getGeneNames(rowGenes) ) %>% 
    mutate( colGenes = getGeneNames(colGenes) ) %>% 
    spread(key = rowGenes,value = rnks)
}


paralogPCC <- function(geneIDs, expMat, rename=T){
  if(!rename)
    nameMat <- function(x){x}
  nameMat(cor(t(expMat[geneIDs[geneIDs %in% rownames(expMat)], ])))
}


#### load metadata #############
#
#
metaOs <- readxl::read_excel("indata/oryza_sativa_sample.xlsx",col_types = rep("text",15))
metaAt <- readxl::read_excel("indata/arabidopsis_thaliana_sample.xlsx",col_types = rep("text",15))
PO <- data_frame( 
  Tissue = c("PO:0005026","PO:0000025","PO:0003011","PO:0009005","PO:0000256","PO:0000263",
         "PO:0009006","PO:0025034","PO:0020103","PO:0000014",
         "PO:0009049","PO:0009046","PO:0000056",
         "PO:0000003"),
  generalTissue = c(rep("root",6),rep("leaf",4),rep("flower",3),"whole"),
  color = c(rep("blue",6),rep("darkgreen",4),rep("red",3),"green"))

metaOs %>%  left_join(PO,by="Tissue") -> metaOs
metaOs$color[is.na(metaOs$color)] <- "gray"
metaAt %>%  left_join(PO,by="Tissue") -> metaAt
metaAt$color[is.na(metaAt$color)] <- "gray"

# reorder the At metadata to match the expression matrix
metaAt <- metaAt[match(colnames(expAt),metaAt$Run),]




#### Analyse tree #########


#
# plot tree
#

plot.phylo(nameTree(tr),show.node.label = T)

#
# Investigate ORS and PRS
#

orthoPairRnks(tr$tip.label,CCS)
orthoPairRnks(tr$tip.label,CCSOs)
orthoPairRnks(tr$tip.label,CCSAt)

paralogPCC( tr$tip.label, expOs,rename = F)
paralogPCC( tr$tip.label, expAt)


plot(expAt["AT5G42650", ],expAt["AT4G15440", ], col=metaAt$color)
plot(as.data.frame(t(expOs[tr$tip.label[tr$tip.label %in% rownames(expOs)], ])), 
     pch=20,cex=0.1,col=metaOs$color)

# ORS indicates that OS02G0110200 and to a lesser degree OsCYP74A4 has the most 
# similar expression pattern with both AtCYP74A and AtCYP74B2.
# Paralog rank score (PRS) indicates that OsCYP74A1/2/3 are similar.
# Also to a lesser degree OS02G0110200 and OsCYP74A4 are similar.
# interrestinlgy OS02G0110200 has a mutually exclusive expression with CYP74A3 and to some degree the others

# OS03G0767000 CYP74A1
# AT5G42650 CYP74A

# OS02G0218800/OS02G0218700 a.k.a. CYP74A3/4 are tandem duplicates but their expression does not correlate very well
# OS03G0225900 a.k.a. CYP74A2 


# AT4G15440 a.k.a. CYP74B2
# OS02G0110200 ?

