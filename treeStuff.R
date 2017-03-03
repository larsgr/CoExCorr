source("R/phyloXML.R")
source("R/getGeneNames.R")
source("R/calcRanks.R")
library(tidyverse)
library(ape)


### Load tree and taxonomy data ########
#
#

treeData <- readRDS("data/treeData/treeData.RDS")

# give the trees names..
names(treeData) <- sprintf("tree%05i",seq_along(treeData))

# load ncbi taxon ID mappings
taxidTbl <- read_tsv(file = "indata/taxid2taxnames.txt", col_types = cols())


# tax ids for our species
spc2taxid <- c(
  At=3702,
  Os=39947,
  Gm=3847,
  Zm=4577,
  Sl=4081)

taxidTbl$spc <- names(spc2taxid)[match(taxidTbl$taxonID,spc2taxid)]


# make a table of taxon counts per tree
taxCount <-
  treeData %>% 
  map( ~ as.character(.x$tip.data$taxonID)) %>% 
  data_frame( taxon = ., treeID = names(.)) %>% 
  unnest() %>% 
  with(., table(treeID,taxon)) %>% 
  as.matrix()


# Find out which trees that contain at least two of our species
hasTwoSpc <- rowSums(taxCount[ ,as.character(spc2taxid)] > 0) >= 2

# make geneID to treeID mapping 
# For each species
lapply(spc2taxid, function(taxid){
  
  # for each tree get geneIDs for this species
  lapply(treeData[hasTwoSpc], function(td){
    td$tip.label[td$tip.data$taxonID == taxid]
  }) -> tree_geneID_list
  
  
  setNames( rep(names(tree_geneID_list),sapply(tree_geneID_list,length)),
            unlist(tree_geneID_list))
  
}) -> geneID2treeID

## Functions ####
# Draw  triangles on duplication nodes
markDupNodes <- function(p4d){
  isDup <- nodeData(p4d)$isDuplication
  nodelabels(node = nodeId(p4d,type="internal")[isDup],pch=2, col="red",
             cex=pmax(nodeData(p4d)$duplicationConfidence[isDup],0.2))
}

# change tip.label from geneIDs to geneNames
nameTree <- function(phylo){
  phylo$tip.label <- getGeneNames(geneIDs = phylo$tip.label)
  return(phylo)
}

# change column and rownames from geneIDs to geneNames
nameMat <- function(m){
  dimnames(m) <- lapply(dimnames(m), getGeneNames)
  return(m)
}



paralogPCC <- function(geneIDs, expMat, rename=T){
  if(!rename)
    nameMat <- function(x){x}
  nameMat(cor(t(expMat[geneIDs[geneIDs %in% rownames(expMat)], ])))
}

# get phylo from treeData
getPhylo <- function(p4d){
  suppressWarnings(
    as(p4d,"phylo")
  )
}

fixTreeData <- function(td){
  # no factors please
  td$tip.data$taxonID <- as.integer(as.character(td$tip.data$taxonID))
  td$node.data$taxonID <- as.integer(as.character(td$node.data$taxonID))
  
  # add spc and taxonName to data
  td$tip.data <- 
    td$tip.data %>% left_join(taxidTbl, by="taxonID")
  
  td$node.data <- 
    td$node.data %>% left_join(taxidTbl, by="taxonID")
  
  # set row names as phylobase likes it
  row.names(td$tip.data) <- td$tip.label
  row.names(td$node.data) <- 1:nrow(td$node.data) + nrow(td$tip.data)
  
  return(td)
}

orthoPairRnks <- function(geneIDs, ccs){
  
  # get all combinations of ortholog pairs
  allOrthoPairs <- 
    expand.grid(rowGenes=geneIDs[geneIDs %in% rownames(ccs)],
                colGenes=geneIDs[geneIDs %in% colnames(ccs)],
                stringsAsFactors = F)
  
  # get ranks for each pair  
  allOrthoPairs$rnks <- quickRanks(ccs,allOrthoPairs$rowGenes,allOrthoPairs$colGenes)
  
  # spread out like a matrix
  rankTbl <-
    as_tibble(allOrthoPairs) %>%
    spread(key = colGenes,value = rnks)
  rankMat <- as.matrix(rankTbl[,-1])
  rownames(rankMat) <- rankTbl$rowGenes
  
  return(rankMat)
}

getRankMat <- function(p4d){
  geneIDs <- tipLabels(p4d) # !! order of the tips in the plot doesnt match
  
  M <- matrix(.0,ncol=nTips(p4d),nrow=nTips(p4d),dimnames=list(geneIDs,geneIDs))
  
  for(ccs in CCS){
    m <- orthoPairRnks(geneIDs, ccs)
    M[rownames(m),colnames(m)] <- m
    M[colnames(m),rownames(m)] <- t(m)
  }
  return(M)
}

### load CCS #######

# CCS_old <- readRDS("data/subsets/treesWithAtOs/AtOs_CCS.RDS")
system.time({
  CCS <-
    dir("data/CCS",pattern="CCS",full.names = T) %>% 
    set_names( sub(".*/(.{4})_.*","\\1",.)) %>% 
    map( readRDS )
})




 




### Candidate tree ###################
# find some interresting candidate trees
# Find a tree with several 1:1 genes

AtOs11orthos <- lapply(
  c(At="data/CCS/At_AtOs11_geneIDs.RDS",
    Os="data/CCS/Os_AtOs11_geneIDs.RDS"),
  readRDS)

as_tibble(table(geneID2treeID$At[AtOs11orthos$At])) %>% 
  filter(n>1) %>% .$Var1 -> treesWithSeveral11


# alternative smaller tree for testing
# treeID <- which(rowSums(taxCount) < 15 &
#              rowSums(taxCount[,as.character(spc2taxid)]>0) > 2)[1]




# select tree of interrest
treeID <- treesWithSeveral11[1]
# convert to p4d
p4d <- treeDataToPhylo4d(fixTreeData(treeData[[treeID]]))
# remove tips from species without data
p4d <- p4d[tipData(p4d)$taxonID %in% spc2taxid[c("At","Os","Zm","Sl")]]

td <- fixTreeData(treeData[[treeID]])
dim(td$m)
dim(td$tip.data)
str(td)


# get the phylo version of the tree
tr <- getPhylo(p4d)

plot.phylo( nameTree(phylo =  tr))
nodelabels(nodeData(p4d)$taxonName,frame = "none",adj=c(0,0.5),cex=0.6)
markDupNodes(p4d)


## Extract rank matrix ####
#

M <- getRankMat(p4d)
heatmap(nameMat(M)^2,Rowv = NA,Colv = NA, scale="none")
heatmap(nameMat(M)^2,scale="none")







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

plot.phylo( nameTree(phylo =  tr))
nodelabels(nodeData(p4d)$taxonName,frame = "none",adj=c(0,0.5),cex=0.6)
markDupNodes(p4d)



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

