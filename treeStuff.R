source("R/phyloXML.R")
source("R/getGeneNames.R")
source("R/calcRanks.R")
library(tidyverse)
library(ape)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(plotly)



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

ggplotlyHeatmap <- function(mData){
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  
  p <-
    ggplot(melt(mData),aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colours = myPalette(100),limits=c(0, 4)) +
    coord_equal() +
    theme_bw()
  
  ggplotly(p)
}



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


# orthoPairRnks <- function(geneIDs, ccs){
#   
#   # get all combinations of ortholog pairs
#   allOrthoPairs <- 
#     expand.grid(rowGenes=geneIDs[geneIDs %in% rownames(ccs)],
#                 colGenes=geneIDs[geneIDs %in% colnames(ccs)],
#                 stringsAsFactors = F)
#   
#   # get ranks for each pair  
#   allOrthoPairs$rnks <- quickRanks(ccs,allOrthoPairs$rowGenes,allOrthoPairs$colGenes)
#   
#   # spread out like a matrix
#   rankTbl <-
#     as_tibble(allOrthoPairs) %>%
#     spread(key = colGenes,value = rnks)
#   rankMat <- as.matrix(rankTbl[,-1])
#   rownames(rankMat) <- rankTbl$rowGenes
#   
#   return(rankMat)
# }

# getRankMat <- function(p4d){
#   geneIDs <- tipLabels(p4d) # !! order of the tips in the plot doesnt match
#   
#   M <- matrix(NA_real_,ncol=nTips(p4d),nrow=nTips(p4d),dimnames=list(geneIDs,geneIDs))
#   
#   for(ccs in CCS){
#     m <- orthoPairRnks(geneIDs, ccs)
#     M[rownames(m),colnames(m)] <- m
#     M[colnames(m),rownames(m)] <- t(m)
#   }
#   return(M)
# }

# get ranks from precalculated rank tables
getRankMat <- function(p4d){
  geneIDs <- tipLabels(p4d) # !! order of the tips in the plot doesnt match
  
  M <- matrix(NA_real_,ncol=nTips(p4d),nrow=nTips(p4d),dimnames=list(geneIDs,geneIDs))
  
  for(rnkTbl in orthoRankTables){
    rowGenes <- geneIDs[geneIDs %in% rnkTbl$SPC1]
    colGenes <- geneIDs[geneIDs %in% rnkTbl$SPC2]
    
    if( length(rowGenes)>0 & length(colGenes)>0){
      M[rowGenes,colGenes] <- sapply(colGenes,function(spc2gene){
        rnkTblTmp <- rnkTbl[rnkTbl$SPC2 %in% spc2gene,]
        rnkTblTmp$rnks[match(rowGenes,rnkTblTmp$SPC1)]
      })
      M[colGenes,rowGenes] <- sapply(rowGenes,function(spc1gene){
        rnkTblTmp <- rnkTbl[rnkTbl$SPC1 %in% spc1gene,]
        rnkTblTmp$rnksT[match(colGenes,rnkTblTmp$SPC2)]
      })
    }
  }
  return(M)
}

orthoPairRnks <- function(geneIDs, rnkTbl, transposed=F){
  
  # make matrix for ranks
  rowGenes <- geneIDs[geneIDs %in% rnkTbl$SPC1]
  colGenes <- geneIDs[geneIDs %in% rnkTbl$SPC2]
  
  sapply(colGenes,function(spc2gene){
    rnkTblTmp <- rnkTbl[rnkTbl$SPC2 %in% spc2gene,]
    if(transposed){
      rnkTblTmp$rnksT[match(rowGenes,rnkTblTmp$SPC1)]
    } else {
      rnkTblTmp$rnks[match(rowGenes,rnkTblTmp$SPC1)]
    }
  }) -> rankMat
  colnames(rankMat) <- colGenes
  rownames(rankMat) <- rowGenes
  
  return(rankMat)
}


getCCSMat <- function(p4d){
  geneIDs <- tipLabels(p4d) # !! order of the tips in the plot doesnt match
  
  M <- matrix(NA_real_,ncol=nTips(p4d),nrow=nTips(p4d),dimnames=list(geneIDs,geneIDs))
  
  for(ccs in CCS){
    m <- ccs[geneIDs[geneIDs %in% rownames(ccs)],geneIDs[geneIDs %in% colnames(ccs)]]
    M[rownames(m),colnames(m)] <- m
    M[colnames(m),rownames(m)] <- t(m)
  }
  return(M)
}

### load CCS #######

# CCS_old <- readRDS("data/subsets/treesWithAtOs/AtOs_CCS.RDS")

myReadRDS <- function(filename){
  cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),"Reading",filename,"...\n")
  readRDS(file = filename)
}

# takes about 35mins to load and uses 130GB of RAM
system.time({
  CCS <-
    dir("data/CCS",pattern="CCS",full.names = T) %>% 
    set_names( sub(".*/(.{4})_.*","\\1",.)) %>% 
    map( myReadRDS )
})

## does not belong here ####
#  
# ### load 1:1 ref.orthos
# refOrthos <- 
#   dir("data/CCS",pattern="11_geneIDs.RDS",full.names = T) %>% 
#   set_names( sub(".*/(..)_(..)(..)11_geneIDs.RDS","\\2\\3_\\1",.)) %>% 
#   map(readRDS) %>% 
#   split(substr(names(.),1,4)) %>%   # split by species pair
#   map( ~ setNames(.x, substr(names(.x),6,7) ) ) %>% 
#   map( ~ as_data_frame(.x) )  # make data_frame of ortholog-pairs for each species pair
# 
# rnks <- map2( .x = CCS, .y = refOrthos, 
#               ~ quickRanks(CCS = .x,spc1genes = .y[[1]], spc2genes = .y[[2]]))
# 
# rnks <- list()
# rnksT <- list()
# for( spcsPair in names(CCS)){
#   cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),"Calculating ranks for",spcsPair,"...\n")
#   rnks[[spcsPair]] <- 
#     quickRanks(CCS = CCS[[spcsPair]],
#                spc1genes = refOrthos[[spcsPair]][[1]],
#                spc2genes = refOrthos[[spcsPair]][[2]])
#   rnksT[[spcsPair]] <- 
#     quickRanksT(CCS = CCS[[spcsPair]],
#                spc1genes = refOrthos[[spcsPair]][[2]],
#                spc2genes = refOrthos[[spcsPair]][[1]])
# }
# 
# spcsPairColors <- setNames(rainbow(10),names(CCS))
# newRanksPlot()
# for( spcsPair in names(rnks)){
#   drawRanksCurve(rnks[[spcsPair]],col=spcsPairColors[spcsPair])
#   drawRanksCurve(rnksT[[spcsPair]],col=spcsPairColors[spcsPair],lty=2)
# }
# barplot(rbind(map_dbl(rnks, median),map_dbl(rnksT, median)),beside = T,las=2)
# legend("bottomright",legend = names(rnks),lty=1,col=spcsPairColors[names(rnks)])
# 
# plot(map_int(refOrthos[names(rnks)], nrow),map_dbl(rnks, median),col=spcsPairColors,pch=20)
# points(map_int(refOrthos[names(rnks)], nrow),map_dbl(rnksT, median),col=spcsPairColors,pch=21)
# text(map_int(refOrthos[names(rnks)], nrow),map_dbl(rnks, median),labels=names(rnks))
# 
# 
# # get 1:1:1:1:1 orthologs
# idx <- which(rowSums(taxCount[,as.character(spc2taxid)]==1) == 5)
# 
# orthoSingletons <-
#   treeData[names(idx)] %>% 
#   map_df( ~ .x$tip.label[ match(spc2taxid,.x$tip.data$taxonID )]) %>% 
#   t() %>% 
#   as.data.frame(stringsAsFactors=F) %>%
#   setNames(names(spc2taxid))
#   
# 
# dim(orthoSingletons)
# 
# rnks2 <- list()
# rnksT2 <- list()
# for( spcsPair in names(CCS)){
#   cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),"Calculating ranks for",spcsPair,"...\n")
#   spc1 <- substr(spcsPair,1,2)
#   spc2 <- substr(spcsPair,3,4)
#   rnks2[[spcsPair]] <- 
#     quickRanks(CCS = CCS[[spcsPair]],
#                spc1genes = orthoSingletons[[spc1]],
#                spc2genes = orthoSingletons[[spc2]])
#   rnksT2[[spcsPair]] <- 
#     quickRanksT(CCS = CCS[[spcsPair]],
#                 spc1genes = orthoSingletons[[spc2]],
#                 spc2genes = orthoSingletons[[spc1]])
# }
# 
# newRanksPlot()
# for( spcsPair in names(rnks)){
#   drawRanksCurve(rnks2[[spcsPair]],col=spcsPairColors[spcsPair])
#   drawRanksCurve(rnksT2[[spcsPair]],col=spcsPairColors[spcsPair],lty=2)
# }
# plot(map_dbl(rnks2, median),map_dbl(rnks, median))
# 
# newRanksPlot()
# for( spcsPair in names(rnks)){
#   drawRanksCurve(rnks[[spcsPair]],col=spcsPairColors[spcsPair])
#   drawRanksCurve(rnks[[spcsPair]][orthoSingletons[[substr(spcsPair,1,2)]]],col=spcsPairColors[spcsPair],lty=2)
# }
# drawRanksCurve(rnks$OsZm)
# drawRanksCurve(rnks$OsZm[orthoSingletons$Os],col="red")
# 
# max(rnks$AtOs)
# plot(-log10(1.001-rnks2$AtZm),-log10(1.001-rnks2$AtSl))
# 
# x <- rnks2[!grepl("Z|O",names(rnks))] %>% reduce( cbind ) %>% rowMeans
# y <- rnks2[grepl("At",names(rnks))] %>% reduce( cbind ) %>% rowMeans
# plot(-log10(1.001-x),-log10(1.001-rnks2$OsZm))
# 


### Candidate tree ###################
# find some interresting candidate trees
# Find a tree with several 1:1 genes

AtOs11orthos <- lapply(
  c(At="data/CCS/At_AtOs11_geneIDs.RDS",
    Os="data/CCS/Os_AtOs11_geneIDs.RDS"),
  readRDS)

treesWithSeveral11 <-
  as_tibble(table(geneID2treeID$At[AtOs11orthos$At])) %>% 
  filter(n>1) %>% .$Var1


# alternative smaller tree for testing
# treeID <- which(rowSums(taxCount) < 15 &
#              rowSums(taxCount[,as.character(spc2taxid)]>0) > 2)[1]


# CBF genes
treeID <- geneID2treeID$At["AT1G63030"]


# select tree of interrest
treeID <- sample(treesWithSeveral11,1)
length(treeData[[treeID]]$tip.label)

# convert to p4d
p4d <- treeDataToPhylo4d(fixTreeData(treeData[[treeID]]))
# remove tips from species without data
p4d <- p4d[tipData(p4d)$taxonID %in% spc2taxid]

# get the phylo version of the tree
tr <- getPhylo(p4d)

plot.phylo( nameTree(phylo =  tr))
nodelabels(nodeData(p4d)$taxonName,frame = "none",adj=c(0,0.5),cex=0.6)
markDupNodes(p4d)

## Extract rank matrix ####
#


# load ranks
orthoRankTables <-
  dir("data/ranksMR", pattern="...._ranks.RDS",full.names=T)  %>%
  set_names(sub(".*/(....)_ranks.RDS","\\1", .)) %>% 
  map(readRDS) %>% 
  map( function( orthoTable ){
    # anonymize the spc1 and spc2 columns in the table
    names(orthoTable)[1:2] <- c("SPC1", "SPC2")
    
    return(orthoTable)
  })



rnksMat <- getRankMat(p4d)

rnksMatLog <- -log10(1.001-rnksMat)


# and CCS matrix
ccsMat <- getCCSMat(p4d)




#### load metadata #############
#
#
source("R/PODCfiles.R")
metaOs <- loadSampleMeta("Os")
metaAt <- loadSampleMeta("At")
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


expAt <- readRDS("data/expMat/PODC_At.RDS")
expOs <- readRDS("data/expMat/PODC_Os.RDS")


#### Analyse tree #########


#
# plot tree
#

plot.phylo( nameTree(phylo =  tr))
nodelabels(nodeData(p4d)$taxonName,frame = "none",adj=c(0,0.5),cex=0.6)
markDupNodes(p4d)


#
# plot ranks/CCS matrix
#

ggplotlyHeatmap(nameMat(rnksMatLog))

ggplotlyHeatmap(nameMat(ccsMat))

#
# plot expression correlation
#

AtGenes <- tipLabels(p4d)[tipData(p4d)$taxonID == spc2taxid["At"]]

plot(as.data.frame(t(expAt[AtGenes, ])),col=metaAt$color,pch=20,cex=0.5)

OsGenes <- tipLabels(p4d)[tipData(p4d)$taxonID == spc2taxid["Os"]]

plot(as.data.frame(t(expOs[OsGenes, ])),col=metaOs$color,pch=20,cex=0.5)
