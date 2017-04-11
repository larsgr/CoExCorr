
# get the At1Os1 geneIDs that also have expresion
getAt1Os1GeneIds <- function(outDir){
  dir.create(outDir,showWarnings = F)
  
  orthoFile <- "indata/orthologs/At_Os_orthologs.txt"
  cat("Reading",orthoFile,"...\n")
  orthos <- read.table(orthoFile,header=T,as.is = T)
  
  # get 1:1 orthos
  orthos11 <- orthos[orthos$otype == "1:1",]
    
  # filter only orthos with expression in both species
  for(spc in c("At","Os")){
    # load expression matrix
    expMatFile <- paste0("data/expMat/PODC_",spc,".RDS")
    cat("Reading",expMatFile,"...\n")
    expMat <- readRDS(expMatFile)
    orthos11 <- orthos11[orthos11[[spc]] %in% rownames(expMat), ]
  }
  
  # write geneIDs to file
  for(spc in c("At","Os")){
    geneIDs <- orthos11[[spc]]
    outFile <- file.path(outDir,paste0(spc,"_At1Os1_geneIDs.RDS"))
    cat("Writing",outFile,"...\n")
    saveRDS(geneIDs,outFile)
  }
}


calcCor <- function( expMatFile, geneIDsubsetFile, outFile, method = c("pearson", "spearman")){
  expMat <- readRDS(expMatFile)
  geneIDs <- readRDS(geneIDsubsetFile)
  M <- cor(t(expMat[geneIDs, ]),method=method)
  saveRDS(M, outFile)
}


calcMutualRank <- function(M){
  # M is correlation matrix
  # rank
  t(apply(M,1, rank))*apply(M,2, rank)
}

# calculate rank scores using different CLR methods
calcScore <- function(M1file,M2file,CLRmethod=c("z-score","rank","none"), outFile){
  source("R/loadTriMatrix.R")
  source("R/CLR.R")
  source("R/calcRanks.R")
  
  # load matrices
  lapply(list( M1file, M2file ), function(Mfile){
    # if filename ends with .mi then load as trimatrix
    if(grepl(".mi$",Mfile)){
      geneIDs <- readRDS(sub(".mi$","_At1Os1_geneIDs.RDS",Mfile))
      return(loadTriMatrix(geneIds, filename = Mfile))
    } else {
      return(readRDS(Mfile))
    }
  }) -> M
  
  if( CLRmethod == "z-score" ){
    cat("CLRmethod: z-score\n")
    CLR <- lapply(M,calcCLR)
  } else if(CLRmethod == "rank"){
    cat("CLRmethod: rank\n")
    CLR <- lapply(M,calcMutualRank)
  } else{
    cat("CLRmethod: none\n")
    CLR <- M
  }
  
  # calc CCS
  CCS <- cor(CLR[[1]],CLR[[2]])
  
  # calc ranks of diagonals, i.e. 1:1 orthos
  rnks <- list(
    quickRanks(CCS, spc1genes = rownames(CCS), spc2genes = colnames(CCS)),  
    quickRanksT(CCS, spc1genes = colnames(CCS), spc2genes = rownames(CCS))
  )
  
  saveRDS(rnks,outFile)
}
