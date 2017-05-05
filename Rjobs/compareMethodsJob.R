
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


calcCor <- function( expMatFile, geneIDsubsetFile, outFile, method){
  library(BSplineMI)
  expMat <- readRDS(expMatFile)
  geneIDs <- readRDS(geneIDsubsetFile)
  if( method == "MI"){
    M <- calcSplineMI(expMat[geneIDs, ],nBins = 7,splineOrder = 3)
  } else {
    M <- cor(t(expMat[geneIDs, ]),method=method)
  }
  saveRDS(M, outFile)
}


calcMutualRank <- function(M){
  # M is correlation matrix
  N <- ncol(M)
  MR <- log2(N/sqrt(t(apply(-M,1, rank))*apply(-M,2, rank)))
  diag(MR) <- 0
  return(MR)
}

calcTOM <- function(M){
  M <- abs(M)
  if(max(M)>1)
    M <- M/max(M)
  
  TOM <- TOMsimilarity(M^6)
  dimnames(TOM) <- dimnames(M)
  diag(TOM) <- 0
  return(TOM)
}

myLog <- function(...){
  cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),...)
}

# calculate rank scores using different CLR methods
# uses 8 processes
calcScore <- function(M1file,M2file,outFile){
  library(parallel)
  library(WGCNA)
  source("R/loadTriMatrix.R")
  source("R/CLR.R")
  source("R/calcRanks.R")
  source("R/mc_cor.R")
  
  Mfiles <- c(M1file,M2file)
  spcs <- sub("(..).*","\\1",basename(Mfiles))
  names(Mfiles) <- spcs
  
  myLog("Loading",Mfiles,"...\n")
  M <- lapply( Mfiles, readRDS )
  
  mclapply( mc.cores = 4, 
            c(CLR=calcCLR,MR=calcMutualRank,TOM=calcTOM, none=function(x){x}),
            function(CLRmethod){
    myLog("Calculating CLR..\n")
    CLR <- mclapply(M,CLRmethod,mc.cores = 2)
    myLog("Calculating CCS..\n")
    CCS <- mc_cor2(CLR$At,CLR$Os,2)
    myLog("Calculating Ranks..\n")
    rnks <- list(
      AtOs = quickRanks(CCS, spc1genes = rownames(CCS), spc2genes = colnames(CCS)),  
      OsAt = quickRanksT(CCS, spc1genes = colnames(CCS), spc2genes = rownames(CCS))
    )
    return(rnks)
  }) -> rnksAll # each correlation method
  
  saveRDS(rnksAll,outFile)
  
  if(any(sapply(rnksAll,methods::is,"try-error"))){
    stop("Error inside mclapply. See ",outFile," for error message.")
  }
}


