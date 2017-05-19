
sp2expMatFiles <- c(
  At = "data/expMat/PODC_At.RDS",
  Os = "data/expMat/PODC_Os.RDS",
  Sl = "data/expMat/PODC_Sl.RDS",
  Gm = "data/expMat/EBI_Gm.RDS",
  Zm = "data/expMat/EBI_Zm.RDS")

# get the At1Os1 geneIDs that also have expresion
get11GeneIds <- function(spc1,spc2,outDir){
  source("R/orthoUtils.R")
  dir.create(outDir,showWarnings = F,recursive = T)
  
  cat("Reading",spc1,spc2,"orthologs...\n")
  
  orthos <- loadOrthoTable(spc1,spc2)
  
  # get 1:1 orthos
  orthos11 <- orthos[orthos$otype == "1:1",]
    
  # filter only orthos with expression in both species
  for(spc in c(spc1,spc2)){
    # load expression matrix
    expMatFile <- sp2expMatFiles[spc]
    cat("Reading",expMatFile,"...\n")
    expMat <- readRDS(expMatFile)
    orthos11 <- orthos11[orthos11[[spc]] %in% rownames(expMat), ]
  }
  
  # write geneIDs to file
  for(spc in c(spc1,spc2)){
    geneIDs <- orthos11[[spc]]
    outFile <- file.path(outDir,paste0(spc,"_11_geneIDs.RDS"))
    cat("Writing",outFile,"...\n")
    saveRDS(geneIDs,outFile)
  }
}


calcCor <- function( spc, geneIDsubsetFile, outFile, method){
  library(BSplineMI)
  expMat <- readRDS(sp2expMatFiles[spc])
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
  # use absolute correlation
  M <- abs(M)
  
  # quick fix for NA values
  M[is.na(M)] <- 0
  
  # scale to get a value between 0..1 (needed for MI)
  if(max(M)>1)
    M <- M/max(M)
  
  # Use exponent 6 as soft threshold (as advised on WGCNA webpage)
  TOM <- TOMsimilarity(M^6)
  dimnames(TOM) <- dimnames(M)
  diag(TOM) <- 0
  return(TOM)
}

removeNAs <- function(M){
  # quick fix for NA values
  M[is.na(M)] <- 0
  
  return(M)
}


calcCLR_diag0 <- function(M){
  diag(M) <- NA
  calcCLR(M)
}

myLog <- function(...){
  cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),...)
}

# calculate rank scores using different CLR methods
calcScore <- function(M1file,M2file,outFile){
  library(WGCNA)
  source("R/CLR.R")
  source("R/calcRanks.R")

  Mfiles <- c(M1file,M2file)
  spcs <- sub("(..).*","\\1",basename(Mfiles))
  names(Mfiles) <- spcs
  
  myLog("Loading",Mfiles,"...\n")
  M <- lapply( Mfiles, readRDS )
  
  lapply( c(CLR=calcCLR_diag0,MR=calcMutualRank,TOM=calcTOM, none=removeNAs),
          function(CLRmethod){
    myLog("Calculating CLR..\n")
    CLR <- lapply(M,CLRmethod)
    myLog("Calculating CCS..\n")
    CCS <- cor(CLR[[1]],CLR[[2]])
    myLog("Calculating Ranks..\n")
    rnks <- list(
      rnk = quickRanks(CCS, spc1genes = rownames(CCS), spc2genes = colnames(CCS)),  
      rnkT = quickRanksT(CCS, spc1genes = colnames(CCS), spc2genes = rownames(CCS))
    )
    return(rnks)
  }) -> rnksAll # each correlation method
  
  saveRDS(rnksAll,outFile)
  
  # if(any(sapply(rnksAll,methods::is,"try-error"))){
  #   stop("Error inside mclapply. See ",outFile," for error message.")
  # }
}


