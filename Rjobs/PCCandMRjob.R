library(bigmemory)
library(parallel)
library(WGCNA)

options(bigmemory.allow.dimnames=TRUE)


myLog <- function(...){
  cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),...)
}


bigPCCandMR <- function(expMat, bmOutFile, cores = 1){
  
  nGenes <- nrow(expMat)
  geneIDs <- rownames(expMat)
  
  # transpose expMat
  expMat <- t(expMat)
  
  
  
  myLog("Creating ranks big.matrix...\n")
  
  # create output big.matrix
  bigRanks <- big.matrix(nrow = nGenes, ncol = nGenes)
  
  
  # split genes into chunks to be processed by each core
  idx <- split(1:nGenes,cut(1:nGenes,cores,labels = FALSE))
  
  # calculate PCC in paralell
  # Note: WGCNA's cor function has a thread parameter but it doesn't seem to 
  # work. Hence the mclapply.
  myLog("Calculating PCC and ranks...\n")
  mclapply(idx,mc.cores = cores,FUN=function(i){
    # calculate pearson correlation
    pcc <- WGCNA::cor(y = expMat[ ,i],x = expMat)
    # calculate rank (reverse)
    bigRanks[ ,i] <- apply(pcc,2,function(x){log2(rank(-x))})
    invisible()
  })
  
  myLog("Creating big.matrix",bmOutFile,"...\n")
  
  # create output big.matrix
  bigMR <- big.matrix( nrow = nGenes, ncol = nGenes, 
                       dimnames = list(geneIDs,geneIDs),
                       backingfile = basename(bmOutFile),
                       backingpath = dirname(bmOutFile),
                       descriptorfile = paste0(basename(bmOutFile),".desc"))
  
  myLog("Calculating mutual ranks...\n")
  
  # # let each process work on different part of the matrix
  idxChunks <- split(1:nGenes,cut(1:nGenes,cores,labels = FALSE))
  mclapply(idxChunks,mc.cores = cores,FUN=function(idxChunk){
    for(i in idxChunk){
      # for(i in 1:nGenes){
      # calculate MR = log2(N) - 0.5*(logR+t(logR))
      bigMR[,i] <- log2(nGenes)- 0.5*(bigRanks[ ,i] + bigRanks[i, ])
    }
  })
  
  myLog("Flushing data to file...\n")
  flush(bigMR)
  
  return()
}


bigPCCandMRjob <- function(expMatFile,bmOutFile,cores){
  
  # create output directory
  dir.create(dirname(bmOutFile),recursive = T,showWarnings = F)
  
  myLog("loading",expMatFile,"...\n")
  myLog("calculate PCC+MR for all genes\n")
  bigPCCandMR( expMat = readRDS(expMatFile),
               bmOutFile = bmOutFile, 
               cores = cores)
  myLog("Done.\n")
}




get11Orthos <- function(spc1,spc2,geneIDs1, geneIDs2){
  # load orthos
  source("R/orthoUtils.R")
  orthos <- loadOrthoTable(spc1,spc2)
  
  # which are 1:1 and in both matrices
  is11 <-
    orthos$otype == "1:1" &
    orthos[[spc1]] %in% geneIDs1 &
    orthos[[spc2]] %in% geneIDs2
  
  return( orthos[is11,c(spc1,spc2)] )
}


# multi-core ranks
mcQuickRanks <- function(CCS,spc1genes,spc2genes,cores,rnkT=F){
  spc1idx <- match(spc1genes,rownames(CCS))
  spc2idx <- match(spc2genes,colnames(CCS))
  
  # split into chunks
  n <- length(spc1idx)
  chunkIdx <- split(1:n,cut(1:n,cores,labels = FALSE))
  mclapply(chunkIdx,mc.cores = cores,FUN = function(i){
    if(rnkT){
      quickRanksT(CCS, spc2idx[i], spc1idx[i])
    } else {
      quickRanks(CCS, spc1idx[i], spc2idx[i])
    }
  }) -> rnks
  
  return(unlist(rnks,use.names = F))
}

allOrthoRanks <- function(CCS,spc1,spc2,outFile,cores){

  source("R/orthoUtils.R")
  source("R/calcRanks.R")
  
  myLog("loading orthologs...\n")
  orthos <- loadOrthoTable(spc1,spc2)

  # remove geneIDs that is not in the CCS matrix
  inCCS <- 
    orthos[[1]] %in% rownames(CCS) &
    orthos[[2]] %in% colnames(CCS)
  
  myLog(sum(!inCCS),"of",length(inCCS),"ortholog pairs are missing from CCS...\n")
  
  orthos$rnks <- NA_real_
  orthos$rnksT <- NA_real_
  
  spc1idx <- match(orthos[[1]][inCCS],rownames(CCS))
  spc2idx <- match(orthos[[2]][inCCS],colnames(CCS))
  
  myLog("calculating ranks...\n")
  orthos$rnks[inCCS] <- quickRanks(CCS, spc1idx, spc2idx)
  # orthos$rnks[inCCS] <- mcQuickRanks(CCS, 
  #                                    spc1genes = orthos[[1]][inCCS], 
  #                                    spc2genes = orthos[[2]][inCCS],
  #                                    cores = cores,rnkT=FALSE)

  myLog("calculating transposed ranks...\n")
  orthos$rnksT[inCCS] <- quickRanksT(CCS, spc2idx, spc1idx)
  # orthos$rnksT[inCCS] <- mcQuickRanks(CCS, 
  #                                     spc1genes = orthos[[1]][inCCS], 
  #                                     spc2genes = orthos[[2]][inCCS],
  #                                     cores = cores, rnkT=TRUE)
  
  myLog("writing",outFile,"...\n")
  
  # create dir in case it hasn't been created yet
  dir.create(dirname(outFile),recursive = T,showWarnings = F)
  
  saveRDS(orthos,outFile)
}

calcBigCCS <- function(bigMs,ref.orthos,cores){
  # number of genes
  nGenes1 <- ncol(bigMs[[1]])
  nGenes2 <- ncol(bigMs[[2]])
  
  # get row indexes for ref.orthos
  refIdx1 <- match(ref.orthos[[1]],rownames(bigMs[[1]]))
  refIdx2 <- match(ref.orthos[[2]],rownames(bigMs[[2]]))
  
  
  # reserve shared memory for the results
  bigCCS <- big.matrix(nrow = nGenes1, ncol = nGenes2)
  
  # split into cores^2 chunks
  idx1 <- split(1:nGenes1,cut(1:nGenes1,cores,labels = FALSE))
  idx2 <- split(1:nGenes2,cut(1:nGenes2,cores,labels = FALSE))
  idxCombo <- expand.grid(i1=1:cores,i2=1:cores)
  
  mcmapply(i1=idxCombo$i1,i2=idxCombo$i2,
           mc.cores = cores, #mc.preschedule = F,
           FUN=function(i1,i2){
    
    myLog("Calculating CCS chunk",i1,"-",i2,"...\n")
    # get submatrices
    m1 <- bigMs[[1]][refIdx1,idx1[[i1]]]
    m2 <- bigMs[[2]][refIdx2,idx2[[i2]]]

    # set self MR equal to mean of column to ignore it from correlations
    rowi <- match(colnames(m1),rownames(m1))
    for(i in which(!is.na(rowi))){
      m1[rowi[i],i] <- mean(m1[-rowi[i],i])
    }
    
    rowi <- match(colnames(m2),rownames(m2))
    for(i in which(!is.na(rowi))){
      m2[rowi[i],i] <- mean(m2[-rowi[i],i])
    }
    
    # calculate pearson correlation
    bigCCS[idx1[[i1]],idx2[[i2]]] <- WGCNA::cor(m1,m2)
  })
  
  dimnames(bigCCS) <- list(rownames(bigMs[[1]]),rownames(bigMs[[2]]))
  return(bigCCS)
}

bigCCSandRanksJob <- function(spc1,spc2,bmFile1,bmFile2,rnksOutfile,cores){

  # for each spc
  #   load (attach) big.matrix co-expression data
  
  bmFiles <- setNames(c(bmFile1,bmFile2),c(spc1,spc2))
  bigMs <- lapply(bmFiles, attach.big.matrix )
  
  # get ref.orthos (1:1 orthos)
  ref.orthos <- get11Orthos(spc1,spc2,
                            geneIDs1 = rownames(bigMs[[spc1]]),
                            geneIDs2 = rownames(bigMs[[spc2]]))
  
  myLog( "calculate CCS...\n")
  bigCCS <- calcBigCCS(bigMs,ref.orthos,cores)
  
  # get ranks
  allOrthoRanks(CCS = bigCCS,spc1, spc2,outFile = rnksOutfile, cores = cores)
}