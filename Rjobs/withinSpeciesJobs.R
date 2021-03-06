source("R/halfStudiesSubset.R")

myLog <- function(...){
  cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),...)
}

withinSpeciesPrepJob <- function(numberOfGenes = 4000,
                                 reps = 10,
                                 outDir = "data/subsets/withinSpecies"){
  # load expression matrices
  expMat <-
    dir("data/expMat",pattern="(EBI|PODC)_..\\.RDS",full.names = T) %>% 
    set_names( sub(".*(EBI|PODC)_(..)\\.RDS","\\2",.) ) %>% 
    map(readRDS)
  
  # get geneIDs for all genes with expression data
  hasExpGeneIDs <-
    map( expMat, rownames ) %>% 
    unlist()
  
  # get subset GeneIDs and refOrthos
  l <- getSubsetGeneIDsAndRefOrthos(numberOfGenes, hasExpGeneIDs)
  

  # generate sample subsets
  sampleSubsets <- generateHalfStudiesSubsets(reps)
  
  # for each species
  for(spc in spcs){
    # create a sub-directory for each species
    spcOutDir <- file.path(outDir,spc)
    dir.create(spcOutDir,recursive = T)
    
    # save gene expression subset
    object.size(expMat[[spc]][l$subsetGeneIDs[[spc]],])
    saveRDS( object = expMat[[spc]][l$subsetGeneIDs[[spc]],],
             file = file.path(spcOutDir,"expMat.RDS"))
    
    # save ref.orthos
    saveRDS( object = l$refOrthos[[spc]],
             file = file.path(spcOutDir,"refOrthos.RDS"))
    
    # save sample subset
    saveRDS( object = sampleSubsets[[spc]],
             file = file.path(spcOutDir,"sampleSubset.RDS"))
  }
}

# For each replicate (2 cpu, 4 GB job)
withinSpeciesJob <- function(spc, repNr, outDir = "data/subsets/withinSpecies", cores = 2){
  
  library(parallel)
  library(BSplineMI)
  source("R/loadTriMatrix.R")
  source("R/CLR.R")
  source("R/calcRanks.R")
  
  spcOutDir <- file.path(outDir,spc)
  
  
  #   Read expression
  expMat <- readRDS( file.path(spcOutDir, "expMat.RDS" ) )
  #   Read sample  subsets
  ss <- readRDS( file.path(spcOutDir, "sampleSubset.RDS" ) )
  
  myLog("Starting parallel MI calculation for",spc,"rep nr.",repNr,"\n")
  
  #  Calculate MI+CLR for each of the two halves: (2 cpu)
  mclapply(c(TRUE, FALSE), mc.cores = max(2,cores), mc.silent = F,
           FUN = function( isFirstHalf ){
    # get sample IDs for this specific rep of subset
    studyIDsubset <- names(which(ss$subsetMat[ ,repNr] == isFirstHalf))
    sampleIDsubset <- ss$sampleID[ ss$studyID %in% studyIDsubset]
    
    # Extract subset from expression and calc MI
    mi <- calcSplineMI(expMat[ ,sampleIDsubset],nBins = 7,splineOrder = 3)
    
    # Calc CLR
    calcCLR(mi)
    
  }) -> CLR
  
  # stop on warnings
  if( length(warnings()) > 0){
    warnings()
    stop("Warning from mclapply")
  } 
  
  myLog("CLR calculation complete\n")
  
  refOrthosAll <- readRDS( file.path(spcOutDir, "refOrthos.RDS" ) )
  
  refOrthos <- refOrthosAll[[1]]
  # For each set of ref.orthos
  mclapply( refOrthosAll, mc.cores = cores, mc.silent = F, FUN = function( refOrthos ){
    
    myLog("Calculating CCS...\n")
    
    # Calc CCS
    CCS <- cor(CLR[[1]][refOrthos, ], CLR[[2]][refOrthos, ])
    
    # Calc Ranks
    rnks <- quickRanks(CCS, spc1genes = rownames(CCS),spc2genes = colnames(CCS))
    rnksT <- quickRanksT(CCS, spc1genes = rownames(CCS),spc2genes = colnames(CCS))
    
    return(list(rnks,rnksT))
  }) -> ranksAll

  # stop on warnings
  if( length(warnings()) > 0){
    warnings()
    stop("Warning from mclapply")
  } 
  
  myLog("CCS calculation complete\n")
  
  #   Save ranks to file
  rnkFile <- file.path(spcOutDir,paste0("rnks_",repNr,".RDS"))
  saveRDS(ranksAll,rnkFile)
}


# For each replicate (2 cpu, 4 GB job)
withinSpeciesPCCMRJob <- function(spc, repNr, outDir = "data/subsets/withinSpecies", cores = 2){
  
  library(parallel)
  source("R/calcRanks.R")
  source("R/calcMR.R")
  
  
  spcOutDir <- file.path(outDir,spc)
  
  
  #   Read expression
  expMat <- readRDS( file.path(spcOutDir, "expMat.RDS" ) )
  #   Read sample  subsets
  ss <- readRDS( file.path(spcOutDir, "sampleSubset.RDS" ) )
  
  myLog("Starting PCC+MR calculation for",spc,"rep nr.",repNr,"\n")
  
  # extract both half expression matrix subsets
  lapply(c(TRUE, FALSE),function( isFirstHalf ){
    # get sample IDs for this specific rep of subset
    studyIDsubset <- names(which(ss$subsetMat[ ,repNr] == isFirstHalf))
    sampleIDsubset <- ss$sampleID[ ss$studyID %in% studyIDsubset]
    return(expMat[ ,sampleIDsubset])
  }) -> expMatSubsets
  
  # Find genes with no variance in either subset
  do.call( "|", lapply(expMatSubsets,function(x){
    apply(x,1,max)==apply(x,1,min) # max=min -> no variance!
  })) -> noVar
  
  #  Calculate PCC+MR for each of the two subset halves: (2 cpu)
  mclapply(expMatSubsets, mc.cores = max(2,cores),FUN = function( x ){
    # filter noVar genes and calculate PCC
    pcc <- WGCNA::cor(t(x[!noVar, ]))
    # Calculate MR
    calcLogMR(pcc)
  }) -> MR

  # stop on warnings
  if( length(warnings()) > 0){
    warnings()
    stop("Warning from mclapply")
  } 
  
  myLog("PCC+MR calculation complete\n")
  
  refOrthosAll <- readRDS( file.path(spcOutDir, "refOrthos.RDS" ) )
  
  # For each set of ref.orthos
  mclapply( refOrthosAll, mc.cores = cores, mc.silent = F, FUN = function( refOrthos ){
    
    myLog("Calculating CCS...\n")
    
    # filter noVar genes from ref.orthologs
    filteredRefOrthos <- refOrthos[!(refOrthos %in% names(noVar)[noVar])]
    
    # Calc CCS
    CCS <- calcCCS(MR[[1]], MR[[2]],filteredRefOrthos, filteredRefOrthos)

    # Calc Ranks
    return(
      list(
        rnks = quickRanks(CCS, spc1genes = rownames(CCS),spc2genes = colnames(CCS)),
        rnksT = quickRanksT(CCS, spc1genes = colnames(CCS),spc2genes = rownames(CCS))
        ))
  }) -> ranksAll
  
  # stop on warnings
  if( length(warnings()) > 0){
    warnings()
    stop("Warning from mclapply")
  } 
  
  myLog("CCS calculation complete\n")
  
  #   Save ranks to file
  rnkFile <- file.path(spcOutDir,paste0("rnks_",repNr,".RDS"))
  saveRDS(ranksAll,rnkFile)
}

