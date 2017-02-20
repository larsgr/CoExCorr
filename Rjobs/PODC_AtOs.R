# PODC all At:Os orthologs

# write expression matrix in the format used for mi calculation
writeTempExpMat <- function(expMat, filename){
  cat("writing",filename,"...\n")
  write.table(expMat,file = filename,sep = " ",row.names = F, col.names = F)  
}

source("R/calcRanks.R")

#### job step 1
#
# Prepare genesubset
#
# Get genes that have orthologs in At and Os
# Save expression matrix for MI calculation
#
PODC_AtOs_prepMI <- function(outDir="data/subsets/PODC_AtOs"){
  dir.create(outDir,showWarnings = F)
  
  At_Os_ortho_file <- "indata/At_Os_orthologs.txt"

  cat("Reading",At_Os_ortho_file,"...\n")
  At_Os_ortho <- read.table(At_Os_ortho_file,header=T,as.is = T)
  
  for(spc in c("At","Os")){
    # load expression matrix
    expMatFile <- paste0("data/expMat/PODC_",spc,".RDS")
    cat("Reading",expMatFile,"...\n")
    expMat <- readRDS(expMatFile)
    
    # get geneIDs of genes with orthologs
    geneIDs <- unique(At_Os_ortho[[spc]])
    geneIDs <- geneIDs[geneIDs %in% rownames(expMat)]
    outFile <- file.path(outDir,paste0(spc,"_geneIDs.RDS"))
    cat("Writing",outFile,"...\n")
    saveRDS(geneIDs,outFile)

    
    # save expression matrix subset
    tmpExpMatFile <- file.path(outDir,paste0("expMat_",spc))
    writeTempExpMat(expMat[geneIDs, ],tmpExpMatFile)
  }
}

multiMI <- function(tmpExpMatFile, arrayIdx, arraySize){
  cat("Calculating number of genes and samples in",tmpExpMatFile,"...\n")
  x <- readLines(tmpExpMatFile)
  ngenes <- length(x)
  nexp <- length(strsplit(x[1],split = " ")[[1]])
  rm(x)
  cat("genes:",ngenes,"  samples:",nexp,"\n")
  
  totOps <- ((ngenes-1)*ngenes)/2
  opsPerJob <- ceiling(totOps/arraySize)
  
  ops <- 0
  togene <- integer(arraySize)
  idx <- 1
  for( i in 1:ngenes){
    ops <- ops + i - 1 
    if( ops >= opsPerJob ){
      togene[idx] <- i
      ops <- ops - opsPerJob
      idx <- idx + 1
    }
  }
  togene[arraySize] <- ngenes
  
  fromgene <- c(0, togene[-arraySize])
  
  miFile <- sprintf("%s.mi.%03i",tmpExpMatFile,arrayIdx)
  

  cmd <- paste("~/genepair/genepair 1", 
               tmpExpMatFile, # microarray_file.txt
               miFile, # mi_file.txt
               ngenes, # ngenes
               nexp, # nexp
               fromgene[arrayIdx], # fromgene
               togene[arrayIdx], # togene
               7, # num_bins
               3) # spline_order
  
  cat("CMD:",cmd,"\n")
  
  # run genepair to calculate MI
  system(cmd)
}

#### job step 2
#
# calc MI in parallel
#
PODC_AtOs_calcMI <- function(arrayIdx, arraySize, spc, outDir="data/subsets/PODC_AtOs"){
  tmpExpMatFile <- file.path(outDir,paste0("expMat_",spc))
  multiMI(tmpExpMatFile, arrayIdx, arraySize)
}

#### job step 3
#
# merge MI matrix
# clean up temp files
#
PODC_AtOs_mergeMI <- function(spc, outDir="data/subsets/PODC_AtOs"){
  tmpExpMatFile <- file.path(outDir,paste0("expMat_",spc))
  
  miFile <- file.path(outDir,paste0(spc,".mi"))
  
  # list mi partial files
  miParts <- dir(dirname(tmpExpMatFile), full.names = T,
                 pattern = paste0(basename(tmpExpMatFile),"\\.mi\\.[0-9]{3}$"))
  
  cat("merging mi files...\n")
  cmd <- paste0("cat ",paste(miParts,collapse = " ")," > ",miFile)

  cat("CMD:",cmd,"\n")
  system(cmd)
  
  # remove temp. files
  unlink(dir(dirname(tmpExpMatFile), full.names = T,
             pattern = paste0(basename(tmpExpMatFile),"\\.mi\\..+")))
  unlink(tmpExpMatFile)
  
}

####
#
# calculate CCS and ranks using only 1:1 genes
#
PODC_AtOs_11orthos <- function(){
  dataPath <- "data/subsets/PODC_AtOs"
  source("R/CLR.R")
  source("R/LoadComPlExData.R")
  
  
  # Get geneIDs
  lapply(c(At="At",Os="Os"),function(spc){
    readRDS(file.path(dataPath,paste0(spc,"_geneIDs.RDS")))
  }) -> geneIDs
  
  # reading ortholog data
  At_Os_ortho_file <- "indata/At_Os_orthologs.txt"
  At_Os_orthos <- read.table(At_Os_ortho_file,header=T,as.is = T)
  
  At_Os_11orthos <- At_Os_orthos[ At_Os_orthos$otype=="1:1",] 
  # remove genes not in the expression data
  keepThese <- (At_Os_11orthos$At %in% geneIDs$At) &
    (At_Os_11orthos$Os %in% geneIDs$Os)
  At_Os_11orthos <- At_Os_11orthos[keepThese, ]
  # save to file
  saveRDS(At_Os_11orthos,"data/subsets/PODC_AtOs/At_Os_11orthos.RDS")
  
  # load MI matrix but keep only the 1:1 orthos
  lapply(c(At="At",Os="Os"),function(spc){
    miFile <- file.path(dataPath,paste0(spc,".mi"))
    cat("Reading",miFile,"...\n")
    mi <- loadTriMatrix(geneIDs[[spc]], miFile)
    # keep only the 1:1 orthos
    mi <- mi[At_Os_11orthos[[spc]],At_Os_11orthos[[spc]]]
    gc() # collect garbage
    return(mi)
  }) -> MI11
  
  cat("calcCLR...\n")
  CLR <- lapply(MI11,calcCLR)
  cat("calcCCS...\n")
  CCS <- cor(CLR[[1]],CLR[[2]])
  saveRDS(CCS,file.path(dataPath, "CCS11.RDS"))

  rnks <- calcRanks(CCS)
  saveRDS(rnks,file.path(dataPath, "rnks11.RDS"))
}

####
#
# calculate CCS for all genes using 1:1 genes as ref.orthos
#
PODC_AtOs_11refs <- function(){
  dataPath <- "data/subsets/PODC_AtOs"
  source("R/CLR.R")
  source("R/LoadComPlExData.R")
  
  
  # Get geneIDs
  lapply(c(At="At",Os="Os"),function(spc){
    readRDS(file.path(dataPath,paste0(spc,"_geneIDs.RDS")))
  }) -> geneIDs
  
  # reading ortholog data
  At_Os_11orthos <- readRDS("data/subsets/PODC_AtOs/At_Os_11orthos.RDS")
  
  # load MI matrix
  lapply(c(At="At",Os="Os"),function(spc){
    miFile <- file.path(dataPath,paste0(spc,".mi"))
    cat("Reading",miFile,"...\n")
    mi <- loadTriMatrix(geneIDs[[spc]], miFile)
    gc() # collect garbage
    return(mi)
  }) -> MI
  
  cat("calcCLR...\n")
  CLR <- lapply(MI,calcCLR)
  
  # clear some memory
  rm(MI)
  gc()
  
  cat("calcCCS...\n")
  CCS <- cor(CLR$At[At_Os_11orthos$At, ],CLR$Os[At_Os_11orthos$Os, ])
  saveRDS(CCS,file.path(dataPath, "CCS11refs.RDS"))

  # clear some memory
  rm(CLR)
  gc()
  
  cat("calcRanks\n")
  rnks <- calcRanks(CCS[At_Os_11orthos$At,At_Os_11orthos$Os])
  saveRDS(rnks,file.path(dataPath, "rnks11refs.RDS"))
}

makePODC_AtOs_flow <- function(arraySize=20){
  flowbind(
    # PODC_At_HalfSubsetsJob1
    Rjob(source = "Rjobs/PODC_AtOs.R",
         fun = "PODC_AtOs_prepMI" ),
    
    # PODC_AtOs_calcMI_At
    Rjob(jobName = "PODC_AtOs_calcMI_At",
         prev_jobs = "PODC_AtOs_prepMI", dep_type = "burst", sub_type = "scatter",
         source = "Rjobs/PODC_AtOs.R",
         fun = "PODC_AtOs_calcMI",
         paramMat = data.frame(arrayIdx=1:arraySize, arraySize, spc="At")),
    
    # PODC_AtOs_mergeMI_At
    Rjob( jobName = "PODC_AtOs_mergeMI_At",
          prev_jobs = "PODC_AtOs_calcMI_At", dep_type = "gather", sub_type = "serial",
          source = "Rjobs/PODC_AtOs.R",
          fun = "PODC_AtOs_mergeMI",
          paramMat = data.frame(spc="At")),
  
    # PODC_AtOs_calcMI_Os
    Rjob(jobName = "PODC_AtOs_calcMI_Os",
         prev_jobs = "PODC_AtOs_prepMI", dep_type = "burst", sub_type = "scatter",
         source = "Rjobs/PODC_AtOs.R",
         fun = "PODC_AtOs_calcMI",
         paramMat = data.frame(arrayIdx=1:arraySize, arraySize, spc="Os")),
    
    # PODC_AtOs_mergeMI_Os
    Rjob( jobName = "PODC_AtOs_mergeMI_Os",
          prev_jobs = "PODC_AtOs_calcMI_Os", dep_type = "gather", sub_type = "serial",
          source = "Rjobs/PODC_AtOs.R",
          fun = "PODC_AtOs_mergeMI",
          paramMat = data.frame(spc="Os")),
    
    #PODC_AtOs_11orthos
    Rjob( prev_jobs = "PODC_AtOs_mergeMI_Os,PODC_AtOs_mergeMI_At", 
          dep_type = "gather", sub_type = "serial",
          source = "Rjobs/PODC_AtOs.R",
          memory_reserved = "15G",
          fun = "PODC_AtOs_11orthos"),
    
    #PODC_AtOs_11refs
    Rjob( prev_jobs = "PODC_AtOs_11orthos", 
          dep_type = "gather", sub_type = "serial",
          source = "Rjobs/PODC_AtOs.R",
          memory_reserved = "30G",
          fun = "PODC_AtOs_11refs")
    
  )
}

