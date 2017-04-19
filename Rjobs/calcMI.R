# MI calculation

library(purrr)
library(BSplineMI) # from devtools::install_github("larsgr/BSplineMI")

# write expression matrix in the format used for mi calculation
writeTempExpMat <- function(expMat, filename){
  cat("writing",filename,"...\n")
  saveRDS(expMat,file = filename)  
}


#### job step 1
#
# Prepare genesubset
#
# Load expression matrix
# (optional) extract subset of matrix with given geneIDs
# (alternatively) generate geneID file "{prefix}_geneIDs.RDS"
# Save expression matrix for MI calculation "{prefix}_expmat.RDS"
#
prepMI <- function(expMatFile, geneIDsubsetFile = "", outDir, prefix){
  
  # create the output path (if not already exists)
  dir.create(outDir,showWarnings = F)
  
  # load expression matrix
  expMat <- readRDS(expMatFile)
  
  if(!(geneIDsubsetFile=="")){
    # load geneIDs to subset
    geneIDs <- readRDS(geneIDsubsetFile)
    
    # remove geneIDs that are not in the expression matrix
    geneIDs <- geneIDs[geneIDs %in% rownames(expMat)]
    if(length(geneIDs) == 0){
      stop("prepMI(): geneIDs (",geneIDsubsetFile,")does not match rownames in expMat(",expMatFile,")")
    }
  } else {
    # keep all genes.. 
    geneIDs <- rownames(expMat)

    # save geneIDs as they are not stored in the MI matrix file
    geneIDsFile <- file.path(outDir,paste0(prefix,"_geneIDs.RDS"))
    cat("Writing",geneIDsFile,"...\n")
    saveRDS(geneIDs,geneIDsFile)
  }

  # save expression matrix subset
  tmpExpMatFile <- file.path(outDir,paste0(prefix,".expmat.RDS"))
  writeTempExpMat(expMat[geneIDs, ],tmpExpMatFile)
}

multiMI <- function(tmpExpMatFile, arrayIdx, arraySize){
  cat("Reading",tmpExpMatFile,"...\n")
  x <- readRDS(tmpExpMatFile)
  ngenes <- nrow(x)
  nexp <- ncol(x)
  cat("  genes:",ngenes,"  samples:",nexp,"\n")
  
  # Divide the load equally among the processes (arraySize = No. processes)
  #
  # calcSplineMItoFile generates a triangular matrix so we have to calculate the number of
  # lines (i.e. genes) each process shall generate so that the number of 
  # operations (MI values) per process are approximately equal.
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
  
  fromgene <- 1+c(0, togene[-arraySize])
  
  miFile <- sprintf("%s.mi.%03i",tmpExpMatFile,arrayIdx)
  
  cat("  job part:",arrayIdx,"/",arraySize,"\n")
  cat("  genes:",fromgene[arrayIdx],"-",togene[arrayIdx],"\n")
  
  
  cat("Calculating MI... writing to",miFile,"\n")
  calcSplineMItoFile(x, nBins = 7, splineOrder = 3, 
                     filename = miFile,
                     fromRow = fromgene[arrayIdx],
                     toRow = togene[arrayIdx])
  
}

#### job step 2
#
# calc MI in parallel
#
calcMI <- function(arrayIdx, arraySize, prefix, outDir){
  tmpExpMatFile <- file.path(outDir,paste0(prefix,".expmat.RDS"))
  multiMI(tmpExpMatFile, arrayIdx, arraySize)
}

#### job step 3
#
# merge MI matrix
# clean up temp files
#
mergeMI <- function(prefix, outDir){
  tmpExpMatFile <- file.path(outDir,paste0(prefix,".expmat.RDS"))
  
  miFile <- file.path(outDir,paste0(prefix,".mi"))
  
  # list mi partial files
  miParts <- dir(dirname(tmpExpMatFile), full.names = T,
                 pattern = paste0(basename(tmpExpMatFile),"\\.mi\\.[0-9]{3}$"))
  
  cat("merging mi files...\n")
  cmd <- paste0("cat ",paste(miParts,collapse = " ")," > ",miFile)
  
  cat("CMD:",cmd,"\n")
  system(cmd)
  
  # remove temp. files
  unlink(miParts)
  unlink(tmpExpMatFile)
  
}



#
# generate flow for calculating mi file "[outDir]/[prefix].mi"
#
makeMI_flow <- function(expMatFile,
                        geneIDsubsetFile,
                        outDir,
                        prefix,
                        arraySize=20,
                        mem="4G"){
  subFlow(prefix = prefix, flowList = flowbind(
    # prepMI(expMatFile, geneIDsubsetFile, outDir, prefix)
    Rjob(source = "Rjobs/calcMI.R",
         fun = "prepMI",
         paramMat = data.frame(expMatFile, geneIDsubsetFile, outDir, prefix)),
    
    # calcMI(arrayIdx, arraySize, prefix, outDir)
    Rjob(prev_jobs = "prepMI", dep_type = "burst", sub_type = "scatter",
         source = "Rjobs/calcMI.R",
         fun = "calcMI",
         memory_reserved = mem,
         startSleep = 10, # wait 10 seconds between starting each job
         paramMat = data.frame(arrayIdx=1:arraySize, arraySize,  prefix, outDir)),
    
    # mergeMI(prefix, outDir)
    Rjob( prev_jobs = "calcMI", dep_type = "gather", sub_type = "serial",
          source = "Rjobs/calcMI.R",
          fun = "mergeMI",
          paramMat = data.frame(prefix, outDir))
    
  ))
}

# MI memory requirement calculation (assumes 7 bins):
estMImem <- function(nGenes,nExp){
  # Add 510MB just in...
  paste0((nGenes*nExp*10*8) %/% 2^20 + 510,"M")
}

# Rule of thumb for time estimation:
# minutes = ((genes/1K)^2 * samples/1K) * 2~4
# E.g. 5K genes 0.6K samples => 5^2*0.6 * 2 => 30min (~60min)
# Warning not very accurate...
estMItimeMinutes <- function(nGenes,nExp){
  (nGenes/1000)^2 * (nExp/1000) * 4
}

# Example: estimate time and memory requirements
 
# expMat <- lapply(expMatFiles,readRDS)
# 
# lapply(expMat, function(m){
#   estMImem(nGenes = nrow(m), nExp = ncol(m))
# })
# 
# sapply(expMat, function(m){
#   estMItimeMinutes(nGenes = nrow(m), nExp = ncol(m))
# })


