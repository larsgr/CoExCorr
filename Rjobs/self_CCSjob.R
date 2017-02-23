####
#
# calculate CCS for all genes using 1:1 genes as ref.orthos
#
self_CCSjob <- function(miFile,
                        geneIDfile,
                        refOrthosFile,
                        outFile, cores){
  source("R/CLR.R")
  source("R/mc_cor.R")
  source("R/loadTriMatrix.R")
  
  myLog <- function(...){
    cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),...)
  }
  
  myLog("Started self_CCSjob\n")
  
  # load MI matrix, calculate CLR and filter by ref.orthos
  geneIDs <- readRDS(geneIDfile)
  myLog("Reading",miFile,"...\n")
  mi <- loadTriMatrix(geneIDs, miFile)
  myLog("calcCLR...\n")
  clr <- calcCLR(mi)
  
  # clear some memory
  rm(mi)
  gc()
  
  myLog("Reading", refOrthosFile, "and filtering...\n")
  clr[readRDS(refOrthosFile), ]

  gc() # garbage collect
  
  myLog("calcCCS...\n")
  CCS <- mc_cor2(clr,
                 clr,
                 cores = cores)
  
  myLog("writing",outFile,"...\n")
  saveRDS(CCS, outFile)
  
}