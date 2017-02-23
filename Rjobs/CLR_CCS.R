####
#
# calculate CCS for all genes using 1:1 genes as ref.orthos
#
MI_CLR_CCSjob <- function(mi_file1, mi_file2,
                          geneIDfile1, geneIDfile2,
                          refOrthosFile1, refOrthosFile2,
                          outFile, cores){
  source("R/CLR.R")
  source("R/mc_cor.R")
  source("R/loadTriMatrix.R")
  
  myLog <- function(...){
    cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),...)
  }

  # load MI matrix, calculate CLR and filter by ref.orthos
  mapply(
    miFile=c(mi_file1, mi_file2),
    geneIDfile=c(geneIDfile1, geneIDfile2),
    refOrthosFile=c(refOrthosFile1,refOrthosFile2),
    FUN = function(miFile,geneIDfile,refOrthosFile){
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
    }
  ) -> CLR

  gc() # garbage collect
  
  myLog("calcCCS...\n")
  CCS <- mc_cor2(CLR[[1]],
                 CLR[[2]],
                 cores = cores)

  myLog("writing",outFile,"...\n")
  saveRDS(CCS, outFile)

}


get11orthos <- function(orthoFile,
                        geneIDfile1,
                        geneIDfile2,
                        refOrthosFile1,
                        refOrthosFile2){
  
  orthos <- read.table(orthoFile,sep="\t",header=T, stringsAsFactors = F)
  orthos <- orthos[ orthos$otype == "1:1", ]
  orthos <- orthos[ orthos[[1]] %in% readRDS(geneIDfile1), ]
  orthos <- orthos[ orthos[[2]] %in% readRDS(geneIDfile2), ]
  
  saveRDS(orthos[[1]],refOrthosFile1)
  saveRDS(orthos[[2]],refOrthosFile2)
}

