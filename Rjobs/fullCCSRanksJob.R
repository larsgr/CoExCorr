source("R/calcRanks.R")
library(readr)

myLog <- function(...){
  cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),...)
}

# # DEPRECATED.. use results from allOrthoRanks instead.
# calcRanksJob <- function(spc1,spc2){
#   # load CCS
#   CCS <- readRDS(paste0("data/CCS/",spc1,spc2,"_CCS.RDS"))
#   
#   # get ref.orthos
#   refOrthos1 <- readRDS(paste0("data/CCS/",spc1,"_",spc1,spc2,"11_geneIDs.RDS"))
#   refOrthos2 <- readRDS(paste0("data/CCS/",spc2,"_",spc1,spc2,"11_geneIDs.RDS"))
#   
#   rnks <- 
#     list( quickRanks(CCS,refOrthos1,refOrthos2),
#           quickRanksT(CCS,refOrthos2,refOrthos1) )
#   
#   saveRDS(rnks, paste0("data/CCS/",spc1,spc2,"11_rnks.RDS"))
# }


allOrthoRanks <- function(CCSfile,orthoTblFile,outFile){
  
  # create dir in case it hasn't been created yet
  dir.create(dirname(outFile),recursive = T,showWarnings = F)
  
  myLog("loading",CCSfile,"...\n")
  CCS <- readRDS(CCSfile)
  myLog("loading",orthoTblFile,"...\n")
  orthoTable <- read_tsv(orthoTblFile, col_types = cols(), progress = F)
  
  # remove geneIDs that is not in the CCS matrix
  inCCS <- 
    orthoTable[[1]] %in% rownames(CCS) &
    orthoTable[[2]] %in% colnames(CCS)
  
  myLog(sum(!inCCS),"of",length(inCCS),"ortholog pairs are missing from CCS...\n")

  orthoTable$rnks <- NA_real_
  orthoTable$rnksT <- NA_real_
  
  myLog("calculating ranks...\n")
  orthoTable$rnks[inCCS] <- quickRanks(CCS, orthoTable[[1]][inCCS], orthoTable[[2]][inCCS])
  myLog("calculating transposed ranks...\n")
  orthoTable$rnksT[inCCS] <- quickRanksT(CCS, orthoTable[[2]][inCCS], orthoTable[[1]][inCCS])

  myLog("writing",outFile,"...\n")
  saveRDS(orthoTable,outFile)
}