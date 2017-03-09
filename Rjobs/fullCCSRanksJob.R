source("R/calcRanks.R")

calcRanksJob <- function(spc1,spc2){
  # load CCS
  CCS <- readRDS(paste0("data/CCS/",spc1,spc2,"_CCS.RDS"))
  
  # get ref.orthos
  refOrthos1 <- readRDS(paste0("data/CCS/",spc1,"_",spc1,spc2,"11_geneIDs.RDS"))
  refOrthos2 <- readRDS(paste0("data/CCS/",spc2,"_",spc1,spc2,"11_geneIDs.RDS"))
  
  rnks <- 
    list( quickRanks(CCS,refOrthos1,refOrthos2),
          quickRanksT(CCS,refOrthos2,refOrthos1) )
  
  saveRDS(rnks, paste0("data/CCS/",spc1,spc2,"11_rnks.RDS"))
}