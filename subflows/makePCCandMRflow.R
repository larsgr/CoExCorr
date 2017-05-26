library(tidyverse)
source("R/flowrUtils.R")

spc2expMatFiles <- c(
  At = "data/expMat/PODC_At.RDS",
  Os = "data/expMat/PODC_Os.RDS",
  Sl = "data/expMat/PODC_Sl.RDS",
  Gm = "data/expMat/EBI_Gm.RDS",
  Zm = "data/expMat/EBI_Zm.RDS")

makePCCandMRjob <- function(spc,memReq,cores,bmOutFile){
  Rjob(source = "R/bigcor.R",
       fun="bigPCCandMRjob",
       jobName = paste0(spc,"_PCCMR"),
       cpu_reserved = cores,
       memory_reserved = memReq,
       paramMat = data.frame( expMatFile=spc2expMatFiles[spc],
                              bmOutFile=bmOutFile,
                              cores=cores))
}


makeCCSandRanksJob <- function(spc1, spc2, memReq, cores,
                               bmFile1,bmFile2,rnksOutfile){
  
  # bigCCSandRanksJob(spc1,spc2,bmFile1,bmFile2,rnksOutfile,cores)
  Rjob(source = "Rjobs/PCCandMRjob.R",
       fun="bigCCSandRanksJob",
       jobName = paste0(spc1,spc2,"_CCSRanks"),
       # start after the MR matrix is ready
       prev_jobs = paste0(c(spc1,spc2),"_PCCMR",collapse=","),
       cpu_reserved = cores,
       memory_reserved = memReq,
       paramMat = data.frame( spc1=spc1,
                              spc2=spc2,
                              bmFile1=bmFile1,
                              bmFile2=bmFile2,
                              rnksOutfile=rnksOutfile,
                              cores=cores))
  
}



makePCC_MR_CCS_RanksFlow <- function(PCCMRFlowDef, CCSRanksFlowDef, outDirMR, outDirRnks){
  flowbind(
    
    
    PCCMRFlowDef %>%
      mutate( bmOutFile = file.path(outDirMR,paste0(spc,"_PCCMR.bin"))) %>% 
      pmap(makePCCandMRjob) %>% # create job for each species
      reduce( flowbind ),
  
    # makeCCSandRanksJob <- function(spc1, spc2, memReq, cores,
    #                                bmFile1,bmFile2,rnksOutfile)
    CCSRanksFlowDef %>%
      mutate( bmFile1 = file.path(outDirMR,paste0(spc1,"_PCCMR.bin.desc")),
              bmFile2 = file.path(outDirMR,paste0(spc2,"_PCCMR.bin.desc")),
              rnksOutfile = file.path(outDirRnks,paste0(spc1,spc2,"_ranks.RDS"))) %>% 
      pmap(makeCCSandRanksJob) %>% # create job for each species
      reduce( flowbind )  # combine to single flowlist
  ) 
}

# # approximate memory requirements
# nGenes <- map_dbl(spc2expMatFiles, ~ nrow(readRDS(.x)))
# # 4x the size of the co-expression matrix
# (nGenes^2*8*4)/2^30
# 
# expMatSize <- map_int(spc2expMatFiles, ~ length(readRDS(.x)))
# 
# CCSmemEstimate <- function(spc1,spc2,cores=1){
#   # big.matrix data:
#   #   size of co-expression matrix for both species plus the CCS matrix :
#   Gb <- (nGenes[spc1]^2+nGenes[spc2]^2+nGenes[spc2]*nGenes[spc1])*8/2^30+
#     # Sum of memory used by each process:
#     (nGenes[spc1]^2+nGenes[spc2]^2+nGenes[spc1]*nGenes[spc2]/cores)*8/2^30
#   return(paste0(ceiling(Gb),"G"))
# }


# spcPairs %>% 
#   mutate( memReq = pmap_chr(., CCSmemEstimate),
#           cores=20 ) %>% 
#   dput_tribble()
# 
# tribble(
#   ~spc, ~memReq, ~cores, ~maxRSS,
#   "At",   "30G",     20, 56,
#   "Os",   "60G",     20, 61,
#   "Sl",   "34G",     10, 38,
#   "Gm",   "86G",     30, 114,
#   "Zm",   "50G",     20, 120 ) %>% 
#   mutate(memBm = (nGenes[spc]^2*8*2)/2^30,
#          overUsePerCPU = (maxRSS-memBm)/cores,
#          xSize = expMatSize[spc]*8/2^30,
#          memEst = memBm+xSize*cores*10) %>% 
#   ggplot( aes(x=xSize,y=overUsePerCPU, label=spc)) + geom_label()
