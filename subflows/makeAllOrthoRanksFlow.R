source("R/flowrUtils.R")
library(purrr)
library(tidyr)
library(dplyr)

# make job for running all reps for all spepcies
makeAllOrthoRanksFlow <- function(outDir = "data/ranks",
                                  mem = "20G",
                                  spcs = c("At","Gm","Os","Sl","Zm")){
  
    
  # allOrthoRanks(CCSfile,orthoTblFile,outFile) 
  Rjob( 
    sub_type = "scatter",
    memory_reserved = mem, 
    startSleep = 60, # 1 minute wait between start of each job
    source = "Rjobs/fullCCSRanksJob.R", 
    fun = "allOrthoRanks",
    paramMat = 
      crossing(spc1 = spcs,spc2 = spcs) %>%  # get all possible spc pairs
      filter( spc1 < spc2 ) %>% # remove redundant pairs
      transmute( CCSfile = paste0("data/CCS/",spc1,spc2,"_CCS.RDS"),
                 orthoTblFile = paste0("indata/orthologs/",spc1,"_",spc2,"_orthologs.txt"),
                 outFile = paste0(outDir,"/",spc1,spc2,"_orthoRanks.RDS"))
  )
}

