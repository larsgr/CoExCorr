# Generate workflow for generating all full MI matrices and pairwise CCS
#

library(purrr) # map2(), pmap(), reduce()
library(tidyr) # crossing()
library(dplyr) # transmute()

source("R/flowrUtils.R") # Rjob(), flowbind(), subFlow()
source("Rjobs/calcMI.R") # makeMI_flow()

# memory required by the CCS calculation is approx 5 times the size of the 
# largest MI matrix.
# Gm is the largest matrix ~22GB so let's reserve plenty
makeCCSflow <- function(spc1,spc2,memReq="150G",cpuReq=20){
  flowbind(
    Rjob(
      jobName = "get11refs",
      source = "Rjobs/CLR_CCS.R", 
      fun = "get11orthos",
      paramMat = data.frame( 
        orthoFile = paste0("indata/orthologs/",spc1,"_",spc2,"_orthologs.txt"),
        geneIDfile1 = paste0("data/MI/",spc1,"_geneIDs.RDS"),
        geneIDfile2 = paste0("data/MI/",spc2,"_geneIDs.RDS"),
        refOrthosFile1 = paste0("data/CCS/",spc1,"_",spc1,spc2,"11_geneIDs.RDS"),
        refOrthosFile2 = paste0("data/CCS/",spc2,"_",spc1,spc2,"11_geneIDs.RDS")
      )
    ),
    Rjob(
      jobName = "calcCCS",
      prev_jobs = "get11refs",
      dep_type = "serial",
      source = "Rjobs/CLR_CCS.R",
      fun = "MI_CLR_CCSjob",
      memory_reserved = memReq,
      queue = "hugemem",
      cpu_reserved = cpuReq,
      paramMat = data.frame(
        mi_file1 = paste0("data/MI/",spc1,".mi"),
        mi_file2 = paste0("data/MI/",spc2,".mi"),
        geneIDfile1 = paste0("data/MI/",spc1,"_geneIDs.RDS"),
        geneIDfile2 = paste0("data/MI/",spc2,"_geneIDs.RDS"),
        refOrthosFile1 = paste0("data/CCS/",spc1,"_",spc1,spc2,"11_geneIDs.RDS"),
        refOrthosFile2 = paste0("data/CCS/",spc2,"_",spc1,spc2,"11_geneIDs.RDS"),
        outFile = paste0("data/CCS/",spc1,spc2,"_CCS.RDS") ,
        cores = cpuReq)
    )
  )
}


#
# Generate workflow for generating all full MI matrices
#


makeMIandCCSflow <- function(MIFlowDef){

  #
  # make MI flows for each species
  #
  
  miFlows <- 
    MIFlowDef %>% # use this to make parameters for makeMI_flow()
    transmute( prefix = spc,
               expMatFile = expMatFile,
               arraySize = cpuReq,
               mem = memReq,       
               outDir = "data/MI",
               geneIDsubsetFile = "" ) %>%
    pmap(makeMI_flow) %>% # make MI flowlist for each species
    reduce( flowbind )  # combine to single flowlist
  
  #
  # Make CCS flows for each combination of species
  #
  
  spcs <- MIFlowDef$spc
  
  spcPairs <-
    crossing(spc1 = spcs,spc2 = spcs) %>%  # get all possible spc pairs
    filter( spc1 < spc2 ) # remove redundant pairs
  
  # Helper function that add specified dependency to first job in flowlist
  addDepsToFlowList <- function(fl, prev_jobs, dep_type="gather"){
    fl$flowdef$prev_jobs[1] <- prev_jobs
    fl$flowdef$dep_type[1] <- dep_type
    return(fl)
  }
  
  CCSflows <- 
    pmap(spcPairs, makeCCSflow) %>% 
    map2(pmap_chr(spcPairs,paste0), subFlow)  %>% # add spc names to job names
    map2(paste0(spcPairs$spc1,"mergeMI,", # add dependencies to MI jobs
                spcPairs$spc2,"mergeMI")
         ,addDepsToFlowList) %>% 
    reduce( flowbind )  # combine to single flowlist
  
  return( flowbind(miFlows, CCSflows) )
}

