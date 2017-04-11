source("R/flowrUtils.R") # Rjob(), flowbind(), subFlow()
source("Rjobs/calcMI.R") # makeMI_flow()
library(tidyverse)

makeCompareMethodsFlow <- function(outDir = "data/subsets/compareMethods"){
  spc <- c("At","Os")
  
  addDepsToFlowList <- function(fl, prev_jobs, dep_type="gather"){
    fl$flowdef$prev_jobs[1] <- prev_jobs
    fl$flowdef$dep_type[1] <- dep_type
    return(fl)
  }  
  
  flowbind(
    Rjob(source = "Rjobs/compareMethodsJob.R",
         fun = "getAt1Os1GeneIds",
         paramMat = data.frame( outDir = outDir )),
    # MI jobs:
    list( prefix = spc,
          expMatFile = paste0("data/expMat/PODC_",spc,".RDS"),
          arraySize = 2, # use 2 cpus per MI job
          outDir = outDir,
          geneIDsubsetFile = paste0(outDir,"/",spc,"_At1Os1_geneIDs.RDS")) %>%
      pmap(makeMI_flow) %>% # make MI flowlist for each species
      map2(spc, subFlow) %>%  # add spc names to job names
      map( addDepsToFlowList, prev_jobs="getAt1Os1GeneIds", dep_type="serial") %>% 
      reduce( flowbind ),  # combine to single flowlist
    
    # PCC and SCC jobs
    crossing( method=c("PCC","SCC"), spc=spc) %>% 
      transmute( expMatFile = paste0("data/expMat/PODC_",spc,".RDS"), 
                 geneIDsubsetFile = paste0(outDir,"/",spc,"_At1Os1_geneIDs.RDS"),
                 outFile = paste0(outDir,"/",spc,"_",method,".RDS"), 
                 method = c(PCC="pearson", SCC="spearman")[method]) %>% 
      as.data.frame() %>% 
      Rjob(paramMat = .,
           prev_jobs = "getAt1Os1GeneIds",
           dep_type = "burst",
           sub_type = "scatter",
           source = "Rjobs/compareMethodsJob.R",
           fun = "calcCor")
    
    
  )
}