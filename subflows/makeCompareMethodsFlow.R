source("R/flowrUtils.R") # Rjob(), flowbind(), subFlow()
library(tidyverse)


  
makeCompareMethodsFlow <- function(spc1="At",spc2="Os", 
                                   outDir = "data/subsets/compareMethods",
                                   memory_reserved = "6G"){
  spc <- c(spc1,spc2)
  
  methodSuffix <- c(MI="_MI.RDS",PCC="_PCC.RDS",SCC="_SCC.RDS")
  
  
  addDepsToFlowList <- function(fl, prev_jobs, dep_type="gather"){
    fl$flowdef$prev_jobs[1] <- prev_jobs
    fl$flowdef$dep_type[1] <- dep_type
    return(fl)
  }  
  
  flowbind(
    Rjob(source = "Rjobs/compareMethodsJob.R",
         fun = "get11GeneIds",
         paramMat = data.frame( spc1 = spc[1], spc2 = spc[2], outDir = outDir )),
    
    # PCC, SCC and MI jobs (calcCor)
    crossing( method=c("PCC","SCC","MI"), spc=spc) %>% 
      transmute( spc = spc, 
                 geneIDsubsetFile = paste0(outDir,"/",spc,"_11_geneIDs.RDS"),
                 outFile = paste0(outDir,"/",spc,"_",method,".RDS"), 
                 method = c(PCC="pearson", SCC="spearman",MI="MI")[method]) %>% 
      as.data.frame() %>% 
      Rjob(paramMat = .,
           prev_jobs = "get11GeneIds",
           source = "Rjobs/compareMethodsJob.R",
           fun = "calcCor"),
    
    # CLR+CCS+ranks job for each co-expression matrix
    Rjob(
      prev_jobs = "calcCor",
      source = "Rjobs/compareMethodsJob.R",
      memory_reserved = memory_reserved,
      fun = "calcScore",
      paramMat = data.frame( M1file = paste0(outDir,"/",spc[1],methodSuffix),
                             M2file = paste0(outDir,"/",spc[2],methodSuffix),
                             outFile = paste0(outDir,"/",spc[1],spc[2],"_rnks",methodSuffix))
    )
    
  )
}

makeCompareMethodsAllFlow <- function(outDir = "data/subsets/compareMethods"){
  #spcs = c("At","Os","Sl","Gm","Zm")
  
  # define memory requirements
  tribble(
    ~spc1,  ~spc2, ~GB,
    "At",   "Gm",   1,
    "At",   "Os",   6,
    "At",   "Sl",   9,
    "At",   "Zm",   3,
    "Gm",   "Os",   1,
    "Gm",   "Sl",   1,
    "Gm",   "Zm",   1,
    "Os",   "Sl",   4,
    "Os",   "Zm",  15,
    "Sl",   "Zm",   3
  ) %>% 
    # set output directory
    mutate( outDir = file.path(outDir,paste0(spc1,spc2))) %>%
    # convert memory requirement to string
    mutate( memory_reserved = paste0(GB,"G"),
            GB = NULL) %>% 
    # generate workflow for each species pair 
    by_row(do.call, what=makeCompareMethodsFlow, .to="fl") %>% 
    # make job names unique
    with( pmap(list(flowList=fl,prefix=paste0(spc1,spc2,"_")), subFlow) ) %>% 
    # put them into a single flowList
    do.call(what=flowbind)
    
  
  # # generate all combinations of species
  # crossing( spc1=spcs, spc2=spcs) %>% 
  #   filter(spc1 < spc2) %>% 
  #   # set output directory
  #   mutate( outDir = file.path(outDir,paste0(spc1,spc2))) %>%
  #   # generate workflow for each species pair 
  #   by_row(do.call, what=makeCompareMethodsFlow, .to="fl") %>% 
  #   # make job names unique
  #   with( pmap(list(flowList=fl,prefix=paste0(spc1,spc2,"_")), subFlow) ) %>% 
  #   # put them into a single flowList
  #   do.call(what=flowbind)
}
