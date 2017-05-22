library(tidyverse) # tribble, pmap
source("R/flowrUtils.R")

makeProcessExpMatFlow <- function(){
  # make individual jobs for each species using either the PODC or EBI procedure
  tribble(
    ~spc, ~fun,
    "Os", "processExpMat_PODC",
    "At", "processExpMat_PODC",
    "Sl", "processExpMat_PODC",
    "Zm", "processExpMat_EBI",
    "Gm", "processExpMat_EBI"
  ) %>% 
    pmap( function(spc, fun){
      Rjob(
        jobName = paste0(spc,"_",fun),
        source = "Rjobs/processExpMat.R",
        fun = fun,
        paramMat = data.frame( spc = spc )
      )    
    }) %>% 
    do.call(what = flowbind)
}
  