# Need to load the slurm environment to be able to submit jobs.
# But first need to init the environment modules manager
library(RLinuxModules) # devtools::install_github("larsgr/RLinuxModules")
moduleInit()
module("load slurm")

library(flowr)

# Set the module_cmds option
opts_flow$set(module_cmds = "source /etc/profile.d/modules.sh\nmodule load R/3.3.2")

source("R/flowrUtils.R")
library(tidyverse)


## processExpMat ####
#
# 
#


fobj <- Rflow(
  flowname = "CoExCorr",
  
  # #processExpMat_PODC(spc)
  # Rjob(
  #   jobName = "processExpMat_Os",
  #   source = "Rjobs/processExpMat.R",
  #   fun = "processExpMat_PODC",
  #   paramMat = data.frame( spc = "Os" )
  # ),
  # 
  # # processExpMat_PODC(spc)
  # Rjob(
  #   jobName = "processExpMat_At",
  #   source = "Rjobs/processExpMat.R",
  #   fun = "processExpMat_PODC",
  #   paramMat = data.frame( spc = "At" )
  # ),
  #
  # # processExpMat_PODC(spc)
  # Rjob(
  #   jobName = "processExpMat_Sl",
  #   source = "Rjobs/processExpMat.R",
  #   fun = "processExpMat_PODC",
  #   paramMat = data.frame( spc = "Sl" )
  # ),
  

  # processExpMat_ebi( spc)
  Rjob(
    jobName = "processExpMat_Zm",
    source = "Rjobs/processExpMat.R",
    fun = "processExpMat_ebi",
    paramMat = data.frame( spc = "Zm" )
  ),
  # processExpMat_ebi( inPath, outFile)
  Rjob(
    jobName = "processExpMat_Gm",
    source = "Rjobs/processExpMat.R",
    fun = "processExpMat_ebi",
    paramMat = data.frame( spc = "Gm" ) 
  )  
)


####
#
# convertTrees ####
#

arraySize = 12
fobj <- Rflow(flowname = "convertTrees",
              
              # extractTreeDataJob
              Rjob(sub_type = "scatter",
                   source = "Rjobs/extractTreeDataJob.R",
                   fun = "extractTreeDataJob",
                   paramMat = data.frame(arrayIdx = 1:arraySize,
                                         arraySize = arraySize,
                                         outDir = 'data/treeData') ),
              
              # mergeTreeDataJob
              Rjob(prev_jobs = "extractTreeDataJob",
                   dep_type = "gather",
                   sub_type = "serial",
                   source = "Rjobs/extractTreeDataJob.R",
                   fun = "mergeTreeDataJob",
                   paramMat = data.frame(outDir = 'data/treeData') ) )



## Full MI+CLR+CCS ####
#
# 
#
source("subflows/makeMIandCCSflow.R")
library(tibble)

dir.create("data/CCS",showWarnings = F)

MIFlowDef <- tribble(
  ~spc, ~cpuReq, ~memReq, ~expMatFile,
  "Gm", 46,      "4G",    "data/expMat/EBI_Gm.RDS",
  "Zm", 50,      "10G",   "data/expMat/EBI_Zm.RDS",
  "At", 30,      "3.5G",  "data/expMat/PODC_At.RDS",
  "Os", 22,      "2G",    "data/expMat/PODC_Os.RDS",
  "Sl", 22,      "2.5G",  "data/expMat/PODC_Sl.RDS"
)

fl <- makeMIandCCSflow(MIFlowDef)

fobj <- Rflow(flowname = "MI_CCS", fl)

fobj <- 
  fl %>% 
  startFromJob( startJob = c("GmprepMI", "ZmprepMI"))  %>% 
  Rflow(flowname = "RestCCS" ) %>% 
  plot_flow()

## calculate ranks from full CCS ####
#
#
#

source("subflows/makeAllOrthoRanksFlow.R")

fl <- makeAllOrthoRanksFlow()

fobj <-
  fl %>% 
  Rflow(flowname = "CalcOrthoRnks")

# 
# nGenes <- c(At=30919, Gm=54174, Os=45513, Sl=33624, Zm=41768)
# 
# spcs <- names(nGenes) %>% set_names(.,.)
# 
# spcPairs <-
#   crossing(spc1 = spcs,spc2 = spcs) %>%  # get all possible spc pairs
#   filter( spc1 < spc2 ) %>% # remove redundant pairs
#   # calculate memory requirements: !(turned out to be underestimated by up to 5G!)
#   mutate( memReq = paste0( ceiling(1.0*nGenes[spc1]*nGenes[spc2]*8/2^30 + 1),"G"))
# 
# fl <- 
#   pmap(spcPairs, function(spc1,spc2,memReq){
#     Rjob(jobName = paste0(spc1,spc2,"_calcRanks"),
#          source = "Rjobs/fullCCSRanksJob.R",
#          fun = "calcRanksJob",
#          memory_reserved = memReq,
#          paramMat = data.frame(spc1=spc1,spc2=spc2)
#          )
#   }) %>% 
#   reduce( flowbind )
# 



## withinSpecies subsets ####
#
#
#

source("subflows/makeWithinSpeciesFlow.R")

fl <- makeWithinSpeciesFlow()

fobj <- 
  fl %>% 
  # startFromJob( startJob = "SlwithinSpeciesJob") %>% 
  Rflow(flowname = "WithinSpecies" ) %>% 
  plot_flow()


## self CCS ####
#
# 
#

Rflow(
  flowname = "CoExCorr",
  Rjob(
    jobName = "At_selfCCS",
    source = "Rjobs/self_CCSjob.R", 
    fun = "self_CCSjob",
    memory_reserved = "20G",
    cpu_reserved = 10,
    paramMat = data.frame(
      miFile = "data/subsets/treesWithAtOs/At.mi",
      geneIDfile = "data/subsets/treesWithAtOs/At_geneIDs.RDS",
      refOrthosFile =  "data/subsets/treesWithAtOs/At_AtOs11_geneIDs.RDS",
      outFile = "data/subsets/treesWithAtOs/At_selfCCS.RDS" ,
      cores = 10
    )
  ),
  Rjob(
    jobName = "Os_selfCCS",
    source = "Rjobs/self_CCSjob.R", 
    fun = "self_CCSjob",
    memory_reserved = "20G",
    cpu_reserved = 10,
    paramMat = data.frame(
      miFile = "data/subsets/treesWithAtOs/Os.mi",
      geneIDfile = "data/subsets/treesWithAtOs/Os_geneIDs.RDS",
      refOrthosFile =  "data/subsets/treesWithAtOs/Os_AtOs11_geneIDs.RDS",
      outFile = "data/subsets/treesWithAtOs/Os_selfCCS.RDS" ,
      cores = 10
    )
  )
) -> fobj




##  submit flow ####
#
# 
#

opts_flow$set(verbose=2)
opts_flow$get()

fobj <- submit_flow(fobj,execute = F)
fobj <- submit_flow(fobj,execute = T, plot=F)

status(fobj)
# kill(fobj)
system("sacct")

plot_flow(fobj)

## resurrect fobj ####
#
#
fobj <- flowr:::read_fobj("/mnt/users/lagr/flowr/runs/RestCCS-foo-20170302-12-46-30-KdEy19MJ/")
fobj <- flowr:::read_fobj("/mnt/users/lagr/flowr/runs/ZmRestMI_CCS-foo-20170301-17-51-04-DYc2i1Zu/")

## example get jobIDs ####
#
# Example how to get jobIDs from jobNames
#

# get jobNames from flowList
jobs <- setdiff( getDownstreamFlowJobs(fl,c("ZmcalcMI")), "ZmcalcMI")

# get jobNames from fobj
jobs <- names(fobj@jobs) %>% .[grepl("Gm",.)]

# get jobIDs for jobs with specific names
jobIDs <- map(fobj@jobs[jobs], ~ .x@id)

# get all jobIDs
jobIDs <- map(fobj@jobs, ~ .x@id)


sacctFrmt<-"JobID,JobName,AllocTRES,MaxRSS,State,Elapsed,NodeList"

# sacct for specified jobIDs
jobStatus <- 
  jobIDs %>% 
  unlist() %>% 
  paste(collapse=",") %>% 
  paste0("sacct -P --format=",sacctFrmt," -j ",.) %>% 
  system(intern=T) %>% 
  paste(collapse = "\n") %>% 
  read_delim(delim="|")

jobStatus %>% 
  mutate( JobName = sub(".*_[0-9]{3}\\.(.*-[0-9]+$)","\\1",JobName))

## Example of rerun ####
#
# Example of rerun part of a job
#

# if fobj is available:
jobs <- names(fobj@jobs) %>% .[grepl("Gm",.)] # get all jobs with Gm in name
fobj <- flowr::rerun(fobj, select = jobs,kill = F)

####
#
#  example: using the startFromJob function
#

# library(magrittr)
# 
# startFromJob(flowList = fl,startJob = "Os_prepMI")
# fobj <- 
#   fl %>%
#   startFromJob("Os_prepMI") %>%
#   Rflow( flowname = "CoExCorr" )

####
#
#  example: resubmit from job 3
#

# fobj2 <- submit_flow(fobj,execute = T, .start_jid = 3) # execute

# it works even with the dependencies because the fobj contains the jobIDs for
# the previous job..

# However.. If I didn't have the fobj, would it work?.. the job ID's could have
# been acquired from the directory, but i'm not sure if that has been
# implemented

# Using the rerun command, the flow is loaded from the run directory
# WARNING! it kills all jobs by default unless kill=FALSE
# It does not seem to care about dependencies... So cannot be used to 

# flowr::rerun("/mnt/users/lagr/flowr/runs/PODC_At_Subsets-foo-20170113-11-11-43-FsVnOK8I",
#              start_from="PODC_At_HalfSubsetsJob3", kill=F) -> fobj

