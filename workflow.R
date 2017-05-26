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
source("R/flowrSlurmUtils.R")


## processExpMat ####
#
# 
#
source("subflows/makeProcessExpMatFlow.R")

fl <- makeProcessExpMatFlow()

fobj <- Rflow(flowname = "CoExCorr",fl)


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


## Full PCC+MR+CCS ####
#
# 
#
source("subflows/makePCCandMRflow.R")

PCCMRFlowDef <-
  tribble(
    ~spc, ~memReq, ~cores,
    "At",   "30G",     20,
    "Os",   "60G",     20,
    "Sl",   "34G",     10,
    "Gm",   "86G",     30,
    "Zm",   "50G",     20 )

CCSRanksFlowDef <- tribble(
  ~spc1, ~spc2, ~cores, ~memReq,
  "At", "Gm", 20, "71G",
  "At", "Os", 20, "57G",
  "At", "Sl", 20, "40G",
  "At", "Zm", 20, "50G",
  "Gm", "Os", 20, "94G",
  "Gm", "Sl", 20, "75G",
  "Gm", "Zm", 20, "86G",
  "Os", "Sl", 20, "60G",
  "Os", "Zm", 20, "71G",
  "Sl", "Zm", 20, "53G")


fl <- makePCC_MR_CCS_RanksFlow(PCCMRFlowDef,CCSRanksFlowDef,
                               outDirMR="data/MR", outDirRnks="data/ranksMR")


fobj <- Rflow(flowname = "FullPCCMR",fl)

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
  "At", 30,      "4G",    "data/expMat/PODC_At.RDS",
  "Os", 22,      "2G",    "data/expMat/PODC_Os.RDS",
  "Sl", 22,      "3G",    "data/expMat/PODC_Sl.RDS"
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


## Compare methods ####
#
#
#

source("subflows/makeCompareMethodsFlow.R")

fl <- makeCompareMethodsAllFlow()

fobj <- Rflow(flowname = "cm",fl)


##  nSamples ####
#
# 
#

source("subflows/makeNSamplesFlow.R")

fl <- makeNSamplesFlow()

fobj <- Rflow(flowname = "nSamples", fl)

##  SUBMIT FLOW ####
#
# 
#

opts_flow$set(verbose=2)
opts_flow$get()

fobj <- submit_flow(fobj,execute = F)
fobj <- submit_flow(fobj,execute = T, plot=F)

# check status
status(fobj)

# check status (extended version)
knitr::kable(qstatus(fobj,use_cache = F))

# kill(fobj)
system("sacct")

plot_flow(fobj, detailed=FALSE)

## resurrect fobj ####
#
#

fobj <- flowr:::read_fobj("/mnt/users/lagr/flowr/runs/MI_CCS-foo-20170419-17-55-29-PjHcJ5XN/")
fobj <- flowr:::read_fobj("/mnt/users/lagr/flowr/runs/cm-foo-20170504-10-53-46-LMDX0Gph/")

## example get jobIDs ####
#
# Example how to get jobIDs from jobNames
#

# get jobNames from flowList
jobs <- setdiff( getDownstreamFlowJobs(fl,c("ZmcalcMI")), "ZmcalcMI")

# get jobNames from fobj
jobs <- names(fobj@jobs) %>% .[grepl("calcCCS",.)]

# get jobIDs for jobs with specific names
jobIDs <- unlist(map(fobj@jobs[jobs], ~ .x@id))

jobIDs <- unlist(fobj@jobs$AtcalcMI@id)

# get all jobIDs
jobIDs <- map(fobj@jobs, ~ .x@id)

jobIDs <- jobIDs[!sapply(jobIDs,"==","character(0)")]

getJobStatus(jobIDs)

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
# the previous job.. But it does not let you select specific jobs




####
#
# example: using rerun
#
# WARNING! it kills all jobs by default unless kill=FALSE
# rerun() removes dependancies from jobs that are not selected for rerunning, so 
# cannot be used for jobs that depend on running jobs.

jobs <- getRestJobs(fobj) # get name of jobs that needs to be rerun
fobj <- flowr::rerun(fobj, select=jobs, kill=F) # resubmit




####
#
# example: resume after failed subjobs
#

# first stop all pending jobs (assume they are waiting for the failed jobs)
killPendingJobs(fobj)

# resubmitted only failed subjobs of a specified job:
fobj <- resubmitFailedSubJobs(fobj,jobNames = "calcCor")
# Note.. resubmitFailedSubJobs does not reset the trigger files so the status

# note: is there any way to submit the downstream jobs?




####
#
# example: fix for (launch failed requeued held)
#

qState <- system("squeue -u lagr -h",intern=T)
ids <- sub("[ ]+([0-9]+).*","\\1",qState)[grepl("(launch failed requeued held)",qState)]

system(paste("scontrol hold",paste(ids,collapse=",")))
system(paste("scontrol release",paste(ids,collapse=",")))

for(id in ids[11:12]){
  # cmd <- paste("scontrol hold",id)
  cmd <- paste("scontrol release",id)
  cat(cmd,"\n")
  system(cmd)
  # Sys.sleep(1)
}