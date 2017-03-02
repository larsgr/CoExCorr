# Need to load the slurm environment to be able to submit jobs.
# But first need to init the environment modules manager
library(RLinuxModules) # devtools::install_github("larsgr/RLinuxModules")
moduleInit()
module("load slurm")

library(flowr)

# Set the module_cmds option
opts_flow$set(module_cmds = "source /etc/profile.d/modules.sh\nmodule load R/3.3.2")

source("R/flowrUtils.R")


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
  # # processExpMat_ebi( inPath, outFile)
  # Rjob(
  #   jobName = "processExpMat_Zm",
  #   source = "Rjobs/processExpMat.R",
  #   fun = "processExpMat_ebi",
  #   paramMat = data.frame( 
  #     inPath = "indata/RNAseqZm", 
  #     outFile = "data/expMat/EBI_Zm.RDS"
  #   ) 
  # ),
  
  # processExpMat_PODC(spc)
  Rjob(
    jobName = "processExpMat_Sl",
    source = "Rjobs/processExpMat.R",
    fun = "processExpMat_PODC",
    paramMat = data.frame( spc = "Sl" )
  ),
  
  # processExpMat_ebi( inPath, outFile)
  Rjob(
    jobName = "processExpMat_Gm",
    source = "Rjobs/processExpMat.R",
    fun = "processExpMat_ebi",
    paramMat = data.frame( 
      inPath = "indata/RNAseqGm", 
      outFile = "data/expMat/EBI_Gm.RDS"
    ) 
  )  
)

## PODC_At_Subsets ####
#
# 
#

fobj <- Rflow(flowname = "PODC_At_Subsets",
              
              # PODC_At_HalfSubsetsJob1
              Rjob(source = "Rjobs/PODC_At_Subsets.R",
                   fun = "PODC_At_HalfSubsetsJob1" ),
              
              # PODC_At_HalfSubsetsJob2
              Rjob(prev_jobs = "PODC_At_HalfSubsetsJob1", dep_type = "burst", sub_type = "scatter",
                   source = "Rjobs/PODC_At_Subsets.R",
                   fun = "PODC_At_HalfSubsetsJob2",
                   paramMat = expand.grid(subset_idx = 1:10, invertSubset = c(T,F)) ),
              
              # PODC_At_HalfSubsetsJob3
              Rjob(prev_jobs = "PODC_At_HalfSubsetsJob2", dep_type = "gather", sub_type = "scatter",
                   source = "Rjobs/PODC_At_Subsets.R",
                   fun = "PODC_At_HalfSubsetsJob3",
                   paramMat = data.frame(subset_idx = 1:10)),
              
              # PODC_At_HalfSubsetsJob4
              Rjob(prev_jobs = "PODC_At_HalfSubsetsJob1", dep_type = "burst", sub_type = "scatter",
                   source = "Rjobs/PODC_At_Subsets.R",
                   fun = "PODC_At_HalfSubsetsJob4",
                   paramMat = data.frame(subset_idx = 1:10))
)


flowList = makePODC_AtOs_flow()
flowList$flowdef$jobname
startJob = c("PODC_AtOs_mergeMI_At","PODC_AtOs_calcMI_Os")
makePODC_AtOs_flow() %>% startFromJob(startJob) -> x
x <- startFromJob(flowList,startJob)
x <- onlyTheseJobs(flowList,startJob)
fobj <- Rflow(flowname = "test",x)
plot_flow(Rflow(flowname = "test",x))



source("Rjobs/PODC_OsAt_Subsets.R")
fobj <- Rflow( flowname = "PODC_At_Subsets", 
               makePODC_At_Subsets_flow() %>% startFromJob("PODC_At_HalfSubsetsJob2") )
fobj <- Rflow( flowname = "PODC_Os_Subsets", makePODC_Os_Subsets_flow() )

# # example of resubmitting starting from specified jobs
# x <- makePODC_Os_Subsets_flow() %>% 
#   startFromJob(c("PODC_Os_HalfSubsetsJob3","PODC_Os_HalfSubsetsJob4"))
# fobj <- Rflow( flowname = "rerunPODC_Os_Subsets", x )

## PODC_AtOs ####
#
# 
#

source("Rjobs/PODC_AtOs.R")
fobj <- Rflow( flowname = "PODC_AtOs", makePODC_AtOs_flow())

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

## treesWithAtOs ####
#
# 
#

source("Rjobs/calcMI.R")


fobj <- 
  Rflow(
    flowname = "CoExCorr",
    subFlow(
      prefix="At_",
      makeMI_flow(
        expMatFile = "data/expMat/PODC_At.RDS",
        outDir = "data/subsets/treesWithAtOs",
        geneIDsubsetFile = "data/subsets/treesWithAtOs/At_geneIDs.RDS",
        prefix = "At",
        arraySize = 32,
        mem = "3G")
    ), 
    subFlow(
      prefix="Os_",
      makeMI_flow(
        expMatFile = "data/expMat/PODC_Os.RDS",
        outDir = "data/subsets/treesWithAtOs",
        geneIDsubsetFile = "data/subsets/treesWithAtOs/Os_geneIDs.RDS",
        prefix = "Os",
        arraySize = 16,
        mem = "1G")
    )
  )

#
# At Os CCS
#

Rflow(
  flowname = "CoExCorr",
  Rjob(
    jobName = "getAtOs11refs",
    source = "Rjobs/CLR_CCS.R", 
    fun = "get11orthos",
    paramMat = data.frame( 
      orthoFile = "indata/At_Os_orthologs.txt",
      geneIDfile1 = "data/subsets/treesWithAtOs/At_geneIDs.RDS",
      geneIDfile2 = "data/subsets/treesWithAtOs/Os_geneIDs.RDS",
      refOrthosFile1 =  "data/subsets/treesWithAtOs/At_AtOs11_geneIDs.RDS",
      refOrthosFile2 = "data/subsets/treesWithAtOs/Os_AtOs11_geneIDs.RDS"
    )
  ),
  Rjob( 
    dep_type = "serial", 
    prev_jobs = "getAtOs11refs", 
    jobName = "calcAtOsCCS",
    source = "Rjobs/CLR_CCS.R",
    fun = "MI_CLR_CCSjob",
    memory_reserved = "50G",
    cpu_reserved = 20,
    paramMat = data.frame(
      mi_file1 = "data/subsets/treesWithAtOs/At.mi",
      mi_file2 = "data/subsets/treesWithAtOs/Os.mi",
      geneIDfile1 = "data/subsets/treesWithAtOs/At_geneIDs.RDS",
      geneIDfile2 = "data/subsets/treesWithAtOs/Os_geneIDs.RDS",
      refOrthosFile1 =  "data/subsets/treesWithAtOs/At_AtOs11_geneIDs.RDS",
      refOrthosFile2 = "data/subsets/treesWithAtOs/Os_AtOs11_geneIDs.RDS",
      outFile = "data/subsets/treesWithAtOs/AtOs_CCS.RDS" ,
      cores = 20)
  )
) -> fobj


## Full MI matrices ####
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
  startFromJob( startJob = c("AtGmcalcCCS", "AtZmcalcCCS", "GmOscalcCCS", "GmSlcalcCCS", 
                             "GmZmcalcCCS", "OsZmcalcCCS", "SlZmcalcCCS"))  %>% 
  Rflow(flowname = "RestCCS" ) %>% 
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
fobj <- flowr:::read_fobj("/mnt/users/lagr/flowr/runs/restMI_CCS-foo-20170301-14-10-33-UOaq6i3n")
fobj <- flowr:::read_fobj("/mnt/users/lagr/flowr/runs/ZmRestMI_CCS-foo-20170301-17-51-04-DYc2i1Zu/")

## example get jobIDs ####
#
# Example how to get jobIDs for donwstream jobs in flowList
#
jobs <- setdiff( getDownstreamFlowJobs(fl,c("ZmcalcMI")), "ZmcalcMI")
jobs <- names(fobj@jobs) %>% .[grepl("calcCCS",.)]
map(fobj@jobs[jobs], ~ .x@id) %>% unlist() %>% paste(collapse=",")


## Example of rerun ####
#
# Example of rerun part of a job
#

# if fobj is available:
rerun(fobj, start_from = "processExpMat_Sl", select = "processExpMat_Sl",kill = F)
# if fobj not available
retVal <- rerun("/mnt/users/lagr/flowr/runs/CoExCorr-foo-20170222-19-52-41-OJ5qPXNG", start_from = "calcAtOsCCS",kill = F)

fobj <- retVal[[1]]
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

