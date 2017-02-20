# Need to load the slurm environment to be able to submit jobs.
# But first need to init the environment modules manager
library(RLinuxModules) # devtools::install_github("larsgr/RLinuxModules")
moduleInit()
module("load slurm")

library(flowr)

# Set the module_cmds option
opts_flow$set(module_cmds = "source /etc/profile.d/modules.sh\nmodule load R/3.3.1")



source("R/flowrUtils.R")


####
#
# processExpMat PODC
#

fobj <- Rflow(flowname = "processExpMat",
              # processExpMat_PODC_Os
              Rjob(source = "Rjobs/processExpMat.R",
                   fun = "processExpMat_PODC_Os" ),
              
              # processExpMat_PODC_At
              Rjob(source = "Rjobs/processExpMat.R",
                   fun = "processExpMat_PODC_At" )
              )


####
#
# PODC_At_Subsets
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

####
#
# PODC_AtOs
#

source("Rjobs/PODC_AtOs.R")
fobj <- Rflow( flowname = "PODC_AtOs", makePODC_AtOs_flow())

####
#
# convertTrees
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


###
#
# submit flow
#
fobj <- submit_flow(fobj,execute = F)
fobj <- submit_flow(fobj,execute = T)
status(fobj)
# kill(fobj)
system("sacct")

plot_flow(fobj)
####
#
#  example: resubmit from job 3
#

# fobj2 <- submit_flow(fobj,execute = T, .start_jid = 3) # execute

fobj <- Rflow(flowname = "PODC_At_Subsets",

              # PODC_At_HalfSubsetsJob3
              Rjob(#prev_jobs = "PODC_At_HalfSubsetsJob2", dep_type = "gather",
                sub_type = "scatter",
                   source = "Rjobs/PODC_At_Subsets.R",
                   fun = "PODC_At_HalfSubsetsJob3",
                   paramMat = data.frame(subset_idx = 1:10)),

              # PODC_At_HalfSubsetsJob4
              Rjob(#prev_jobs = "PODC_At_HalfSubsetsJob1", dep_type = "burst",
                sub_type = "scatter",
                   source = "Rjobs/PODC_At_Subsets.R",
                   fun = "PODC_At_HalfSubsetsJob4",
                   paramMat = data.frame(subset_idx = 1:10))
)

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

