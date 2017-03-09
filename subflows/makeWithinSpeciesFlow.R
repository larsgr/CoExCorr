source("R/flowrUtils.R")
library(purrr)

# Helper function that add specified dependency to first job in flowlist
addDepsToFlowList <- function(fl, prev_jobs, dep_type){
  fl$flowdef$prev_jobs[1] <- prev_jobs
  fl$flowdef$dep_type[1] <- dep_type
  return(fl)
}

# make jobs for running all reps for one spc
makeWithinSpeciesSpcFlow <- function(spc, 
                                     reps, 
                                     outDir = "data/subsets/withinSpecies",
                                     cores = 2,
                                     mem = "4G"){
  Rjob(
    sub_type = "scatter",
    startSleep = 1, # 1 second wait between starting each job
    cpu_reserved = cores,
    memory_reserved = mem,
    source = "Rjobs/withinSpeciesJobs.R", 
    fun = "withinSpeciesJob",
    paramMat = data.frame(
      spc = spc,
      repNr = 1:reps,
      outDir = outDir,
      cores = cores)
  )
}

# make prep job and jobs for running all reps for all spepcies
makeWithinSpeciesFlow <- function(numberOfGenes = 4000,
                                  reps = 10,
                                  outDir = "data/subsets/withinSpecies",
                                  cores = 2,
                                  mem = "4G",
                                  spcs = c("At","Gm","Os","Sl","Zm")){
  
  flowbind(
    Rjob(
      source = "Rjobs/withinSpeciesJobs.R", 
      fun = "withinSpeciesPrepJob",
      paramMat = data.frame( 
        numberOfGenes = numberOfGenes,
        reps = reps,
        outDir = outDir)
    ),
    spcs %>%                               # for each species
      map( makeWithinSpeciesSpcFlow,       # make species specific job
           reps = reps, outDir = outDir,
           cores = cores,mem = mem) %>% 
      map2( spcs, subFlow) %>%             # add spc to jobname
      map( addDepsToFlowList,              # add dependency to prep job
           prev_jobs = "withinSpeciesPrepJob", 
           dep_type = "burst") %>% 
      reduce( flowbind )                   # combine
  ) 
}