library(tidyverse)

makeNSamplesFlow <- function(outdir = "data/subsets/nSamples", 
                             nStudyReps=10,
                             nStudy=c(25,50,100,200,400,800),
                             nSampleReps=5,
                             nSample=c(12,25,50,100,200,400),
                             nTotalSamples=1363){
  
  flowbind(
    # nSamplesInit(outdir, nStudyReps, nSampleReps)
    Rjob(
      source = "Rjobs/nSamplesJob.R",
      fun = "nSamplesInit",
      paramMat = data.frame( outdir=outdir, 
                             nStudyReps=nStudyReps, 
                             nSampleReps=nSampleReps )
    ),
    
    # nSamplesJob(n, subsetIdx, subsetFile, outfile, mrOs.bm.desc, xAt.bm.desc)
    crossing( n=nStudy, subsetIdx=1:nStudyReps) %>% 
      mutate( subsetFile = file.path(outdir,"rndStudies.RDS"), 
              outfile = paste0(outdir,"/rnks_rndStudies_",subsetIdx,"_",n,".RDS"),
              mrOs.bm.desc = file.path(outdir,"mrOs.bm.desc"),
              xAt.bm.desc = file.path(outdir,"xAt.bm.desc") ) %>% 
      as.data.frame() %>% 
      Rjob(paramMat = .,
           jobName = "nSamplesRndStudies",
           source = "Rjobs/nSamplesJob.R",
           fun = "nSamplesJob",
           prev_jobs = "nSamplesInit"),
    
    
    crossing( n=nSample, subsetIdx=1:nSampleReps) %>% 
      mutate( subsetFile = file.path(outdir,"rndSamples.RDS"), 
              outfile = paste0(outdir,"/rnks_rndSamples_",subsetIdx,"_",n,".RDS"),
              mrOs.bm.desc = file.path(outdir,"mrOs.bm.desc"),
              xAt.bm.desc = file.path(outdir,"xAt.bm.desc") ) %>% 
      as.data.frame() %>% 
      Rjob(paramMat = .,
           jobName = "nSamplesRndSamples",
           source = "Rjobs/nSamplesJob.R",
           fun = "nSamplesJob",
           prev_jobs = "nSamplesInit"),
    
    Rjob(
      paramMat = data.frame(
        n = nTotalSamples,
        subsetIdx = 1, # order doesn't matter when all samples are included
        subsetFile = file.path(outdir,"rndStudies.RDS"), 
        outfile = paste0(outdir,"/rnks_allSamples.RDS"),
        mrOs.bm.desc = file.path(outdir,"mrOs.bm.desc"),
        xAt.bm.desc = file.path(outdir,"xAt.bm.desc")
      ),
      jobName = "nSamplesAllSamples",
      source = "Rjobs/nSamplesJob.R",
      fun = "nSamplesJob",
      prev_jobs = "nSamplesInit" )
    
  )
}


