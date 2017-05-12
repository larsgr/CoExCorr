library(tidyverse)

makeNSamplesFlow <- function(outdir = "data/subsets/nSamples", 
                             nStudyReps=10,
                             nStudy=c(25,50,100,200,400,800),
                             nSampleReps=5,
                             nSample=c(12,25,50,100,200,400),
                             nTotalSamples=1363,
                             threads = 4){
  
  flowbind(
    # nSamplesInit(outdir, nStudyReps, nSampleReps)
    Rjob(
      source = "Rjobs/nSamplesJob.R",
      fun = "nSamplesInit",
      paramMat = data.frame( outdir=outdir, 
                             nStudyReps=nStudyReps, 
                             nSampleReps=nSampleReps )
    ),
    
    # nSamplesJob(n, subsetIdx, subsetFile, outfile, clrOs.bm.desc, xAt.bm.desc, threads)
    crossing( n=nStudy, subsetIdx=1:nStudyReps) %>% 
      mutate( subsetFile = file.path(outdir,"rndStudies.RDS"), 
              outfile = paste0(outdir,"/rnks_rndStudies_",subsetIdx,"_",n,".RDS"),
              clrOs.bm.desc = file.path(outdir,"clrOs.bm.desc"),
              xAt.bm.desc = file.path(outdir,"xAt.bm.desc"),
              threads = threads ) %>% 
      as.data.frame() %>% 
      Rjob(paramMat = .,
           jobName = "nSamplesRndStudies",
           source = "Rjobs/nSamplesJob.R",
           fun = "nSamplesJob",
           prev_jobs = "nSamplesInit",
           cpu_reserved = threads),
    
    
    crossing( n=nSample, subsetIdx=1:nSampleReps) %>% 
      mutate( subsetFile = file.path(outdir,"rndSamples.RDS"), 
              outfile = paste0(outdir,"/rnks_rndSamples_",subsetIdx,"_",n,".RDS"),
              clrOs.bm.desc = file.path(outdir,"clrOs.bm.desc"),
              xAt.bm.desc = file.path(outdir,"xAt.bm.desc"),
              threads = threads ) %>% 
      as.data.frame() %>% 
      Rjob(paramMat = .,
           jobName = "nSamplesRndSamples",
           source = "Rjobs/nSamplesJob.R",
           fun = "nSamplesJob",
           prev_jobs = "nSamplesInit",
           cpu_reserved = threads ),
    
    Rjob(
      paramMat = data.frame(
        n = nTotalSamples,
        subsetIdx = 1, # order doesn't matter when all samples are included
        subsetFile = file.path(outdir,"rndStudies.RDS"), 
        outfile = paste0(outdir,"/rnks_allSamples.RDS"),
        clrOs.bm.desc = file.path(outdir,"clrOs.bm.desc"),
        xAt.bm.desc = file.path(outdir,"xAt.bm.desc"),
        threads = threads
      ),
      jobName = "nSamplesAllSamples",
      source = "Rjobs/nSamplesJob.R",
      fun = "nSamplesJob",
      prev_jobs = "nSamplesInit",
      cpu_reserved = threads )
    
  )
}


