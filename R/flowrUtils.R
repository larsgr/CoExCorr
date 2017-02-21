# Rjob generates multiple commandline calls to function "fun" which is sourced 
# from "source" using the parameters in "paramMat". One call for each row.
Rjob <- function( source, fun, paramMat=data.frame(), jobName=fun,
                  prev_jobs="none", dep_type="none", sub_type="serial",
                  cpu_reserved = 1, 
                  memory_reserved = "4G",
                  walltime = "10-00:00:00",
                  queue = "verysmallmem,cigene",
                  platform = "slurm" ){
  
  # convert any factors to strings in paramMat
  paramMat <- rapply(paramMat, as.character, classes="factor", how="replace")
  
  # Generate R function calls as strings
  #
  # Explanation of the obfuscated code below:
  #   "paramMat" is a data.frame where each row contains parameters to "fun".
  #   The do.call function effectively converts the paramMat to (...) i.e.
  #   a list of arguments, which is sent to mapply, which applies these 
  #   arguments to call(), which generates a call object. sapply() is then used
  #   to apply capture.output() to all the call objects, so they become strings.
  funCallStr <- do.call(args = paramMat, what = function(...){
    sapply(mapply(call,name=fun,..., USE.NAMES=F),capture.output)
  })
  
  # Generate source() function call as string
  sourceCallStr <- capture.output(call("source",source))
  
  # Glue them together and add quotes 
  exprStr <- shQuote(paste(sourceCallStr,funCallStr,sep = "; "))  
  
  flowmat <- data.frame( stringsAsFactors = F,
                         samplename = "foo",
                         jobname = jobName,
                         cmd = paste0("cd ",getwd()," && ",  # execute in project directory
                                      "Rscript -e ", exprStr))
  
  flowdef <- data.frame(stringsAsFactors = F,
                        jobname = jobName,
                        sub_type = sub_type, 
                        prev_jobs = prev_jobs,
                        dep_type = dep_type,
                        platform = platform,
                        cpu_reserved = cpu_reserved,
                        queue = queue,
                        memory_reserved = memory_reserved,
                        walltime = walltime)
  
  return(list(flowmat = flowmat,flowdef = flowdef))  
}

# Include this/these job(s) and all downstream jobs
startFromJob <- function(flowList, startJob){
  # check if startJob exists
  if( !all(startJob %in% flowList$flowdef$jobname) )
    stop("Given startJob(s) not found in flowList: ",
         paste(startJob[!(startJob %in% flowList$flowdef$jobname)],collapse=", "),"")
  
  if(nrow(flowList$flowdef) < 1)
    stop("Empty flowList!")
  
  jobDeps <- setNames( strsplit(flowList$flowdef$prev_jobs,","),
                       flowList$flowdef$jobname)
  
  jobsIncluded <- setNames( flowList$flowdef$jobname %in% startJob,
                            flowList$flowdef$jobname)
  
  # recursive function that walk up the dependency graph
  walkDeps <- function(job,level=0){
    # cat("\n",rep(". ",level),job,"  ",sep = "")
    if(job %in% c("","none")){
      # cat("reached end of dependency graph. No dependencies here\n")
      return(F) # reached end of dependency graph. No dependencies here
    } else {
      if( jobsIncluded[job] ){
        # cat("This job is already included. Found dependency!\n")
        return(T) # reached end of dependency graph. Found dependency
      } else {
        if(any(sapply(jobDeps[[job]],FUN=walkDeps,level=level+1))){
          # there are dependencies upstream
          jobsIncluded[job] <<- T # include this job
          # cat(rep(". ",level),"Include ",job,"\n",sep = "")
          return(T)
        }
      }
    }
    return(F)
  }
  
  
  # get all downstream jobs
  for( job in flowList$flowdef$jobname){
    walkDeps(job)
  }
  
  # remove jobs dependencies that are not to be included
  return( onlyTheseJobs(flowList,jobs = names(jobsIncluded)[jobsIncluded]) )
}

# remove all but these jobs. Also remove dependencies.
onlyTheseJobs <- function(flowList, jobs){
  # check if jobs exists
  if( !all(jobs %in% flowList$flowdef$jobname) )
    stop("Given startJob(s) not found in flowList: ",
         paste(jobs[!(jobs %in% flowList$flowdef$jobname)],collapse=", "),"")
  
  
  # create reduced flowdef and flowmat
  flowdef <- flowList$flowdef[flowList$flowdef$jobname %in% jobs, ]
  flowmat <- flowList$flowmat[flowList$flowmat$jobname %in% jobs, ]
  
  # fix dependencies
  flowdef$prev_jobs[flowdef$prev_jobs=="none"] <- character(0)
  oldDeps <- setNames( strsplit(flowdef$prev_jobs,","),
                       flowdef$jobname)
  # remove non-existant dependencies
  newDeps <- lapply(oldDeps, function(x) {x[x %in% jobs]})
  
  # fix dependency type
  for(i in seq_along(flowdef$jobname)){
    job = flowdef$jobname[i]
    if(!identical(oldDeps[[job]],newDeps[[job]])){
      # dependencies have changed!
      flowdef$prev_jobs[i] <- paste(newDeps[[job]],collapse=",")
      
      # three possibilities:
      
      if(length(newDeps[[job]])==0){
        
        # 1)  burst/gather/serial -> none
        
        cat(job,"changed from",flowdef$dep_type[i],"to none\n")
        flowdef$dep_type[i] <- "none"
        flowdef$prev_jobs[i] <- "none"
        
      } else if(length(newDeps[[job]])==1 & 
                flowdef$sub_type[match(newDeps[[job]],flowdef$jobname)] == "serial"){
        
        # 2)  gather -> serial (several jobs to one serial job)
        
        cat(job,"changed from",flowdef$dep_type[i],"to serial\n")
        flowdef$dep_type[i] <- "serial"
        
      } else {
        
        # 3)  gather -> gather (several jobs to less, but still several jobs)
        
        if(flowdef$dep_type[i] != "gather")
          warning("Job (",job,") has job_type: '",flowdef$dep_type[i],"'. Expected 'gather'! ")
        
      }
    }
  }
  
  return(list(flowdef=flowdef,flowmat=flowmat))
}


# combines lists containing flowmat and flowdef tables by rbind
flowbind <- function(...){
  flowLists <- list(...)
  list( flowmat = do.call(rbind,lapply(flowLists,"[[","flowmat")),
        flowdef = do.call(rbind,lapply(flowLists,"[[","flowdef")))
}

Rflow <- function(flowname, ...){
  flowList <- flowbind(...)
  
  to_flow(x = flowList$flowmat, def = as.flowdef(flowList$flowdef),
          module_cmd = "source /etc/profile.d/modules.sh\nmodule load R/3.3.1",
          flowname = flowname)
}