#' Create flowr job for R script
#' 
#' Rjob generates multiple commandline calls to function "fun" which is sourced 
#' from "source" using the parameters in "paramMat". One call for each row.
#'
#' @param source .R file to be source that contains the function to be called
#' @param fun name of function to be called
#' @param paramMat data.frame with arguments to the function. Each
#' @param jobName name of job
#' @param prev_jobs name of job
#' @param dep_type "none", "serial", "burst" or "merge"
#' @param sub_type "serial": function call will be run after each other, 
#'                 or "scatter": each function call will be run in parallel
#' @param cpu_reserved number of CPU's per job
#' @param memory_reserved memory to reserve for this job. e.g. "4G" or "200M"
#' @param walltime 
#' @param queue 
#' @param platform 
#' @param startSleep number of seconds to delay at the start of each array job
#'
#' @return list with two data.frames (flowmat and flowdef)
#' @export
#'
Rjob <- function( source, fun, paramMat=data.frame(), jobName=fun,
                  prev_jobs="none", dep_type="none", sub_type="serial",
                  cpu_reserved = 1, 
                  memory_reserved = "4G",
                  walltime = "10-00:00:00",
                  queue = "verysmallmem,cigene",
                  platform = "slurm",
                  startSleep = 0){
  
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
    sapply(mapply(call,name=fun,..., USE.NAMES=F),function(funCall){
      # convert call to string
      funCallStr <- capture.output(funCall)
      # long calls will use multiple lines. Force it to be one line by removing 
      # leading spaces and collapsing the string
      paste(sub("^\\s+", "", funCallStr),collapse="")
    })
  })

  
  # Generate source() function call as string
  sourceCallStr <- capture.output(call("source",source))
  
  # Glue them together and add quotes 
  exprStr <- shQuote(paste(sourceCallStr,funCallStr,sep = "; "))
  
  # generate sleep command in the beginning
  if(startSleep == 0){
    sleepCmd <- ""
  } else {
    if(length(startSleep) != 1 | !is.numeric(startSleep) | startSleep < 0){
      stop("startSleep must be positive numeric of length 1")
    }
    sleepCmd <- paste("sleep", (0:(length(exprStr)-1))*startSleep ,"&& ")
  }
  
  flowmat <- data.frame( stringsAsFactors = F,
                         samplename = "foo",
                         jobname = jobName,
                         cmd = paste0(sleepCmd,
                                      "cd ",getwd()," && ",  # execute in project directory
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

# Include this/these job(s) and all downstream jobs, 
# remove up-stream dependencies
startFromJob <- function(flowList, startJob){
  subsetFlowJobs(flowList, getDownstreamFlowJobs(flowList, startJob))
}

getDownstreamFlowJobs <- function(flowList, startJob){
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
  
  # return name of jobs
  return( names(jobsIncluded)[jobsIncluded] )
}



#' Fix dependencies in flowList
#' 
#' Remove dependencies to jobs that are not in the list and alter the 
#' dependency types of affected jobs if needed
#'
#' @param flowList list containing flowmat and flowdef tables
#'
#' @return flowList
fixFlowDependencies <- function(flowList){
  
  flowdef <- flowList$flowdef
  flowmat <- flowList$flowmat
  
  # split the old dependencies and replace "none" with ""
  oldDeps <- setNames( strsplit( sub("^none$","",flowdef$prev_jobs),","),
                       flowdef$jobname)
  # remove non-existant dependencies
  newDeps <- lapply(oldDeps, function(x) {x[x %in% flowdef$jobname]})
  
  # fix dependency types
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



#' Subset jobs in flowList
#' 
#' Removes all jobs from flowList that is not in the given set of jobs
#'
#' @param flowList list containing flowmat and flowdef tables
#' @param jobs character vector with names of jobs to include
#' @param fixDeps logical, if TRUE then fixFlowDependencies() is called after removing jobs
#'
#' @return flowList
subsetFlowJobs <- function(flowList, jobs, fixDeps=T){
  # check if jobs exists
  if( !all(jobs %in% flowList$flowdef$jobname) )
    stop("Given job(s) not found in flowList: ",
         paste(jobs[!(jobs %in% flowList$flowdef$jobname)],collapse=", "),"")
  
  
  # create reduced flowdef and flowmat
  retFlowList = list(
    flowdef = flowList$flowdef[flowList$flowdef$jobname %in% jobs, ],
    flowmat = flowList$flowmat[flowList$flowmat$jobname %in% jobs, ]
  )
  
  if(fixDeps){
    return( fixFlowDependencies(retFlowList) )
  } else {
    return( retFlowList )
  }
}

# combines lists containing flowmat and flowdef tables by rbind
flowbind <- function(...){
  flowLists <- list(...)
  list( flowmat = do.call(rbind,lapply(flowLists,"[[","flowmat")),
        flowdef = do.call(rbind,lapply(flowLists,"[[","flowdef")))
}

# Add prefix to all jobnames and dependencies
subFlow <- function(flowList, prefix){
  
  # first the flowmat (command list)
  flowList$flowmat$jobname <- paste0(prefix,flowList$flowmat$jobname)
  
  # Then the flowdef (dependencies and resource requirements)
  # the jobnames are easy:
  flowList$flowdef$jobname <- paste0(prefix,flowList$flowdef$jobname)
  
  
  # For the dependencies we need to take into account that it could "none" or
  # that there could be multiple dependencies separated by ","
  
  sapply( strsplit(flowList$flowdef$prev_jobs,split = ","), function(deps){
    if( (length(deps) == 1) & (deps[1] == "none") ){
      
      # no dependencies
      return(deps) # just leave it as it is
      
    } else if( length(deps) == 0){
      
      # also allow the possibility that deps==""
      return("") # deps = ""
      
    } else {
      
      # add prefix
      return(paste0(prefix, deps, collapse=","))
      
    }
  }) -> flowList$flowdef$prev_jobs
  
  # return the altered flowList
  return(flowList)
}

Rflow <- function(flowname, ...){
  flowList <- flowbind(...)
  
  to_flow(x = flowList$flowmat, 
          def = as.flowdef(flowList$flowdef),
          flowname = flowname)
}