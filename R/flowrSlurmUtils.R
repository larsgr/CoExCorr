
#
# SLURM specific helper functions that extend flowr functionality 
#

#' Update an already submitted flow
#' 
#' Compares flow definitions and updates the jobs with changes if they have not
#' already completed or is still running/pending.
#'
#' @param fobj flowr object
#' @param fl flow list with new definitions
#' @param dump write changes to disk?
#'
#' @return
#' @export
#'
#' @examples
updateFlow <- function(fobj, fl, dump=F){
  # generate new fobj from flow list
  newfobj <- Rflow("new",fl)
  
  # check that it is the same jobs
  stopifnot(all(names(newfobj@jobs) %in% names(fobj@jobs)))

  
  # check which jobs that are altered
  sapply( names(newfobj@jobs), function(jobName){
    compareJobs(fobj@jobs[[jobName]],newfobj@jobs[[jobName]])
  }) %>% 
    which %>% names -> alteredJobs
  
  if(length(alteredJobs)==0){
    cat("No change!\n")
    return(fobj)
  }
  
  
  # Get the status of the jobs
  # Do not allow altering jobs that are completed, running or pending
  cat("Checking job status...")

  conflictingJobStatus <-
    qstatus(fobj, use_cache = F) %>% # get status
    # complete or in queue?
    filter(completed - withError + pending + running > 0) %>% 
    # altered?
    filter(jobnm %in% alteredJobs) 

  if( nrow(conflictingJobStatus) > 0){
    cat("conflicting jobs:\n")
    kable(conflictingJobStatus)
    return(fobj)
  } else {
    cat("OK\n")
  }
  
  
  # replace altered jobs
  fobj@jobs[alteredJobs] <- newfobj@jobs[alteredJobs]
  
  # write changes to disk
  if(dump){
    flowr:::write_flow_details(fobj@flow_path, fobj = fobj,
                       flow_det = to_flowdet(fobj), 
                       flow_def = to_flowdef(fobj, verbose = 0), 
                       flow_mat = to_flowmat(fobj, verbose = 0))
  }
  return(fobj)
}

#' Get jobs that needs to be resubmitted
#'
#' Detects jobs that have not been submitted or have been completely killed or 
#' failed. Partially completed jobs will result in error.
#'
#' @param fobj flowr object
#'
#' @return character vector of job names
#' @export
getRestJobs <- function(fobj){
  qstat <- qstatus(fobj,use_cache = F)
  
  completedJobs <- 
    qstat %>% filter(completed - withError == total) %>% .$jobnm

  partialJobs <- 
    qstat %>% 
    # not completed
    filter( !(jobnm %in% completedJobs) ) %>% 
    # but has some pending/running/completed sub-jobs
    filter( pending + running + completed - withError > 0 ) %>% 
    .$jobnm
  
  if( length(partialJobs) > 0 ){
    kable(qstat[match(partialJobs, qstat$jobnm), ])
    stop("Flow includes partially completed sub-jobs.")
  }
  
  return( setdiff(qstat$jobnm, completedJobs))
}

#' Resubmit jobs that failed
#' 
#' Resubmit sub-jobs that failed or was cancelled without regenerating script.
#'
#' @param fobj flowr object
#' @param jobNames jobs to resubmit
#' @param dump write flow changes to disk?
#'
#' @return updated flowr object
#' @export
#'
#' @examples
resubmitFailedSubJobs <- function(fobj, jobNames, dump=T){
  for( jobName in jobNames){
    job <- fobj@jobs[[jobName]]  
    jobIDs <- job@id[isValidJobID(job@id)]
    
    # only do if job has been submitted before
    stopifnot(job@status == "submitted")

    # Figure out which task in the job that needs to be resubmitted
    if( length(jobIDs) > 0){
      queueStatus <- 
        getJobStatus(jobIDs) %>% 
        filter(JobID %in% jobIDs) # remove the batch entries
      
      # do not resubmit jobs that are pending, running or completed
      idx <- which(job@id %in% queueStatus$JobID[!(queueStatus$State %in% c("PENDING","RUNNING","COMPLETED"))])
    } else {
      idx <- 1:length(job@cmds)
    }
    cat(jobName,":",length(idx),"of",length(job@cmds),"tasks to resubmit.\n")
    
    # resubmit
    for( i in idx ){
      cmd <- flowr:::render_queue_cmd(jobj = job, fobj = fobj, index=i)
      cat("Resubmitting with command:",cmd,"\n")
      retVal <- system(cmd, intern = TRUE)
      cat(retVal,"\n")
      newJobID <- flowr:::parse_jobids(retVal, platform = job@platform)
      if( isValidJobID(newJobID) ){
        job@id[i] <- newJobID
      } else {
        cat("Invalid jobID! Stopping resubmission of this job\n")
        break
      }
    }
    
    # put the updated job back in the flow object 
    fobj@jobs[[jobName]] <- job
  }
  
  # write changes to disk
  if(dump){
    flowr:::write_flow_details(fobj@flow_path, fobj = fobj,
                               flow_det = to_flowdet(fobj), 
                               flow_def = to_flowdef(fobj, verbose = 0), 
                               flow_mat = to_flowmat(fobj, verbose = 0))
  }  
  return(fobj)
}


#' Get flow status including from slurm
#'
#' @param fobj flowr object
#' @param use_cache 
#' @param summarize If T then returns summary per flow job if F then returns
#' info per individual sub-jub
#' 
#'
#' @return
#' @export
#'
#' @examples
qstatus <- function(fobj, use_cache = T, summarize=T){
  
  # get status from flowr (trigger files)
  status <- get_status(to_flowdet(fobj),verbose = 0, use_cache = use_cache)
  # decode the exit code:
  status$completed <- !(is.na(status$exit_code) | status$exit_code == -1)
  status$withError <- !is.na(status$exit_code) & status$exit_code > 0
  
  # Is there some jobs that should be in the queue
  isSubmittedButNotFinished <- isValidJobID(status$job_sub_id) & 
    !status$completed
  
  status$qstate <- NA

  # get the job status from slurm if needed
  if( any(isSubmittedButNotFinished) ){
    queueStatus <- 
      getJobStatus( status$job_sub_id[isSubmittedButNotFinished]) %>% 
      filter(!grepl("batch",JobID))
    
    status$qstate[match(as.character(queueStatus$JobID),status$job_sub_id)] <- queueStatus$State
  }

  # summarize
  if( summarize ){
    summary <-
      status %>% 
      group_by(jobnm) %>% 
      summarize( total = length(job_sub_id),
                 submitted = sum(isValidJobID(job_sub_id)),
                 started = sum(started),
                 completed = sum(completed),
                 withError = sum(withError),
                 cancelled = sum( grepl("CANCELLED",qstate)),
                 pending = sum( grepl("PENDING",qstate)),
                 running = sum( grepl("RUNNING",qstate)),
                 failed = sum( grepl("FAILED",qstate)))
    return(summary[match(unique(status$jobnm),summary$jobnm),])
  } 
  return(status)
}

isValidJobID <- function(ids){
  !is.na(suppressWarnings(as.integer(ids)))
}

#' Get job status
#' 
#' Uses the sacct command to get status of a set of job ids.
#'
#' @param jobIDs slurm job ids
#'
#' @return data.frame
#' @export
#'
#' @examples
getJobStatus <- function(jobIDs){
  sacctFrmt <- "JobID,JobName,AllocTRES,MaxRSS,State,Elapsed,TotalCPU,CPUTime,NodeList"
  
  # sacct for specified jobIDs
  jobIDs %>% 
    unlist() %>% 
    paste(collapse=",") %>% 
    paste0("sacct -P --format=",sacctFrmt," -j ",.) %>% 
    system(intern=T) %>% 
    paste(collapse = "\n") %>% 
    read_delim(delim="|") %>% 
    mutate( JobName = sub(".*_[0-9]{3}\\.(.*-[0-9]+$)","\\1",JobName))

}
# # TODO: a solution to merge the batch entries (which contains the MaxRSS etc.)
# # with the job entries (which contains the JobName).
#
# jobstats <- getJobStatus(jobIDs)
# 
# jobstats %>% 
#   filter( !grepl("batch",JobID)) %>% 
#   mutate( MaxRSS = jobstats$MaxRSS[match(paste0(JobID,".batch"),jobstats$JobID)])


compareJobs <- function(j1,j2){
  # sanity check
  stopifnot(all(slotNames(j1) == slotNames(j2)))
  
  changed <- F
  
  # compare job definitions
  for(slotName in c("platform","walltime","memory","cpu","nodes","queue",
                    "submission_type","dependency_type","previous_job",
                    "extra_opts")){
    
    v1 <- slot(j1,slotName)
    v2 <- slot(j2,slotName)
    if(!identical(v1, v2)){
      cat(j1@name,":",slotName,"changed from",v1,"to",v2,"\n")
      changed <- T
    }
  }
  
  # check if commands have been changed
  if( !identical(j1@cmds,j2@cmds)){
    cat(j1@name,": Commands have been changed\n")
    changed <- T
  }
  
  return(changed)
}



#' kill all pending jobs
#' 
#' Checks the queue to see if any of the jobs in the flow are pending then
#' cancels them.
#'
#' @param fobj flowr object
#'
#' @return nothing
#' @export
#'
killPendingJobs <- function(fobj){
  # for each job
  for( job in fobj@jobs){
    # get job IDs
    jobIDs <- job@id[isValidJobID(job@id)]
    
    # skip if no IDs available
    if( length(jobIDs)==0 )
      next
    
    # get job state
    queueStatus <- 
      getJobStatus(jobIDs) %>% 
      filter(JobID %in% jobIDs) # remove the batch entries
    
    # get id of pending jobs
    pendingIDs <- queueStatus$JobID[queueStatus$State == "PENDING"]
    
    if( length(pendingIDs)==0 )
      next
    
    cat( job@name,"has",length(pendingIDs),"pending job(s).\n")
    CMD <- paste("scancel",paste(pendingIDs,collapse=" "))
    cat("Executing CMD:",CMD,"\n")
    system(CMD)
  }
}
