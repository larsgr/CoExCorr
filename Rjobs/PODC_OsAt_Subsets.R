#
# Functions used in jobs:
#
source("R/calcRanks.R")

calcPODC_HalfSubsets <- function(PODC_metadata_file,
                                 outDir){

  library(readxl)
  
  PODC_metadata <- read_excel(PODC_metadata_file,col_types = rep("text",15))
  
  studyID <- unique(PODC_metadata$Study)
  
  # set random seed
  set.seed(700) 

  # randomly select half of the studies. 10 times
  subsetsStudies <- replicate(10,sample(studyID,size = length(studyID) %/% 2))
  
  # get the sampleIDs for the selected studies
  apply( subsetsStudies, 2, function(subsetStudies){
    PODC_metadata$Run[PODC_metadata$Study %in% subsetStudies]
  }) -> subsetsSamples
  
  # save the sampleIDs to file (in case we need it later)
  dir.create(outDir,showWarnings = F, recursive = T)
  
  saveRDS(subsetsSamples,file.path(outDir,"subsetsSamples.RDS"))
}


#### job step 1
#
# Prepare samplesubsets and genesubset ExpMat
#
# calc subsets (in a deterministic way)
# load expMat
# save expMat subset with only 1:1 genes
PODC_Os_HalfSubsetsJob1 <- function(outDir = "data/subsets/PODC_Os_HalfSubsets",
                                    spc="Os",
                                    expMatFile="data/expMat/PODC_Os.RDS",
                                    PODC_metadata_file = "indata/oryza_sativa_sample.xlsx"){
  
  calcPODC_HalfSubsets( PODC_metadata_file, outDir)
  
  expMat <- readRDS(expMatFile)
  
  # subset genes, here I use the Ensembl At:Os 1:1 orthologs
  AtOsOrthos <- read.table("indata/At_Os_orthologs.txt",header=T,stringsAsFactors = F)
  geneIDs <- AtOsOrthos[[spc]][AtOsOrthos$otype=="1:1"]
  geneIDs <- geneIDs[geneIDs %in% rownames(expMat)]
  expMat11 <- expMat[geneIDs,]
  
  
  # save to file
  saveRDS(geneIDs,file.path(outDir,"geneIDs.RDS"))
  saveRDS(expMat11, file.path(outDir,"expMat11.RDS"))
}

#### job step 2 (parallel execution for each subset including the inverse subsets)
#
# Calculate MI
#
# run genepair to calculat MI and output as triangular matrix
PODC_Os_HalfSubsetsJob2 <- function(subset_idx, invertSubset, 
                                    outDir = "data/subsets/PODC_Os_HalfSubsets",
                                    expMatFile = file.path(outDir,"expMat11.RDS") ){
  
  subsetsSamples <- readRDS(file.path(outDir,"subsetsSamples.RDS"))
  expMat11 <- readRDS(expMatFile)
  
  samples <- subsetsSamples[[subset_idx]]
    
  # subset samples (if invertSubset select the opposite subset)
  x <- expMat11[ , xor(colnames(expMat11) %in% samples,invertSubset)]
  
  # write expression data in the format genepair understands....
  tmpExpMatFile <- file.path(outDir,paste0("expMat_",subset_idx, 
                                           ifelse(invertSubset,"_2","_1")))
  cat("writing",tmpExpMatFile,"...\n")
  write.table(x,file = tmpExpMatFile,sep = " ",row.names = F, col.names = F)
  
  miFile <- file.path(outDir,paste0("mi_",subset_idx, 
                                    ifelse(invertSubset,"_2","_1")))
  
  # genepair [6|1] microarray_file.txt mi_file.txt ngenes nexp fromgene togene num_bins spline_order
  cmd <- paste("~/genepair/genepair 1", tmpExpMatFile, miFile, nrow(x), ncol(x), 0, nrow(x), 7, 3)

  cat("CMD:",cmd,"\n")
  
  # run genepair to calculate MI
  system(cmd)
  
  # cleanup
  unlink(tmpExpMatFile) # remove expression data
  unlink(paste0(miFile,".mi_log.txt"))
}
  

#### job step 3 (parallel for each pair of subsets)
#
# CCS and ranks based in MI and MI+CLR
#
# calcCLR for both subsets in pair
# calcRanks
# save ranks
PODC_Os_HalfSubsetsJob3 <- function(subset_idx,
                                    outDir = "data/subsets/PODC_Os_HalfSubsets"){
  source("R/LoadComPlExData.R") # loadTriMatrix
  source("R/CLR.R")

  #
  # define functions
  #
  
  mySaveRDS <- function(x, prefix){
    filename <- file.path(outDir,
                          paste0(prefix,"_",subset_idx,".RDS"))
    cat("Saving",filename,"...\n")
    saveRDS(x,file = filename) 
  }
  

  #
  # Loading MI matrix
  #
  
  geneIDs <- readRDS(file.path(outDir,"geneIDs.RDS"))
  
  miFiles <- file.path(outDir,paste0("mi_",subset_idx, "_",1:2))

  lapply(miFiles,function(miFile){
    cat("Loading",miFile,"...\n")
    loadTriMatrix(geneIds = geneIDs, filename = miFile)
  }) -> mi

  
  #
  # Generate and save CCS and ranks using only MI
  #
  
  cat("Calculating MI CCS...\n")
  miCCS <- cor(mi[[1]],mi[[2]])
  
  mySaveRDS(miCCS,"CCS_mi")

  mySaveRDS(calcRanks(miCCS), "ranks_mi")

  
  #
  # calculate CLR
  # 
  
  cat("Calculating CLR...\n")
  clr <- lapply(mi,calcCLR)
  
  #
  # Generate and save CCS and ranks using MI+CLR
  #
  
  cat("Calculating MI+CLR CCS...\n")
  clrCCS <- cor(clr[[1]],clr[[2]])
  
  mySaveRDS(clrCCS,"CCS_clr")

  mySaveRDS(calcRanks(clrCCS), "ranks_clr")

}



#### job step 4 (parallel execution for each subset pair)
#
# Calculate PCC+CCS+ranks
#
PODC_Os_HalfSubsetsJob4 <- function(subset_idx, 
                                    outDir = "data/subsets/PODC_Os_HalfSubsets"){
  #
  # define functions
  #
  
  mySaveRDS <- function(x, prefix){
    filename <- file.path(outDir,
                          paste0(prefix,"_",subset_idx,".RDS"))
    cat("Saving",filename,"...\n")
    saveRDS(x,file = filename) 
  }
  
  
  #
  # load subset definition and expression matrix
  #
  
  subsetsSamples <- readRDS(file.path(outDir,"subsetsSamples.RDS"))
  expMat11 <- readRDS(file.path(outDir,"expMat11.RDS"))
  
  # get logical vector for current subset
  isInSubset <- colnames(expMat11) %in% subsetsSamples[[subset_idx]]

  cat("Calculating PCC for subset",subset_idx,"...\n")  

  pcc1 <- cor(t(expMat11[ ,isInSubset]))

  cat("Calculating PCC for inverse subset",subset_idx,"...\n")  

  pcc2 <- cor(t(expMat11[ ,!isInSubset]))
  
  cat("Calculating PCC CCS...\n")
  pcc_CCS <- cor(pcc1, pcc2)
  
  mySaveRDS(pcc_CCS,"CCS_pcc")
  
  mySaveRDS(calcRanks(pcc_CCS), "ranks_pcc") 
  
}

####
#
# Make At version of the jobs by overriding the default arguments
#
PODC_At_HalfSubsetsJob1 <- function(outDir = "data/subsets/PODC_At_HalfSubsets",
                                    spc="At",
                                    expMatFile="data/expMat/PODC_At.RDS",
                                    PODC_metadata_file = "indata/arabidopsis_thaliana_sample.xlsx"){
  PODC_Os_HalfSubsetsJob1(outDir, spc, expMatFile, PODC_metadata_file )
}

PODC_At_HalfSubsetsJob2 <- function(subset_idx, invertSubset, 
                                    outDir = "data/subsets/PODC_At_HalfSubsets") {
  PODC_Os_HalfSubsetsJob2(subset_idx, invertSubset, outDir)  
}

PODC_At_HalfSubsetsJob3 <- function(subset_idx,
                                    outDir = "data/subsets/PODC_At_HalfSubsets"){
  PODC_Os_HalfSubsetsJob3(subset_idx, outDir )
}

PODC_At_HalfSubsetsJob4 <- function(subset_idx, 
                                    outDir = "data/subsets/PODC_At_HalfSubsets"){
  PODC_Os_HalfSubsetsJob4(subset_idx, outDir)
}


makePODC_Os_Subsets_flow <- function(){
  flowbind(
    # PODC_Os_HalfSubsetsJob1
    Rjob(source = "Rjobs/PODC_OsAt_Subsets.R",
         fun = "PODC_Os_HalfSubsetsJob1" ),
    
    # PODC_Os_HalfSubsetsJob2
    Rjob(prev_jobs = "PODC_Os_HalfSubsetsJob1", dep_type = "burst", sub_type = "scatter",
         source = "Rjobs/PODC_OsAt_Subsets.R",
         fun = "PODC_Os_HalfSubsetsJob2",
         paramMat = expand.grid(subset_idx = 1:10, invertSubset = c(T,F)) ),
    
    # PODC_Os_HalfSubsetsJob3
    Rjob(prev_jobs = "PODC_Os_HalfSubsetsJob2", dep_type = "gather", sub_type = "scatter",
         source = "Rjobs/PODC_OsAt_Subsets.R",
         fun = "PODC_Os_HalfSubsetsJob3",
         paramMat = data.frame(subset_idx = 1:10)),
    
    # PODC_Os_HalfSubsetsJob4
    Rjob(prev_jobs = "PODC_Os_HalfSubsetsJob1", dep_type = "burst", sub_type = "scatter",
         source = "Rjobs/PODC_OsAt_Subsets.R",
         fun = "PODC_Os_HalfSubsetsJob4",
         paramMat = data.frame(subset_idx = 1:10))
  )
}
              
makePODC_At_Subsets_flow <- function(){
  flowbind(
    # PODC_At_HalfSubsetsJob1
    Rjob(source = "Rjobs/PODC_OsAt_Subsets.R",
         fun = "PODC_At_HalfSubsetsJob1" ),
    
    # PODC_At_HalfSubsetsJob2
    Rjob(prev_jobs = "PODC_At_HalfSubsetsJob1", dep_type = "burst", sub_type = "scatter",
         source = "Rjobs/PODC_OsAt_Subsets.R",
         fun = "PODC_At_HalfSubsetsJob2",
         paramMat = expand.grid(subset_idx = 1:10, invertSubset = c(T,F)) ),
    
    # PODC_At_HalfSubsetsJob3
    Rjob(prev_jobs = "PODC_At_HalfSubsetsJob2", dep_type = "gather", sub_type = "scatter",
         source = "Rjobs/PODC_OsAt_Subsets.R",
         fun = "PODC_At_HalfSubsetsJob3",
         paramMat = data.frame(subset_idx = 1:10)),
    
    # PODC_At_HalfSubsetsJob4
    Rjob(prev_jobs = "PODC_At_HalfSubsetsJob1", dep_type = "burst", sub_type = "scatter",
         source = "Rjobs/PODC_OsAt_Subsets.R",
         fun = "PODC_At_HalfSubsetsJob4",
         paramMat = data.frame(subset_idx = 1:10))
  )
}

