myLog <- function(...){
  cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),...)
}




# Generate random subsets and prepare Os MR big.matrix and At expression 
# big.matrix for quick retrieval by subsequent steps
nSamplesInit <- function(outdir = "data/subsets/nSamples", nStudyReps=10, nSampleReps=5){
  source("R/orthoUtils.R")
  source("R/calcMR.R")
  library(bigmemory)
  library(dplyr)

  # create output directory
  dir.create(outdir)
  
  # load expression data
  xAt <- readRDS("data/expMat/PODC_At.RDS")
  xOs <- readRDS("data/expMat/PODC_Os.RDS")
  
  # get As:Os 1:1 orthos
  AtOsOrthos <- 
    loadOrthoTable("At","Os") %>% 
    filter(otype == "1:1") %>%
    filter(At %in% rownames(xAt)) %>% 
    filter(Os %in% rownames(xOs))
  
  # only keep expression of 1:1 orthos
  xOs <- xOs[AtOsOrthos$Os, ]
  xAt <- xAt[AtOsOrthos$At, ]
  
  # save At expression as big.matrix
  xAt.bm <- as.big.matrix(xAt, backingpath = outdir, 
                          backingfile = "xAt.bm", 
                          descriptorfile = "xAt.bm.desc")
  
  # calc Os PCC+MR and save as big.matrix
  cat("calculating PCC+MR for Os...\n")
  mrOs <- calcLogMR(WGCNA::cor(t(xOs)))
  
  cat("Writing MR matrix to file.\n")
  mrOs.bm <- as.big.matrix(mrOs, backingpath = outdir, 
                            backingfile = "mrOs.bm", 
                            descriptorfile = "mrOs.bm.desc")
  
  
  # load sample metadata
  mAt <- readRDS("indata/sampleMeta/At.RDS")
  
  # set seed to be reproducable
  set.seed(3)
  
  # group samples by study
  studySamples <- split(mAt$Run, mAt$Study)
  
  # randomize studies
  rndStudies <- replicate(simplify = F, n = nStudyReps, {
    unlist(sample(studySamples))
  })
  
  saveRDS( rndStudies, file.path(outdir,"rndStudies.RDS") )
  
  # set seed to be reproducable
  set.seed(8135)
  
  # randomize samples
  rndSamples <- replicate(simplify = F, n = nSampleReps, {
    sample(mAt$Run)
  })
  
  saveRDS( rndSamples, file.path(outdir,"rndSamples.RDS") )
  
}

nSamplesJob <- function(n = 100, 
                        subsetIdx = 1, 
                        subsetFile = "data/subsets/nSamples/rndStudies.RDS", 
                        outfile = "data/subsets/nSamples/rnks_rndStudies_1_100.RDS", 
                        mrOs.bm.desc = "data/subsets/nSamples/mrOs.bm.desc", 
                        xAt.bm.desc = "data/subsets/nSamples/xAt.bm.desc"){
  library(bigmemory)
  source("R/calcMR.R")
  source("R/calcRanks.R")

  
  myLog("calculating PCC+MR for first",n,"At samples of subset",subsetIdx,"in",subsetFile,"...\n")
  
  # attach At expression data
  xAt.bm <- attach.big.matrix(xAt.bm.desc)

  sampleIDs <- readRDS(subsetFile)[[subsetIdx]][1:n]
  
  # extract subset
  xAtSubset <- xAt.bm[ ,sampleIDs]
  
  # find genes with no variance
  noVar <- apply(xAtSubset,1,max)==apply(xAtSubset,1,min) # max=min -> no variance!
  
  myLog("Number of genes with no variance:",sum(noVar),"\n")

  # calculate PCC+MR (filter noVar genes)
  mrAt <- calcLogMR(WGCNA::cor(t(xAtSubset[!noVar, ])))

  myLog("calculating CCS...\n")

  # attach Os MR matrix
  mrOs.bm <- attach.big.matrix(mrOs.bm.desc)
  
  # extract submatrix without the noVar genes
  mrOs <- mrOs.bm[!noVar,!noVar]
  
  # calc CCS
  CCS <- calcCCS(mr1 = mrAt, mr2 = mrOs,
                 refOrthos1 = rownames(mrAt),
                 refOrthos2 = rownames(mrOs))
  
  myLog("calculating ranks...\n")
  
  rnks <- list( 
    AtOs = quickRanks(CCS,rownames(CCS),colnames(CCS)),
    OsAt = quickRanksT(CCS,colnames(CCS),rownames(CCS)))
  
  myLog("writing to file",outfile,"\n")
  
  saveRDS(rnks,outfile)
}