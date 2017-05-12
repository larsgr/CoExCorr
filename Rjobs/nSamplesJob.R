myLog <- function(...){
  cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),...)
}


# Generate random subsets and prepare Os clr big.matrix and At expression 
# big.matrix for quick retrieval by subsequent steps
nSamplesInit <- function(outdir = "data/subsets/nSamples", nStudyReps=10, nSampleReps=5){
  source("R/orthoUtils.R")
  library(bigmemory)
  library(dplyr)
  library(BSplineMI)
  source("R/CLR.R")
  
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
  
  # calc Os CLR and save as big.matrix
  cat("calculating MI+CLR for Os...\n")
  clrOs <- calcCLR(calcSplineMI(xOs,nBins = 7,splineOrder = 3))
  
  cat("Writing clr matrix to file.\n")
  clrOs.bm <- as.big.matrix(clrOs, backingpath = outdir, 
                            backingfile = "clrOs.bm", 
                            descriptorfile = "clrOs.bm.desc")
  
  
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
                        clrOs.bm.desc = "data/subsets/nSamples/clrOs.bm.desc", 
                        xAt.bm.desc = "data/subsets/nSamples/xAt.bm.desc",
                        threads = 4){
  library(bigmemory)
  library(BSplineMI)
  source("R/CLR.R")
  source("R/mc_cor.R")
  source("R/calcRanks.R")

  
  myLog("calculating MI+CLR for first",n,"At samples of subset",subsetIdx,"in",subsetFile,"...\n")
  
  # attach At expression data
  xAt.bm <- attach.big.matrix(xAt.bm.desc)

  sampleIDs <- readRDS(subsetFile)[[subsetIdx]][1:n]
  
  clrAt <- calcCLR(calcSplineMI(xAt.bm[ ,sampleIDs],nBins = 7,splineOrder = 3,threads=threads))
  
  myLog("calculating CCS...\n")

  clrOs.bm <- attach.big.matrix(clrOs.bm.desc)
  CCS <- mc_cor2(clrAt,clrOs.bm,cores = threads)
  
  myLog("calculating ranks...\n")
  
  rnks <- list( 
    AtOs = quickRanks(CCS,rownames(CCS),colnames(CCS)),
    OsAt = quickRanksT(CCS,colnames(CCS),rownames(CCS)))
  
  myLog("writing to file",outfile,"\n")
  
  saveRDS(rnks,outfile)
}