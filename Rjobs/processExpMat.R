####
#
# processExpMatConvertIDs
#
# Convert transcript expression matrix to gene expression matrix by summing 
# isoform expression levels. 
# Convert transcript IDs to gene IDs.
# Log transform
#
processExpMatConvertIDs <- function( inFile, outFile, convertID ){
  
  cat("Reading", inFile,"...\n")
  srcExp <- as.matrix(read.table(inFile,header=T,row.names = 1))
  
  
  cat("Processing...\n")
  
  # convert transcript IDs to gene IDs
  srcGenes <- convertID(transcriptIDs = rownames(srcExp))
  dstGenes <- unique(srcGenes) # keep one copy of eachall genes
  
  # add together the FPKM values of all isoforms for each gene
  idx <- match(dstGenes, srcGenes)
  dstExp <- srcExp[idx, ]
  rownames(dstExp) <- dstGenes
  
  repeat{
    # remove the transcripts we have already added
    srcGenes[idx] <- ""
    # find next matching genes
    idx <- match(dstGenes, srcGenes)
    isNotNA <- !is.na(idx)
    idx <- idx[isNotNA]
    # stop if nothing left
    if(length(idx)==0) break
    # add isoform expression
    dstExp[isNotNA, ] <- dstExp[isNotNA, ] + srcExp[idx, ]
  }
  
  # log transform
  dstExp <- log2(0.1+dstExp)
  
  # Make directory if it doesn't exist
  if( !file.exists(dirname(outFile)) ){
    dir.create(dirname(outFile),recursive = T, showWarnings = F)
  }
  
  cat("Saving", outFile, "...\n")
  saveRDS(dstExp,file = outFile)
}


processExpMat_ebi <- function( spc ){
  library(readr)
  source("R/PODCfiles.R")
  
  inPath <- paste0("indata/RNAseq",spc)
  outFile <- paste0("data/expMat/EBI_",spc,".RDS")
  sampleMeta <- loadSampleMeta(spc)
  
  inFiles <- dir(inPath,pattern = "genes.fpkm.tsv",full.names = T)
  expDFs <- lapply(inFiles,read_tsv)
  
  # convert data.frame to matrix using first column as rownames
  lapply(expDFs,function(df){
    m <- as.matrix(df[,-1])
    rownames(m) <- df[[1]]
    return(m)
  }) -> expMats
  
  expMat <- do.call(cbind,expMats)
  
  # sanity check the rownames just in case:
  invisible(
    sapply(expMats,function(x){
      if(!identical(rownames(x),rownames(expMat)))
        stop("expression matrix rownames don't match!")
    }))

  # reorder and remove columns to match metadata
  expMat <- expMat[ , match(sampleMeta$Run, colnames(expMat)) ]
  
  # log transform
  expMat <- log2(0.1+expMat)
  
  # Make directory if it doesn't exist
  if( !file.exists(dirname(outFile)) ){
    dir.create(dirname(outFile),recursive = T, showWarnings = F)
  }
  
  cat("Saving", outFile, "...\n")
  saveRDS(expMat,file = outFile)  
}

####
#
# Convert transcript IDs to gene IDs
# 

convertID <- list(
  # (Oryza sativa) E.g. Os03t0413100-03 -> OS03G0413100
  Os = function(transcriptIDs){
    sub("Os([0-9]+)t([0-9]+)-[0-9]+$","OS\\1G\\2",transcriptIDs) 
  },
  
  # (Arabidopsis thaliana) E.g. AT2G37860.2 -> AT2G37860
  At = function(transcriptIDs){
    sub("\\.[0-9]+$","",transcriptIDs)
  },
  
  # (Solanum lycopersicum) E.g Solyc05g041450.2.1 -> Solyc05g041450.2
  Sl = function(transcriptIDs){
    sub("\\.[0-9]+$","",transcriptIDs)
  }  
)

spc2PODCspcNames <- 
  c(Os = "oryza_sativa",
    At = "arabidopsis_thaliana",
    Sl = "solanum_lycopersicum",
    Sb = "sorghum_bicolor",
    Vv = "vitis_vinifera",
    Mt = "medicago_truncatula",
    St = "solanum_tuberosum",
    Gm = "glycine_max",
    Nt = "nicotiana_tabacum",
    Zm = "zea_mays")

processExpMat_PODC <- function(spc){
  processExpMatConvertIDs( 
    inFile = paste0("indata/",spc2PODCspcNames[spc],"_fpkm.tsv.gz"),
    outFile = paste0("data/expMat/PODC_",spc,".RDS"),
    convertID = convertID[[spc]])
}


# # check the transcript conversion
# 
# x <- system(paste("gunzip -c ","indata/solanum_lycopersicum_fpkm.tsv.gz"," | cut -f 1"),inter=T)[-1]
# 
# Genes <- biomaRt::getBM(
#   attributes=c("ensembl_gene_id", "ensembl_transcript_id"),
#   mart = biomaRt::useMart("plants_mart",dataset="slycopersicum_eg_gene", 
#                           host="plants.ensembl.org"))
# 
# table(x %in% SlGenes$ensembl_transcript_id)
# table(SlGenes$ensembl_transcript_id %in% x)
# 
# z <- convertID$Sl(x)



