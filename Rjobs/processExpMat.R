####
#
# processExpMat_PODC
#
# Convert transcript expression matrix to gene expression matrix by summing 
# isoform expression levels. 
# Convert transcript IDs to gene IDs.
# Log transform
#
processExpMat_PODC <- function( inFile, outFile, convertID ){
  
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
####
#
# Convert transcript IDs to gene IDs
# 

# (Oryza sativa) E.g. Os03t0413100-03 -> OS03G0413100
convertID_Os <- function(transcriptIDs){
  sub("Os([0-9]+)t([0-9]+)-[0-9]+$","OS\\1G\\2",transcriptIDs) 
}

# (Arabidopsis thaliana) E.g. AT2G37860.2 -> AT2G37860
convertID_At <- function(transcriptIDs){
  sub("\\.[0-9]+$","",transcriptIDs)
}


processExpMat_PODC_Os <- function(){
  processExpMat_PODC( inFile = "indata/oryza_sativa_fpkm.tsv.gz",
                      outFile = "data/expMat/PODC_Os.RDS",
                      convertID = convertID_Os)
}

processExpMat_PODC_At <- function(){
  processExpMat_PODC( inFile = "indata/arabidopsis_thaliana_fpkm.tsv.gz",
                      outFile = "data/expMat/PODC_At.RDS",
                      convertID = convertID_At)
}
