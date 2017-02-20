source("R/phyloXML.R")

#' Extract the information of interrest from all the phyloXML files downloaded from Ensembl
#'
#' @param arrayIdx the index of this specific process
#' @param arraySize the number of parallel processes
#' @param outDir directory to write output to
#'
#' @return None
extractTreeDataJob <- function(arrayIdx, arraySize, outDir){
  # list extracted files
  # phyloXMLfiles <- dir("data/Compara.phyloxml_aa_trees.33/016", recursive = T, full.names = T) # test with only 1000 files
  phyloXMLfiles <- dir("data/Compara.phyloxml_aa_trees.33", recursive = T, full.names = T)
  # get name of tree from filename
  names(phyloXMLfiles) <- sub("\\.aa.xml$","",basename(phyloXMLfiles))
  
  idx <- seq(arrayIdx,length(phyloXMLfiles),by = arraySize)
  tl <- lapply(phyloXMLfiles[idx], function(phyloXMLfile){
    cat("Reading",phyloXMLfile,"...\n")
    extractTreeData(loadPhyloXMLasList(phyloXMLfile))
  })
  
  outFile <- file.path(outDir,sprintf("treeData%04i.RDS",arrayIdx))
  cat("Writing",outFile,"\n")
  saveRDS(tl,file=outFile)
}

#' Merge the files output from extractTreeDataJob
#'
#' @param outDir directory where the output from extractTreeDataJob can be found
#'
#' @return None
mergeTreeDataJob <- function(outDir){
  inFiles <- dir(outDir,pattern="treeData[0-9]+.RDS",full.names = T)
  cat("Loading",length(inFiles),"files...\n")
  mergedTreeData <- do.call(c,lapply(inFiles,readRDS))
  outFile <- file.path(outDir,"treeData.RDS")
  cat("Writing",outFile,"\n")
  saveRDS(mergedTreeData,file=outFile)
  cat("Removing",length(inFiles),"temporary files...\n")
  unlink(inFiles)
}

# Command to extract the downloaded tar file:
# tar -zxvf --directory=data indata/Compara.phyloxml_aa_trees.33.tar.gz 

#
# command line:
# cd $PROJECTDIR && Rscript  -e "source('Rjobs/extractTreeDataJob.R'); extractTreeDataJob( $ARRAY_INDEX, $ARRAY_SIZE, 'data/treeData')"