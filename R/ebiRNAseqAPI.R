library(readr)
library(jsonlite)

RNAseq_getRun <- function( runID, minMappingQuality = 70 ){
  url <- paste0("http://www.ebi.ac.uk/fg/rnaseq/api/json/",minMappingQuality,"/getRun/",runID)
  fromJSON(url)
}

RNAseq_getRunsByStudy <- function( studyID, minMappingQuality = 70){
  url <- paste0("http://www.ebi.ac.uk/fg/rnaseq/api/tsv/",minMappingQuality,"/getRunsByStudy/",studyID)
  read_tsv(file = url, col_types = cols())
}

RNAseq_getStudiesByOrganism <- function( organism ){
  url <- paste0("http://www.ebi.ac.uk/fg/rnaseq/api/tsv/getStudiesByOrganism/",organism)
  read_tsv(file = url, col_types = cols())
}


