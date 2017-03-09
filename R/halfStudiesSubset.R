source("R/PODCfiles.R")

library(readxl)
library(purrr)
library(readr)
library(dplyr)

# species
spcs <- c("At","Gm","Os","Sl","Zm") %>% setNames(.,.)


#' Randomly divide samples per species in half by studyID
#'
#' @param reps number of replicates
#' @param seed random seed
#'
#' @return list with each species containg a list with studyIDs, sampleIDs and subset matrix
generateHalfStudiesSubsets <- function(reps, seed = 700){
  
  
  # set random seed
  set.seed(seed) 
  
  # for  each species
  lapply(spcs, function(spc){
    
    # read sample meta data
    metadata <- loadSampleMeta(spc)

    # get studyIDs
    studyID <- unique(metadata$Study)
    
    # get number of studies
    n <- length(studyID)
    
    # create a n x r matrix of logicals such that each column has 50% TRUE's
    m <- replicate( reps, sample(rep(c(T,F),ceiling(n/2)), n))
    
    rownames(m) <- studyID
    
    return( list( studyID = metadata$Study,
                  sampleID = metadata$Run,
                  subsetMat = m ) )
  })
}





# Select the most relevant subset of genes for each species
getSubsetGeneIDsAndRefOrthos <- function(numberOfGenes, hasExpGeneIDs){
  
  # load orthology data and filter by genes with expression data
  orthos <- 
    dir("indata/orthologs",full.names = T) %>%
    set_names(sub(".*(.._..)_orthologs.txt","\\1",.)) %>%
    map(read_tsv, col_types = cols(), progress = F) %>% 
    # only keep orthologs with expression data:
    map( ~ .x[.x[[1]] %in% hasExpGeneIDs & .x[[2]] %in% hasExpGeneIDs, ]) 

  # helper function to get the ortho set
  getOrthos <- function(spc1,spc2){
    if( spc1 < spc2 ){
      return(orthos[[paste0(spc1,"_",spc2)]])
    } else {
      return(orthos[[paste0(spc2,"_",spc1)]])
    }
  }
  
  # prioritize genes that has 1:1 ortholog in any species and any type of ortho in all species
  
  # for each species
  lapply(spcs, function(spc1){
    
    # the other species
    spcs2 <- spcs[spcs != spc1]
    
    # sort genes by how many species it has 1:1 ortholog in
    prioritize11 <-
      map(spcs2,getOrthos,spc1) %>% # get orthologs to all other species
      reduce( bind_rows ) %>%       # combine to one table
      filter(otype=="1:1") %>%      # keep only 1:1 orthologs
      .[[spc1]] %>%                 # select the spc1 geneIDs
      table() %>%                   # count occurences
      sort(decreasing = T)          # sort by number of occurences
    
    table(prioritize11)
    
    # not enough genes with 1:1?
    if( length(prioritize11) < numberOfGenes){
      # sort genes by how many species it has 1:N ortholog in
      prioritize1N <-
        map( spcs2, getOrthos, spc1) %>% 
        map( filter, otype=="1:N" ) %>% 
        map( ~ unique(.x[[spc1]]) ) %>% 
        reduce( c ) %>%
        .[!(. %in% names(prioritize11))] %>% 
        table() %>% 
        sort(decreasing = T)
      
      return(c(names(prioritize11),names(prioritize1N))[1:numberOfGenes])
    } else {
      return(names(prioritize11)[1:numberOfGenes])
    }
  }) -> subsetGeneIDs
  
  # get the 1:1 ref.orthos but limited to the subset of genes
  lapply(spcs, function(spc1){
    
    # the other species
    spcs2 <- spcs[spcs != spc1]
    
    # get 1:1 orthos
    map(spcs2,getOrthos,spc1) %>%  # get orthologs to all other species
      map( filter, otype=="1:1") %>% # only 1:1 orthos
      map( ~ .x[[spc1]]) %>%         # keep only the geneIDs of the spc of interrest
      map( ~ .x[.x %in% subsetGeneIDs[[spc1]]])  # filter only genes In the subset
  }) -> refOrthos
  
  return( list( subsetGeneIDs = subsetGeneIDs,
                refOrthos = refOrthos))
}

