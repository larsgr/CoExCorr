# Utility functions for handling orthology data from Ensembl biomart
library("readr") # read_tsv

loadOrthoTable <- function(spc1,spc2){
  read_tsv(paste0("indata/orthologs/",min(spc1,spc2),"_",max(spc1,spc2),"_orthologs.txt"),
                  col_types = cols(), progress = F)
}


####
#
# getOrthos1N(orthoTable, spc1, spcN)
#
# Get 1:N orthologs as a named list (spc1->spc2)
#
getOrthos1N <- function(orthoTable, spc1, spcN){
  # get 1:N (and N:1) orthos using the otype
  ort1N <- orthoTable[orthoTable$otype == "1:N",c(spc1,spcN)]
  
  # group by the species with single ortholog
  grpSpc1 <- tapply(ort1N[[spcN]],ort1N[[spc1]],c)
  
  # remove the groups that are actually N:1
  return( grpSpc1[sapply(grpSpc1,length)>1] )
}

####
#
# getOrthos1Nmatrix(orthoTable, N, spc1, spcN, spc1_geneIDs, spcN_geneIDs)
#
# Get 1:N orthos where N is specified and filter orthologs based on a list of 
# geneIDs. Return as matrix with N columns.
#
getOrthos1Nmatrix <- function(orthoTable, N, spc1, spcN, 
                              spc1_geneIDs, spcN_geneIDs){
  # get 1:N orthos for any N
  ortho1N <- getOrthos1N(orthoTable, spc1, spcN)
  
  # filter specific N
  ortho1N <- ortho1N[sapply(ortho1N,length)==N]
  
  # convert to matrix
  ortho1N <- matrix(unlist(ortho1N),ncol=N,byrow = T,dimnames=list(names(ortho1N)))
  
  # filter genes
  idxGeneMissing <- which(!(names(ortho1N) %in% spc1_geneIDs))
  idxOrthoMissing <- which(!(apply( apply(ortho1N,2,function(x){x %in% spcN_geneIDs}),1,all)))
  
  idxMissing <- unique(c(idxGeneMissing,idxOrthoMissing))
  if(length(idxMissing) > 0){
    cat("Removing",length(idxMissing),"ortholog groups because of unknown geneIDs.\n")
    ortho1N <- ortho1N[-idxMissing, ]
  }
  
  return(ortho1N)
}

orthoTbl2Array <- function(orthoTable){
  # first group by the genes in spc1 and collapse the spc2 orthologs to a single string
  grpOrt2byTbl1 <- tapply(orthoTable[[2]],orthoTable[[1]], function(x){paste(sort(x),collapse=",")})
  # Then group by the collapsed spc2 orthologs and collapse the spc2 genes which have these orthologs
  # (It is assumed that the Gene IDs were sorted)
  grpOrt1byOrt2 <- tapply(names(grpOrt2byTbl1),grpOrt2byTbl1,paste,collapse=",")
  
  # convert to array of list
  grps <- strsplit(c(as.vector(grpOrt1byOrt2),names(grpOrt1byOrt2)),split=",")
  dim(grps) <- c(length(grpOrt1byOrt2),2)
  colnames(grps) <- colnames(orthoTable)[1:2]
  
  return(grps)
}

####
#
# getOrthosNN(orthoTable, spc1, spc2)
#
# Get N:N orthologs as array of list of geneIDs
#
getOrthosNN <- function(orthoTable, spc1, spc2){
  # get N:N orthos table
  orthoTbl2Array(orthoTable[orthoTable$otype == "N:N",c(spc1,spc2)])
}
