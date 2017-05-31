
# calculate log transformed MR from correlation matrix
calcLogMR <- function(pcc){
  #   MR = sqrt( r(i,j)*r(j,i) )
  #   logMR = 0.5*( log2(r(i,j)) + log2(r(j,i)) )
  
  rankMat <- log2(apply(-pcc,2,rank))
  
  nGenes <- nrow(pcc)
  
  return( log2(nGenes) - 0.5*(rankMat+t(rankMat)) )
  
  # note: adding log2(nGenes) and multiplying with -0.5 has no effect on CCS
  # it simply changes the scale to something similar to CLR
}

# calculate CCS with removal the self co-expression in the MR matrix
calcCCS <- function( mr1, mr2, refOrthos1,refOrthos2){
  # Set self co-expression values to mean to ignore it from correlation
  # This is similar to setting it to NA but the following correlation is faster
  mapply(SIMPLIFY = F, mr=list(mr1,mr2),refOrthos=list(refOrthos1,refOrthos2), FUN=function(mr,refOrthos){
    m <- mr[refOrthos, ]
    rowi <- match(colnames(m),rownames(m))
    for(i in which(!is.na(rowi))){
      m[rowi[i],i] <- mean(m[-rowi[i],i])
    }
    return(m)
  }) -> M
  
  # calculate CCS
  return( WGCNA::cor(M[[1]],M[[2]]) )
}