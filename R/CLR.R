# # Note: This function does not support matrices larger than 2^31 and can use a lot of memory
# calcCLR <- function(mi){
#   miZ <- pmax(scale(mi),0,na.rm = T)^2
#   sqrt(miZ + t(miZ))
# }

# memory efficient implementation, supports matrices larger than 2^31
calcCLR <- function(mi){
  n <- ncol(mi)
  
  # calculate sum and standard deviation for each column
  miColMean <- numeric(n)
  miColSD <- numeric(n)
  # for each column
  for(i in 1:n){
    miColMean[i] <- mean(mi[ ,i],na.rm = T)
    miColSD[i] <- sd(mi[ ,i],na.rm = T)
  }
  
  # create clr matrix with same dimensions and dimnames as mi
  clr <- mi
  
  for(i in 1:n){
    # note: clr[ ,i] is the same as mi[refOrthoIDs,i]
    Zc <- (clr[ ,i] - miColMean[i])/miColSD[i]
    Zr <- (clr[ ,i] - miColMean)/miColSD
    clr[ ,i] <- sqrt(pmax(Zc,0,na.rm = T)^2 + pmax(Zr,0,na.rm = T)^2)
  }
  
  return(clr)
}


# memory efficient implementation
# returns only the rows given by refOrthoIDs
calcCLRref <- function(mi, refOrthoIDs){
  n <- ncol(mi)

  # calculate sum and standard deviation for each column
  miColMean <- numeric(n)
  miColSD <- numeric(n)
  # for each column
  for(i in 1:n){
    miColMean[i] <- mean(mi[ ,i],na.rm = T)
    miColSD[i] <- sd(mi[ ,i],na.rm = T)
  }
  
  # only need the rows with ref.orthos
  clr <- mi[refOrthoIDs, ]
  miColMeanRef <- miColMean[match(refOrthoIDs,rownames(mi))]
  miColSDRef <- miColSD[match(refOrthoIDs,rownames(mi))]

  for(i in 1:n){
    # note: clr[ ,i] is the same as mi[refOrthoIDs,i]
    Zc <- (clr[ ,i] - miColMean[i])/miColSD[i]
    Zr <- (clr[ ,i] - miColMeanRef)/miColSDRef
    clr[ ,i] <- sqrt(pmax(Zc,0,na.rm = T)^2 + pmax(Zr,0,na.rm = T)^2)
  }
  
  return(clr)
}

