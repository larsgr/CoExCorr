library(Rcpp)
sourceCpp("R/BSplineMI.cpp")


calcWeights <-function( x, numBins, splineOrder){
  
  # scales the expressions for each gene to between [0,1]*(numBins-splineOrder+1)
  z <- apply(x,1,scale0to1) * (numBins-splineOrder+1)
  # note z is transposed. i.e. [samples x genes]
  
  knots <- calcKnots(numBins, splineOrder)
  
  w <- SplineBlendAll(z, knots, splineOrder, numBins)
  # add some dimensions
  # should be [bins x samples x genes]
  dim(w) <- c(numBins,dim(z))
  dimnames(w) <- c(list(paste0("bin",1:numBins)),dimnames(z))
  
  return(w)
}


# weights: 3d array [numBins x samples x genes]
# returns: vector of entropy per gene in bits
entropy1d <- function( weights ){
  # p is mean of weights for all samples per bin and gene
  P <- apply(weights, c(1,3),mean)
  # p: array [numBins x genes]

  # for each gene calculate entropy
  apply(P,2,function(p){
    -sum(p * log2(p), na.rm = T)
    # note: na.rm = T because when p=0 0*log2(0) returns NaN
  })
}


calcKnots <- function(numBins, splineOrder){
  c(rep(0,splineOrder),1:(numBins-splineOrder),rep(numBins-splineOrder+1,splineOrder))
}

# scales a vector to between [0,1]
scale0to1 <- function(x){
  Xrange <- range(x,na.rm = T)
  Xmin <- Xrange[1]
  Xfactor <- (Xrange[2]-Xrange[1])
  if( Xfactor == 0){
    # handle genes where all samples have 0 expression
    return( rep(0.0,length(x)) )
  } else{
    return( (x-Xmin)/Xfactor )
  }
}

# main function
calcSplineMI <- function(x, numBins, splineOrder){

  weights <- calcWeights(x, numBins, splineOrder)

  entropy <- entropy1d(weights)

  nGenes <- nrow(x)
  mi <- matrix(0,nrow=nGenes,ncol=nGenes)
  for( i in 2:nGenes){
    for( j in 1:(i-1)){
      mi[i,j] <- mi[j,i] <- entropy[j] + entropy[i] - entropy2dC(weights,i-1,j-1)
    }
  }
  return(mi)
}

library(parallel)

bind3d <- function(arrayList){
  X <- do.call(c,arrayList)
  
  dim(X) <- c(dim(arrayList[[1]])[1:2],rowSums(sapply(arrayList,dim))[3])
  dim3names <- do.call(c,lapply(setNames(arrayList,NULL),function(x){dimnames(x)[[3]]}))
  dimnames(X) <- c(dimnames(arrayList[[1]])[1:2],list(dim3names))
  return(X)
}


calcSplineMImultiCore <- function(x, numBins, splineOrder, cores=5){
  
  nChunks <- cores
  
  chunkIdx <- split(1:nrow(x),cut(1:nrow(x),breaks=nChunks,labels = paste0("chunk",1:nChunks)))
  
  weightChunk <- mclapply(chunkIdx, function(idx){
    calcWeights(x[idx, ], numBins, splineOrder)
  },mc.cores = cores)

  entropyChunk <- mclapply(weightChunk, entropy1d ,mc.cores = cores)

  # concatenate the weights and entropies
  weight <- bind3d(weightChunk)
  #entropy <- do.call(c, setNames(entropyChunk,NULL))
  rm(weightChunk)

  system.time({
    
    chunk_i <- 1:nChunks
    chunk_j <- 1:nChunks
    
    for( i in 2:nChunks){
      for( j in 1:(i-1)){
        chunk_i <- c(chunk_i, i)
        chunk_j <- c(chunk_j, j)
      }
    }
    
    mcmapply( c_i=chunk_i, c_j=chunk_j, FUN= function(c_i,c_j){
      if( c_i == c_j ){
        idx <- chunkIdx[[c_i]]
        nGenes <- length(idx)
        entropy <- entropyChunk[[c_i]]
        miChunk <- matrix(0,nrow=nGenes,ncol=nGenes)
        for( i in 2:nGenes){
          for( j in 1:(i-1)){
            miChunk[i,j] <- miChunk[j,i] <- entropy[j] + entropy[i] - entropy2dC(weights,idx[i]-1,idx[j]-1)
          }
        }
      } else {
        idx_i <- chunkIdx[[c_i]]
        idx_j <- chunkIdx[[c_j]]
        miChunk <- matrix(0,nrow=length(idx_i),ncol=length(idx_j))
        for( i in 1:length(idx_i)){
          for( j in 1:length(idx_j)){
            miChunk[i,j] <-  entropyChunk[[c_j]][j] +  entropyChunk[[c_i]][i] - 
                             entropy2dC(weights,idx_i[i]-1,idx_j[j]-1)
          }
        }
      }
      return(miChunk)
    }, mc.cores=cores, SIMPLIFY = F) -> miChunks
    
  })
  
  # put miChunks together
  nGenes <- nrow(x)
  mi <- matrix(0,nrow=nGenes,ncol=nGenes)
  
  for( i in 1:length(chunk_i)){
    c_i <- chunk_i[i]
    c_j <- chunk_j[i]
    mi[chunkIdx[[c_i]],chunkIdx[[c_j]]] <- miChunks[[i]]
    if( c_i != c_j ){
      mi[chunkIdx[[c_j]],chunkIdx[[c_i]]] <- t(miChunks[[i]])
    }
  }
  
  return(mi)
}

## R versions of functions implemented in Rcpp ####
#
#

# calculate entropy between two genes 
entropy2dR <- function( w1, w2 ){
  # p is mean of weights over all samples per 2d bin
  p <- rowMeans(sapply(1:ncol(w1),function(k){w1[,k] %*% t(w2[,k])}))
  # p: array [numBins^2]
  
  # calculate entropy
  -sum(p * log2(p), na.rm = T)
}


SplineBlendR <- function( curBin, splineOrder, knots, v, numBins){
  
  if (splineOrder == 1) {
    
    if (((knots[curBin] <= v) && (v < knots[curBin+1])) ||
        ( abs(v - knots[curBin+1]) < 1e-10 && (curBin == numBins) ) ){
      value = 1.0; 
    } else {
      value = 0.0;
    }
    
  } else {
    
    d1 = knots[curBin+splineOrder-1] - knots[curBin];
    d2 = knots[curBin+splineOrder] - knots[curBin+1];
    
    if ( (d1 == 0) && (d2 == 0) ){
      value = 0;
    } else if (d1 == 0) {
      value = (knots[curBin+splineOrder] - v) / d2 * Recall(curBin+1,splineOrder-1,knots,v,numBins);
    } else if (d2 == 0) {
      value = (v - knots[curBin]) / d1 * Recall(curBin,splineOrder-1,knots,v,numBins);
    } else {
      value = (v - knots[curBin]) / d1 * Recall(curBin,splineOrder-1,knots,v,numBins) +
        (knots[curBin+splineOrder] - v) / d2 * Recall(curBin+1,splineOrder-1,knots,v,numBins);
    }
  }
  if (value < 0)  {
    #  this does not happen...
    value = 0; # rounding sometimes makes this < 0, e.g. -0.000000001 
  }
  return(value);
}


calcWeightsR <-function( x, numBins, splineOrder){
  
  # scales the expressions for each gene to between [0,1]*(numBins-splineOrder+1)
  z <- apply(x,1,scale0to1) * (numBins-splineOrder+1)
  # z is transposed of x
  
  knots <- calcKnots(numBins, splineOrder)
  
  apply(z,c(1,2),function(zValue){
    sapply(1:numBins,function(curBin){
      SplineBlendR(curBin, splineOrder, knots, zValue, numBins);
    })
  })
}


## Testing ####
#
#

eOs <- readRDS("data/expMat/PODC_Os.RDS")

x = eOs[1:100, ]

system.time(
  z <- apply(x,1,scale0to1)
)

x = eOs[c(1,2,957,1910), 1:100]


w <- calcWeights(x,numBins = 7,splineOrder = 3)

w1 <- w[,,1]
w2 <- w[,,3]

w2[,1:10]

entropy2dR(w1,w2)
entropy2dC(w, i1=0, i2=1)



weights <- calcWeights(x, numBins, splineOrder)

entropy <- entropy1d(weights)


nGenes <- nrow(x)
mi <- matrix(0,nrow=nGenes,ncol=nGenes)
system.time({
  for( i in 2:nGenes){
    for( j in 1:(i-1)){
      mi[i,j] <- mi[j,i] <- entropy[j] + entropy[i] - entropy2dR(weights[,,i],weights[,,j])
    }
  }
})

system.time({
  for( i in 2:nGenes){
    for( j in 1:(i-1)){
      mi[i,j] <- mi[j,i] <- entropy[j] + entropy[i] - entropy2dC(weights,i-1,j-1)
    }
  }
})

entropy2dDummy <- function(x,y){
  return(1)
}


# plot bins
x = matrix(0:100,nrow=1)
w <- calcWeights(x,numBins = 7,splineOrder = 3)
plot(NULL,xlim=c(0,100),ylim=c(0,1),xlab="x",ylab="weight")
for(i in 1:dim(w)[1]){
  lines(x[1, ],w[i, ,1],col=rainbow(dim(w)[1])[i])
}


source("Rjobs/calcMI.R")
source("R/loadTriMatrix.R")


x = eOs[1:500, ]

writeTempExpMat(x,filename = "~/temp/test.expmat")
unlink("~/temp/test.mi")
system.time(
  system(genepairMICmd( tmpExpMatFile = "~/temp/test.expmat",
                        miFile = "~/temp/test.mi",
                        ngenes = nrow(x), nexp = ncol(x)))
)
mi_0 <- loadTriMatrix(rownames(x),filename = "~/temp/test.mi")

system.time(
  mi <- calcSplineMI(x,numBins = 7,splineOrder = 3)
)

range(mi-mi_0)

system.time(
  mi_0 <- calcSplineMImultiCore(x,numBins = 7,splineOrder = 3,cores = 20)
)

range(mi[1:5,1:5]-mi_0[1:5,1:5])
