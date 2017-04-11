

# /* "Lifted" from http://local.wasp.uwa.edu.au/~pbourke/curves/spline/ */
#   /*
#   Calculate the blending value, this is done recursively.
# 
# If the numerator and denominator are 0 the expression is 0.
# If the denominator is 0, the expression is 0
# SplineBlend(curBin, splineOrder, knots, z[curSample], numBins);
# http://en.wikipedia.org/wiki/B-spline
# Computes the basis B-spline of degree t and knot u[k] using the Cox-deBoor recursion formula
# */
sourceCpp("R/BSplineMI.cpp")


calcWeights <-function( x, numBins, splineOrder){
  
  # scales the expressions for each gene to between [0,1]*(numBins-splineOrder+1)
  z <- apply(x,1,scale0to1) * (numBins-splineOrder+1)
  # z is transposed of x
  
  knots <- calcKnots(numBins, splineOrder)
  
  w <- SplineBlendAll(z, knots, splineOrder, numBins)
  dim(w) <- c(numBins,dim(z))
  dimnames(w) <- c(list(paste0("bin",1:numBins)),dimnames(z))
  
  return(w)
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

# calculate entropy between two genes 
entropy2dR <- function( w1, w2 ){
  # p is mean of weights over all samples per 2d bin
  p <- rowMeans(sapply(1:ncol(w1),function(k){w1[,k] %*% t(w2[,k])}))
  # p: array [numBins^2]
  
  # calculate entropy
  -sum(p * log2(p), na.rm = T)
}


entropy2d <- function( w1, w2 ){
  # p is mean of weights over all samples per 2d bin
  p <- hist2d(w1, w2)
  # p: array [numBins x genes]
  
  # calculate entropy
  -sum(p * log2(p), na.rm = T)
}


calcKnots <- function(numBins, splineOrder){
  c(rep(0,splineOrder),1:(numBins-splineOrder),rep(numBins-splineOrder+1,splineOrder))
}

# SplineBlend(curBin, splineOrder, knots, z[curSample], numBins);
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

# x: vector of expression values for one gene
scale0to1 <- function(x){
  #scales the expressions to between [0,1]*(numBins-splineOrder+1)
  Xrange <- range(x,na.rm = T)
  Xmin <- Xrange[1]
  Xfactor <- (Xrange[2]-Xrange[1])
  return((x-Xmin)/Xfactor)
}

# main function
calcSplineMI <- function(x, numBins, splineOrder){

  system.time(
    weights <- calcWeights(x, numBins, splineOrder)
  )
  entropy <- entropy1d(weights)

  nGenes <- nrow(x)
  mi <- matrix(0,nrow=nGenes,ncol=nGenes)
  for( i in 2:nGenes){
    for( j in 1:(i-1)){
      mi[i,j] <- mi[j,i] <- entropy[j] + entropy[i] - entropy2d(weights[,,i],weights[,,j])
    }
  }
  return(mi)
}

## Testing ####
#
#

x = eOs[1:50, ]

# x = eOs[c(1,2,957,1910), 1:100]


w <- calcWeights(x,numBins = 7,splineOrder = 3)

w1 <- w[,,1]
w2 <- w[,,2]
p <- rowMeans(sapply(1:ncol(w1),function(k){w1[,k] %*% t(w2[,k])}))
p
p2 <- hist2d(w1,w2)
dim(p) <- dim(p2)

# plot bins
x = matrix(0:100,nrow=1)
w <- calcWeights(x,numBins = 7,splineOrder = 3)
plot(NULL,xlim=c(0,100),ylim=c(0,1),xlab="x",ylab="weight")
for(i in 1:dim(w)[1]){
  lines(x[1, ],w[i, ,1],col=rainbow(dim(w)[1])[i])
}

source("Rjobs/calcMI.R")
writeTempExpMat(x,filename = "~/temp/test.expmat")
unlink("~/temp/test.mi")
system.time(
  system(genepairMICmd( tmpExpMatFile = "~/temp/test.expmat",
                        miFile = "~/temp/test.mi",
                        ngenes = nrow(x), nexp = ncol(x)))
)
source("R/loadTriMatrix.R")
mi_0 <- loadTriMatrix(rownames(x),filename = "~/temp/test.mi")

system.time(
  mi <- calcSplineMI(x,numBins = 7,splineOrder = 3)
)
range(mi-mi_0)
