#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//




// SplineBlend(curBin, splineOrder, knots, z[curSample], numBins);

// [[Rcpp::export]]
double SplineBlend( int curBin, int splineOrder, NumericVector knots, double v, int numBins){
  double value = 0.0;
  double d1,d2;
  
  if (splineOrder == 1) {
    
    if (((knots[curBin] <= v) && (v < knots[curBin+1])) ||
        ( abs(v - knots[curBin+1]) < 1e-10 && (curBin+1 == numBins) ) ){
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
      value = (knots[curBin+splineOrder] - v) / d2 * SplineBlend(curBin+1,splineOrder-1,knots,v,numBins);
    } else if (d2 == 0) {
      value = (v - knots[curBin]) / d1 * SplineBlend(curBin,splineOrder-1,knots,v,numBins);
    } else {
      value = (v - knots[curBin]) / d1 * SplineBlend(curBin,splineOrder-1,knots,v,numBins) +
        (knots[curBin+splineOrder] - v) / d2 * SplineBlend(curBin+1,splineOrder-1,knots,v,numBins);
    }
  }
  if (value < 0)  {
    //  this does not happen...
    value = 0; // rounding sometimes makes this < 0, e.g. -0.000000001 
  }
  return(value);
}

// [[Rcpp::export]]
NumericVector SplineBlendBins( int splineOrder, NumericVector knots, double v, int numBins){
  NumericVector retVal(numBins);
  for(int curBin = 0; curBin < numBins; curBin++ ){
    retVal[curBin] = SplineBlend( curBin, splineOrder, knots, v, numBins);
  }
  
  return(retVal);
}

// [[Rcpp::export]]
NumericVector SplineBlendAll( NumericVector z, NumericVector knots, int splineOrder, int numBins){
  int len = z.size();
  NumericVector retVal(numBins*len);
  
  for( int i = 0; i < len; i++){
    for(int curBin = 0; curBin < numBins; curBin++ ){
      retVal[curBin+i*numBins] = SplineBlend( curBin, splineOrder, knots, z[i], numBins);
    }
  }
  
  return(retVal);
}

// // w1, w2: weights for two genes [bins x samples]
// // returns: p, the joint probability distribution [bins x bins]
// // [[Rcpp::export]]
// NumericVector hist2d(NumericMatrix w1, NumericMatrix w2){
//   
//   int i,j,k;
//   int nBins = w1.nrow();
//   int nSamples = w1.ncol();
//   
//   NumericMatrix retVal = NumericMatrix( nBins, nBins);
//   
//   for(i = 0; i < nSamples; i++){
//     for(j = 0; j < nBins; j++){
//       for(k = 0; k < nBins; k++){
//         retVal(j,k) += w1(j,i) * w2(k,i);
//       }
//     }
//   }
//   // rowMeans(sapply(1:ncol(w1),function(k){w1[,k] %*% t(w2[,k])}))
//   // p is mean of weights over all samples per 2d bin
//   // p: array [numBins x genes]
//   return( retVal / nSamples);
// }

// // [[Rcpp::export]]
// double entropy2d(NumericMatrix w1, NumericMatrix w2){
//   
//   NumericVector pMat = hist2d(w1, w2);
// 
//   // calculate entropy
//   // -sum(p * log2(p), na.rm = T)
//   double H = 0.0; //returned entropy
//   double p;
//   int psize = pMat.size();
//   for(int i = 0; i < psize; i++){
//     p = pMat[i];
//     // need to skip the 0's as it would generate NaN
//     if( p > 0 ){
//       H -= p*log2(p);
//     }
//   }
//   return(H);
// }
// 

// [[Rcpp::export]]
NumericVector hist2dC(NumericVector weights, int i1, int i2){
  // weights is a 3d array [bins x samples x genes ]
  NumericVector dim = weights.attr("dim");
  int nBins = dim[0];
  int nSamples = dim[1];

  int w1 = nBins*nSamples*i1; // add nBins*nSamples*i1 to get weights for gene i1
  int w2 = nBins*nSamples*i2; // add nBins*nSamples*i2 to get weights for gene i2
  
  int i,j,k;
  
  NumericMatrix retVal = NumericMatrix( nBins, nBins);
  
  for(i = 0; i < nSamples; i++){
    for(j = 0; j < nBins; j++){
      for(k = 0; k < nBins; k++){
        retVal(j,k) += weights[j + i*nBins + w1] * weights[k + i*nBins + w2];
      }
    }
  }
  // rowMeans(sapply(1:ncol(w1),function(k){w1[,k] %*% t(w2[,k])}))
  // p is mean of weights over all samples per 2d bin
  // p: array [numBins x genes]
  return( retVal / nSamples);
}

// [[Rcpp::export]]
double entropy2dC(NumericVector weights, int i1, int i2){
  NumericVector pMat = hist2dC(weights, i1, i2);
  
  // calculate entropy
  // -sum(p * log2(p), na.rm = T)
  double H = 0.0; //returned entropy
  double p;
  int psize = pMat.size();
  for(int i = 0; i < psize; i++){
    p = pMat[i];
    // need to skip the 0's as it would generate NaN
    if( p > 0 ){
      H -= p*log2(p);
    }
  }
  return(H);
}

    