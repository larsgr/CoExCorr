# parallel correlation

library(parallel)

# cross correlate two matrices in parallel
mc_cor2 <- function(M1, M2, cores){
  # split M2 by columns
  idx <- split(1:ncol(M2),cut(1:ncol(M2),cores,labels = FALSE))
  names(idx) <- NULL
  do.call(cbind, mclapply(idx,function(i){
    cor(M1,M2[ ,i])
  },mc.cores = cores))
}
