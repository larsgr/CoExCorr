calcCLR <- function(mi){
  miZ <- scale(mi)
  gc()
  # pmax does not support long vectors so it is applied columnwise
  miZ <- apply(miZ,2,pmax,0)^2
  gc()
  sqrt(miZ + t(miZ))
}