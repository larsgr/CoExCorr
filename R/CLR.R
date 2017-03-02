calcCLR <- function(mi){
  miZ <- pmax(scale(mi),0)^2
  gc()
  sqrt(miZ + t(miZ))
}