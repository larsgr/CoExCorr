calcCLR <- function(mi){
  miZ <- pmax(scale(mi),0)^2
  sqrt(miZ + t(miZ))
}