#' Load triangular matrix
#' 
#' Loads a text file containing space separated values of the lower half of a 
#' matrix with one row per line (i.e. each line has one more value than the last).
#' The diagonal is not included so it will be 0. The upper half of the matrix 
#' will be filled with the transposed of the lower half to produce a symmetrical
#' matrix.
#'
#' @param geneIds Gene identifiers to apply to the rownames and colnames
#' @param filename Path to the triangular matrix file
#'
#' @return A symmetrical matrix
#'
loadTriMatrix <- function(geneIds,filename){
  n <- length(geneIds)
  m <- matrix(0, ncol = n, nrow = n, dimnames = list(geneIds, geneIds) )
  inCon <- file(filename, "rt")
  for(i in 1:(n-1)){
    # read one line at the time
    x <- scan(inCon,nlines = 1,quiet = T)
    # write row-wise
    m[i+1, 1:i] <- x
    # write column-wise
    m[1:i, i+1] <- x
  }
  close(inCon)
  return(m)
}