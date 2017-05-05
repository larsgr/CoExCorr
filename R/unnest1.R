

# flatten one level of a nested list inside a data.frame
unnest1 <- function(x, .id){
  listColIdx <- which(map_lgl(x,is.list))
  if( length(listColIdx) != 1) stop("Must be exactly one list column")
  
  # extract the list column
  l <- x[[listColIdx]]
  x <- x[-listColIdx]
  
  # replace it with a list of names
  x[[.id]] <- lapply(l, function(y){
    if(is.null(names(y))){
      return(as.character(seq_along(y)))
    } else {
      return(names(y))
    }
  })
  # unnest names
  x <- unnest(x)
  # add unnested(flattened) list
  x[[names(listColIdx)]] <- purrr::flatten(l)
  
  return(x)
}