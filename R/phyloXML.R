library(XML)
library(phylobase)

####
#
# Read phyloXML and store as a flat list of clades (leaves + nodes)
#

loadPhyloXMLasList <- function( xmlFile ){
  
  # recursive function to flatten the phyloXML hierarchy
  phyloFlat <- function(parentID=NA, cladeID = "r", cc, leafEnv, nodeEnv){
    
    # Separate subclades from data for this clade
    isSubClade <- names(cc) == "clade"
    subClades <- cc[isSubClade]
    thisCladeData <- cc[!isSubClade]
    
    # add parentID to data
    thisCladeData$parentID <- parentID
    
    # save data about this clade  
    if(length(subClades) > 0){
      # has sub-clades, i.e. it is internal node, so save to nodeEnv
      assign(cladeID, thisCladeData, envir=nodeEnv)
    } else {
      # is leaf node so save to leafEnv
      assign(cladeID, thisCladeData, envir=leafEnv)
    }
    
    # recursively iterate through sub-clades
    for(j in seq_along(subClades)){
      subCladeID <- paste(cladeID,j,sep = ".") # suffix subclades with a number
      Recall(parentID=cladeID, cladeID = subCladeID, cc = subClades[[j]], leafEnv, nodeEnv) # recursive call
    }
  }
  
  # load xml file and convert to list
  xl <- xmlToList(xmlParse(xmlFile))

  # Create two environments that will be used to put the clade data
  leafEnv <- new.env()
  nodeEnv <- new.env()
  
  # Flatten
  phyloFlat(cc=xl$phylogeny$clade,leafEnv=leafEnv,nodeEnv=nodeEnv)

  # convert to lists
  nodes <- as.list(nodeEnv)
  leaves <- as.list(leafEnv)
  # Sort to preserve order
  return(list(nodes = nodes[order(names(nodes))], leaves = leaves[order(names(leaves))]))
}

#####
#
# Get phylogeny without metadata from the flattened phyloXML structure
#

listToTree <- function(leaves,nodes){
  # generate edge matrix
  edges <- sapply(c(leaves,nodes),function(x){x$parentID})
  m <- cbind(ancestor=match(edges,names(edges)),descendant=seq_along(edges))
  m[is.na(m)]<- 0
  
  # get edge lengths
  edge.length <- sapply(c(leaves,nodes),function(x){as.numeric(x$.attrs["branch_length"])})
  
  phylo4(m,
         edge.length = edge.length,
         tip.label = unname(sapply(leaves,function(x){x$name})))
}



#####
#
# Functions to extract data elements from the phyloXML clade structure
#

getBranchLength <- function(clade){
  as.numeric(clade$.attrs["branch_length"])
}

getName <- function(clade){
  clade$name
}

getTaxonID <- function(clade){
  clade$taxonomy$id
}

getDuplicationsEvent <- function(clade){
  i <- which( names(clade) == "events")
  if( length(i) > 0 ){
    if(length(i) > 1){
      stop("multiple events!")
    }
    if( clade[[i]]$type != "speciation_or_duplication"){
      stop("unknown event!")
    }
    if(clade[[i]]$duplications != "1"){
      stop("duplications != 1")
    }
    return(TRUE)
  } else {
    return(FALSE)
  }
}

getBootstrap <- function(clade){
  retVal <- NA_integer_
  i <- which(names(clade) == "confidence")
  if( length(i) > 0){
    j <- which(sapply( clade[i], function(x){x$.attrs["type"] == "bootstrap"}))
    if( length(j) > 0){
      retVal <- as.integer(clade[i][[j]]$text)
    }
  }
  return(retVal)
}

getDuplicationConfidence <- function(clade){
  retVal <- NA_real_
  i <- which(names(clade) == "confidence")
  if( length(i) > 0){
    j <- which(sapply( clade[i], function(x){x$.attrs["type"] == "duplication_confidence_score"}))
    if( length(j) > 0){
      retVal <- as.numeric(clade[i][[j]]$text)
    }
  }
  return(retVal)
}

###
#
# Extract all interresting data from tree
#
extractTreeData <- function(tl){
  with(tl, {
    # generate edge matrix
    edges <- sapply(c(leaves,nodes),function(x){x$parentID})
    m <- cbind(ancestor=match(edges,names(edges)),descendant=seq_along(edges))
    m[is.na(m)]<- 0L
    
    # get tip.label
    tip.label <- unname(sapply(leaves, getName))
    
    # get common data (edge lengths, taxon id)
    edge.length <- unname(sapply(c(leaves,nodes),getBranchLength))
    
    
    
    # get metadata
    tip.data <- data.frame( row.names = seq_along(leaves)+length(leaves),
                            taxonID = sapply(leaves,getTaxonID))
    node.data <- data.frame( row.names = seq_along(nodes)+length(leaves),
                             taxonID = sapply(nodes,getTaxonID),
                             bootstrap = sapply(nodes, getBootstrap),
                             isDuplication = sapply(nodes, getDuplicationsEvent),
                             duplicationConfidence = sapply(nodes, getDuplicationConfidence),
                             stringsAsFactors = F)
    return( list( m = m, edge.length = edge.length, tip.label = tip.label, tip.data = tip.data, node.data = node.data))
  })
}


treeDataToPhylo4 <- function(treeData){
  with(treeData, {
    phylo4(m, edge.length, tip.label)
  })
}

treeDataToPhylo4d <- function(treeData){
  with(treeData, {
    addData( phylo4(m, edge.length, tip.label),
             tip.data=tip.data,
             node.data=node.data)
  })
}
