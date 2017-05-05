# set names according to dimnames of CCS
setRnkNames <- function(rnks, CCS){
  if(!is.null(dimnames(CCS))){
    rnks <- lapply(1:2, function(x){
      setnames(rnks[[1]], dimnames(CCS)[[x]])
    })
  }
  return(rnks)
}


# calcRanks
calcRanks <- function(CCS){
  cat("Calculating ranks...\n")
  nGenes <- ncol(CCS)
  lapply(1:2, function(x){
    diag(apply(CCS,x, rank))/nGenes
  })
}

# calc rank of sum CCS in both directions
calcRanks2 <- function(CCS){
  cat("Calculating ranks...\n")
  nGenes <- ncol(CCS)
  diag(apply(CCS+t(CCS),1, rank))/nGenes
}

# Calculate ranks for specified ortholog pairs
# Ranks are relative to spc1 which are on the rows of CCS
quickRanks <- function(CCS,spc1genes,spc2genes){
  # for each pair of orthologs
  # get the number of genes with lower or equal CCS
  mapply(spc1genes,spc2genes, FUN=function(spc1gene,spc2gene){
    mean(CCS[spc1gene, ]<=CCS[spc1gene, spc2gene],na.rm = T)
  })
}

# Calculate ranks for specified ortholog pairs (Transposed version)
# Ranks are relative to spc1 which are on the columns of CCS
quickRanksT <- function(CCS,spc1genes,spc2genes){
  # for each pair of orthologs
  # get the number of genes with lower or equal CCS
  mapply(spc1genes,spc2genes, FUN=function(spc1gene,spc2gene){
    mean(CCS[ ,spc1gene]<=CCS[spc2gene, spc1gene],na.rm = T)
  })
}

# Automatically detect if quickRanks or quickRanksT should be used
# based on geneIDs
quickRanksA <- function(CCS,spc1genes,spc2genes){
  if( spc1genes[1] %in% rownames(CCS) ){
    return( quickRanks(CCS,spc1genes,spc2genes) )
  } else{
    return( quickRanksT(CCS,spc1genes,spc2genes) )
  }
}


# Calculates and draws a ROC like curve from a rank vector
drawRanksCurve <- function(rnk,by=1/length(rnk),...){
  h <- hist(rnk,breaks=seq(0,1,by=by),plot=F)
  xy <- xy.coords(x=1-h$breaks,y=c(0,cumsum(rev(h$counts))/length(rnk)))
  plot.xy(xy,type="S",...)
}

# Creates a blank plot for drawing rank curves
newRanksPlot <- function(main='Ortholog ranks',
                         ylab="Proportion of genes with rank above threshold",
                         xlab="Rank threshold",...){
  plot(x=NULL, xlim=c(1,0),ylim=c(0,1),
       ylab=ylab,xlab=xlab,
       main=main,...)
}
