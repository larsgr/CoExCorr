library(biomaRt)

downloadGeneNames <- function(){
  lapply(
    c(Os="osativa_eg_gene",
      At="athaliana_eg_gene"),
    function(dataset){
      getBM(
        attributes=c("ensembl_gene_id", "external_gene_id"), 
        mart = useMart("plants_mart",
                       dataset=dataset, 
                       host="plants.ensembl.org"))
    })
}


getGeneNames <- function(geneIDs){
  spcIDs <- c(AT="At",OS="Os")[substr(geneIDs,1,2)]
  for( spc in unique(spcIDs) ){
    idx <- spcIDs==spc
    geneNames_ <- geneNames[[spc]]$external_gene_id[match(geneIDs[idx],geneNames[[spc]]$ensembl_gene_id)]
    geneIDs[idx] <- ifelse(geneNames_=="",geneIDs[idx],paste0(spc, geneNames_) )
  }
  return(geneIDs)
}