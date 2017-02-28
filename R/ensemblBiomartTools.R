library(biomaRt)

#
# The names of the relevant organisms
#
spc2ensemblSpcID <- c(
  Os = "osativa",
  At = "athaliana",
  Gm = "gmax",
  Sl = "slycopersicum",
  # Sb = "sbicolor",
  # Vv = "vvinifera",
  # Mt = "mtruncatula",
  # St = "stuberosum",
  # Nt = "ntabacum",
  Zm = "zmays")

downloadOrthosEnsembl <- function(spc1, 
                                  spc2, 
                                  biomart = "plants_mart", 
                                  host = "plants.ensembl.org"){
  # get spc1 gene id and corresponding spc2 orthologs. Filter: only genes with spc2 homolog
  orthos <- getBM(
    attributes = c("ensembl_gene_id", 
                   paste0( spc2ensemblSpcID[spc2],
                           c( "_eg_homolog_ensembl_gene",
                              "_eg_homolog_orthology_type",
                              "_eg_homolog_subtype",
                              "_eg_homolog_perc_id",
                              "_eg_homolog_perc_id_r1" ) ) ), 
    filters = paste0("with_",spc2ensemblSpcID[spc2],"_eg_homolog"), values = T, 
    mart = useMart(
      biomart = biomart,
      dataset = paste0(spc2ensemblSpcID[spc1],"_eg_gene"),
      host = host))
  
  # change the column names
  colnames(orthos) <- c(spc1,spc2,"otype","MRCA","PID","PIDr")  
  # convert orthology type
  oTypes <- c(ortholog_many2many = "N:N", ortholog_one2many = "1:N", ortholog_one2one = "1:1")
  orthos$otype <- oTypes[orthos$otype]
  
  return(orthos)
}


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