library(dplyr)

geneNames <- 
  dir("indata/geneInfo", pattern="[A-Z][a-z]_geneInfo.txt",full.names = T) %>% 
  setNames( sub(".*([A-Z][a-z])_geneInfo.txt","\\1",.)) %>% 
  map( read_tsv, col_types=cols() ) %>%
  bind_rows( .id = "spc")

getGeneNames <- function(geneIDs){
  data_frame(ensembl_gene_id = geneIDs) %>% 
    left_join(geneNames, by = "ensembl_gene_id") %>% 
    transmute( geneName = ifelse(is.na(external_gene_id),
                                 ensembl_gene_id,
                                 paste0(spc, external_gene_id) ) ) %>% 
    .$geneName
}