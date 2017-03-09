# library(readxl)

spc2PODCspcNames <- 
  c(Os = "oryza_sativa",
    At = "arabidopsis_thaliana",
    Sl = "solanum_lycopersicum",
    Sb = "sorghum_bicolor",
    Vv = "vitis_vinifera",
    Mt = "medicago_truncatula",
    St = "solanum_tuberosum",
    Gm = "glycine_max",
    Nt = "nicotiana_tabacum",
    Zm = "zea_mays")

getPathPODCmeta <- function(spc){
  file.path("indata",paste0(spc2PODCspcNames[spc],"_sample.xlsx"))
}

loadPODCmeta <- function(spc){
  readxl::read_excel( getPathPODCmeta(spc),col_types = rep("text",15))
}

getPathPODCexpMat <- function(spc){
  file.path("indata",paste0(spc2PODCspcNames[spc],"_fpkm.tsv.gz"))
}

getPathSampleMeta <- function(spc){
  file.path("indata/sampleMeta",paste0(spc,".RDS"))
}

loadSampleMeta <- function(spc){
  readRDS(getPathSampleMeta(spc))
}
