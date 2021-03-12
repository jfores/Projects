
#Get drugbank LINCS signature mapping from All_LINCS object.

get_mapping_LINCS_DB <- function(x_All_LINCS){
  x_All_LINCS <- x_All_LINCS[,c("sig_id","pert_id","pert_iname","pert_type","cell_id","pert_idose","pert_itime","Canonical_smiles","primary_key","name","Phase")]
  x_All_LINCS <- x_All_LINCS[!is.na(x_All_LINCS$primary_key),]
  x_All_LINCS <- unique(x_All_LINCS)
  return(x_All_LINCS)
}

