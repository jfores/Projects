
#Get drugbank LINCS signature mapping from All_LINCS object.

get_mapping_LINCS_DB <- function(x_All_LINCS){
  x_All_LINCS <- x_All_LINCS[,c("sig_id","pert_id","pert_iname","pert_type","cell_id","pert_idose","pert_itime","Canonical_smiles","primary_key","name","Phase")]
  x_All_LINCS <- x_All_LINCS[!is.na(x_All_LINCS$primary_key),]
  x_All_LINCS <- unique(x_All_LINCS)
  return(x_All_LINCS)
}

#Given an ICD code or a set of ICD codes porvided as a strig separated by a bar | (example "331.0|692.71") it returns the Drugbank codes for its drug indications reported in MEDI-AN Hig precission set.

from_ICD_to_Indicated_DB <- function(icd_code,x_main_an_hp_drugs_all){
  x_main_an_hp_drugs_all <- x_main_an_hp_drugs_all[grepl(icd_code,x_main_an_hp_drugs_all$ICD9),]
  out_list <- list(unique(x_main_an_hp_drugs_all$primary_key))
  names(out_list) <- paste(unique(x_main_an_hp_drugs_all$ICD9),collapse = "-")
  return(out_list)
}

#Ginen a set of drugbank IDs this function returns a list of the available LINCS L1000 signatures generated using this drug as a perturbation.


