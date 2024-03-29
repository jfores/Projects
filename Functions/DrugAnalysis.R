
# Get drugbank LINCS signature mapping from All_LINCS object.

get_mapping_LINCS_DB <- function(x_All_LINCS){
  x_All_LINCS <- x_All_LINCS[,c("sig_id","pert_id","pert_iname","pert_type","cell_id","pert_idose","pert_itime","Canonical_smiles","primary_key","name","Phase")]
  x_All_LINCS <- x_All_LINCS[!is.na(x_All_LINCS$primary_key),]
  x_All_LINCS <- unique(x_All_LINCS)
  return(x_All_LINCS)
}

# Given an ICD code or a set of ICD codes porvided as a strig separated by a bar | (example "331.0|692.71") it returns the Drugbank codes for its drug indications reported in MEDI-AN Hig precission set.

from_ICD_to_Indicated_DB <- function(icd_code,x_main_an_hp_drugs_all){
  x_main_an_hp_drugs_all <- x_main_an_hp_drugs_all[grepl(icd_code,x_main_an_hp_drugs_all$ICD9),]
  out_list <- list(na.omit(unique(x_main_an_hp_drugs_all$primary_key)))
  names(out_list) <- paste(unique(x_main_an_hp_drugs_all$ICD9),collapse = "-")
  return(out_list)
}

# Given a set of drugbank IDs this function returns a list of the available LINCS L1000 signatures generated using this drug as a perturbation.

get_LINCS_profiles_from_DB_IDs <- function(x_all_lincs,vec_ids){
  x_all_lincs <- x_all_lincs[x_all_lincs$primary_key %in% vec_ids,]
  x_all_lincs <- unique(x_all_lincs[,c("primary_key",c("sig_id"))])
  x_all_lincs <- split(x_all_lincs,as.factor(x_all_lincs$primary_key))
  return(x_all_lincs)
}

########################################################
#Functions to generate stouffers consensus signatures.##
########################################################

#Given a parturbation retrieve the level 5 expression matrix containing the data for all signatures generated employing it in LINCS Pahses I and II.

retrieve_signatures <- function(perturbation,phase_I_path,phase_II_path,all_lincs_to_stouff_data){
  require(cmapR)
  all_lincs_to_stouff_data <- all_lincs_to_stouff_data[perturbation][[1]]
  #View(all_lincs_to_stouff_data)
  sig_phases <- unique(all_lincs_to_stouff_data$Phase)
  print(sig_phases)
  if(length(sig_phases) == 2){
    sigs_to_pert_pI <- all_lincs_to_stouff_data[all_lincs_to_stouff_data$Phase == "Phase_I",]
    pI_mat <- parse_gctx(phase_I_path,cid = sigs_to_pert_pI$sig_id)
    pI_mat <- pI_mat@mat
    sigs_to_pert_pII <- all_lincs_to_stouff_data[all_lincs_to_stouff_data$Phase == "Phase_II",]
    pII_mat <- parse_gctx(phase_II_path,cid = sigs_to_pert_pII$sig_id)
    pII_mat <- pII_mat@mat
    full_mat <- cbind(pI_mat,pII_mat)
  }else if(sig_phases == "Phase_I"){
    sigs_to_pert_pI <- all_lincs_to_stouff_data[all_lincs_to_stouff_data$Phase == "Phase_I",]
    pI_mat <- parse_gctx(phase_I_path,cid = sigs_to_pert_pI$sig_id)
    full_mat <- pI_mat@mat
  }else if(sig_phases == "Phase_II"){
    sigs_to_pert_pII <- all_lincs_to_stouff_data[all_lincs_to_stouff_data$Phase == "Phase_II",]
    pII_mat <- parse_gctx(phase_II_path,cid = sigs_to_pert_pII$sig_id)
    full_mat <- pII_mat@mat
  }
  return(full_mat)
}

#Compute Stouffer's consensus signature from a particular matrix. We have to pass the data always in df format.

compute_Stouffers_Consensus <- function(exp_mat){
  if(ncol(exp_mat) == 1){
    print("One column...")
    return(exp_mat[,1])
  }else if(ncol(exp_mat) == 2){
    print("Two columns...")
    return(rowSums(exp_mat)/2)
  }else{
    print("More than 2 columns...")
    #Computing weights for each signature.
    cor_mat <- cor(exp_mat,method = "spearman")
    diag(cor_mat) <- NA
    corr_to_all_other <- rowSums(cor_mat,na.rm = T)/(ncol(cor_mat)-1)
    weights <- corr_to_all_other/sum(corr_to_all_other)
    #Compute consensus.
    numer_sig <- rowSums(t(t(exp_mat) * weights))
    denom_sig <- sqrt(sum(weights^2))
    consensus_sig <- numer_sig/denom_sig
    return(consensus_sig)
  }
}

##########################
##Computing correlations##
##########################

compute_correlations <- function(x,y){
  intersected_rows <- intersect(rownames(x),rownames(y))
  x <- x[intersected_rows,]
  y <- y[intersected_rows,]
  cor_values <- apply(y,2,function(x,y) cor.test(x,y,method = "spearman")$estimate,x$zval)
  p_values <- apply(y,2,function(x,y) cor.test(x,y,method = "spearman")$p.value,x$zval)
  return(list(cor_values,p_values))
}


#########################
##Enrichment Functions###
#########################

#Enrich multiple signatures.

enrich_clusterprof_mult <- function(x,format_rownames = T,num_chars = 30){
  require(clusterProfiler)
  c2_canonical <- read.gmt("/home/jaume/Projects/mol_sig/c2.cp.v7.2.entrez.gmt")
  if(format_rownames){
    c2_canonical$term <- gsub("_", " ",gsub("^ST","(S)",gsub("^WP","(W)",gsub("^PID","(P)",gsub("^NABA","(N)",gsub("^BIOCARTA","(B)",gsub("^KEGG","(K)",gsub("^REACTOME","(R)",c2_canonical$term))))))))
    c2_canonical$term[nchar(as.character(c2_canonical$term)) > num_chars] <- paste(substr(c2_canonical$term[nchar(as.character(c2_canonical$term)) > num_chars],1,30),"...")
    
  }
  list_out <- list()
  print(ncol(x))
  for(i in 1:ncol(x)){
    print(i)
    temp_prof <- x[,i]
    names(temp_prof) <- rownames(x)
    temp_prof <- temp_prof[order(temp_prof,decreasing = T)]
    list_out[[i]] <- GSEA(temp_prof,TERM2GENE = c2_canonical)
  }
  names(list_out) <- colnames(x)
  return(list_out)
}

#Filter prior to ridgeplot.

filter_top_enrich_to_plot <- function(x,n_up = 20,n_down = 20){
  top_up <- which(x@result$NES > 0)[1:20]
  top_down <- which(x@result$NES < 0)[1:20]
  all  <- c(top_up,top_down)
  x@result <- x@result[all,]
  return(x)
}

