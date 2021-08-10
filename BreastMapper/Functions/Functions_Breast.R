#Functions For Breast TDA prject.

##########################
##Working with arrays#####
##########################

#Lectura de datos phenotípicos desde fichero soft

library(plyr)

#Auxiliar function to work with read_pheno_from_soft

get_pData <- function(x,platform){
  if(Meta(x)$platform_id == platform){
    return(as.data.frame(Meta(x)))
  }
}


#Getting contained platforms from soft readed file.

get_Platforms <- function(x){
  return(unlist(lapply(GPLList(x),function(x) Meta(x)$geo_accession)))
}


read_pheno_from_soft <- function(x,platform){
  list_GSMs <- GSMList(x)
  pData <- lapply(list_GSMs,get_pData,platform)
  pData <- do.call("rbind.fill",pData)
  rownames(pData) <- pData$geo_accession
  return(pData)
}

#Divide columns with multiple characteristics information.

divide_column <- function(x,column,sep=","){
  pData_B <- do.call("rbind",lapply(lapply((lapply(strsplit(x[,column],sep),trimws)),as.data.frame),t))
  colnames(pData_B) <- paste(rep("pTemp_",ncol(pData_B)),1:ncol(pData_B),sep="")
  return(cbind(x,pData_B))
}

#Get scan date for affymetrix data

get_scan_date <- function(x_e,x_p){
  pCh_Scan_Date <- pData(protocolData(x_e))
  pCh_Scan_Date <- pCh_Scan_Date[,1]
  pCh_Scan_Date <-gsub(" .*","",pCh_Scan_Date)
  pCh_Scan_Date <- pCh_Scan_Date[match(col_names,rownames(x_p))]
  return(pCh_Scan_Date)
}

#Merge raw and process phenoData

merge_raw_and_processed <- function(x_p_r,x_p_p){
  x_p_r <- x_p_r[rownames(x_p_p),]
  x_p_d <- cbind(x_p_p,x_p_r)
  return(x_p_d)
}


#Create a numeric vector from a character vector.

char_to_num <- function(vector){
  vector <- as.character(vector)
  unique_items <- unique(vector)
  unique_items <- sort(unique_items)
  for(i in 1:length(unique_items)){
    vector[vector == unique_items[i]] = as.character(i)
  }
  vector = as.numeric(vector)
  return(vector)
}

#Transforming values of PhenoData to numeric and factor.

pCh_to_Num_and_Fac <- function(pData){
  library(RColorBrewer)
  df_num <- data.frame(matrix(nrow = nrow(pData),ncol = ncol(pData)))
  rownames(df_num) <- rownames(pData)
  df_fac <- data.frame(matrix(nrow = nrow(pData),ncol = ncol(pData)))
  rownames(df_fac) <- rownames(pData)
  df_color <- data.frame(matrix(nrow = nrow(pData),ncol = ncol(pData)))
  rownames(df_color) <- rownames(pData)
  for( i in 1:ncol(pData)){
    df_num[,i] <- char_to_num(pData[,i])
    df_fac[,i] <- factor(pData[,i])
    if(length(levels(factor(pData[,i]))) <=8){
      print(colnames(pData)[i])
      niveles <- length(levels(factor(pData[,i])))
      print(niveles)
      brewer_pal <- brewer.pal(length(levels(factor(pData[,i]))),"Dark2")
      print(brewer_pal)
      print(brewer_pal[1:niveles])
      df_temp <- data.frame(levels(factor(pData[,i])),brewer.pal(length(levels(factor(pData[,i]))),"Dark2")[1:niveles])
      colnames(df_temp) <- c("var","color")
      print(df_temp)
      #print(head(df_temp))
      #print(dim(df_temp))
      df_temp_2 <- data.frame(pData[,i,drop = FALSE],rownames(pData))
      colnames(df_temp_2) <- c("var","rownames")
      #print(head(df_temp_2))
      #print(dim(df_temp_2))
      df_temp_merged <- merge(df_temp_2,df_temp,by = "var",all.x = TRUE)
      rownames(df_temp_merged) <- df_temp_merged[,"rownames"]
      df_temp_merged <- df_temp_merged[rownames(pData),]
      print(head(df_temp_merged))
      #print(dim(df_temp_merged))
      #print(head(df_temp_merged))
      df_color[,i] <- df_temp_merged[,"color"]
    }else{
      df_color[,i] <- rep("NSFTP",nrow(pData))
    }
  }
  colnames(df_num) <- gsub("pCh_","pNm_",colnames(pData))
  colnames(df_fac) <- gsub("pCh_","pFac_",colnames(pData))
  colnames(df_color) <- gsub("pCh_","pColor_",colnames(pData))
  out <- cbind(pData,df_num,df_fac,df_color)
  return(out)
}

#Creating structures

# Create folder structure

Create_Structure <- function(root_dir,list_of_GSE){
  for(i in 1:length(list_of_GSE)){
    if(!dir.exists(file.path(root_dir,list_of_GSE[i]))){
      dir.create(file.path(root_dir,list_of_GSE[i]))
      dir.create(file.path(root_dir,list_of_GSE[i],"Data"))
      dir.create(file.path(root_dir,list_of_GSE[i],"Scripts"))
      dir.create(file.path(root_dir,list_of_GSE[i],"PaperNotesTables"))
      dir.create(file.path(root_dir,list_of_GSE[i],"Results"))
    }else{
      print(paste("Directory",list_of_GSE[i],"already exists"))
    }
  }
}


#Auxiliar function for get_pData_table

unfold_charact <- function(GSM){
  GSM_meta <- Meta(GSM)
  characteristics <- GSM_meta$characteristics_ch1
  for(i in 1:length(characteristics)){
    if(length(strsplit(characteristics[i],":")[[1]]) > 1){
      GSM_meta[[paste("pCh_",gsub(" ","_",trimws(strsplit(characteristics[i],":")[[1]][1])),sep="")]] <- gsub(" ","_",trimws(strsplit(characteristics[i],":")[[1]][2]))
    }
  }
  GSM_meta[["characteristics_ch1"]] <- NULL
  for(i in 1:length(GSM_meta)){
    if(length(GSM_meta[[i]])>1){
      GSM_meta[[i]] <- paste(GSM_meta[[i]],collapse="/")
    }
  }
  return(GSM_meta)
}

#getting pheno data from soft

get_pData_table <- function(GSE,GPL){
  counter_temp <- 0
  for(i in 1:length(GSMList(GSE))){
    if(Meta(GSMList(GSE)[[i]])$platform_id == GPL){
      if(counter_temp == 0){
        df_sal <- data.frame(unfold_charact(GSMList(GSE)[[i]]))
      }else{
        df_sal <- rbind.fill(df_sal,data.frame(unfold_charact(GSMList(GSE)[[i]])))
      }
      counter_temp = counter_temp + 1
    }
  }
  return(df_sal)
}


#Expression data from soft.

get_eData <- function(GSE,GPL){
  counter_temp <- 0
  for(i in 1:length(GSMList(GSE))){
    if(Meta(GSMList(GSE)[[i]])$platform_id == GPL){
      #print(counter_temp)
      if(counter_temp == 0){
        e_sal <- Table(GSMList(GSE)[[i]])
        e_sal <- e_sal[-ncol(e_sal)]
        colnames(e_sal)[ncol(e_sal)] <- Meta(GSMList(GSE)[[i]])$geo_accession
        #print(head(e_sal))
      }else{
        e_sal_temp <- Table(GSMList(GSE)[[i]])
        e_sal_temp <- e_sal_temp[-ncol(e_sal_temp)]
        e_sal <- merge(e_sal,e_sal_temp,by="ID_REF")
        colnames(e_sal)[ncol(e_sal)] <- Meta(GSMList(GSE)[[i]])$geo_accession
        #print(head(e_sal))
      }
      counter_temp = counter_temp + 1
    }
  }
  rownames(e_sal) <- e_sal[,1]
  e_sal <- e_sal[,-1]
  return(e_sal)
}


###########################
###Array quality metrics###
###########################



compute_max_connectivity <- function(connectivity){
  return(max(connectivity))
}

compute_scaled_connectivity <- function(connectivity,max_connectivity){
  scaled_connectivity <- connectivity/max_connectivity
  scaled_connectivity[is.na(scaled_connectivity)] <- 0
  return(scaled_connectivity)
}

compute_standardized_connectivity <- function(scaled_connectivity){
  standardized_connectivity <- (scaled_connectivity - mean(scaled_connectivity))/(sqrt(var(scaled_connectivity)))
  return(standardized_connectivity)
}

compute_standardized_connectivity_fadj <- function(adj){
  scaled_connectivity <- compute_scaled_connectivity(compute_connectivity(adj),compute_max_connectivity(compute_connectivity(adj)))
  standardized_connectivity <- (scaled_connectivity - mean(scaled_connectivity))/(sqrt(var(scaled_connectivity)))
  return(standardized_connectivity)
}


#Compute triplets. First computes the number of common neighbours between 2 nodes. Then if this twe neighbours are neighbours between them.
#Matrix multiplication of adjacency matrix tells us how many common neighbours have two nodes. Then by performing element-wise multiplication between
#the previous matrix and the original adjacency matrix we determine if the pair of neighbours sharing a partircular number of common neighbours are neighbours 
#between them and can be cosidered connected triplets.

compute_n_common_neigbours <- function(adjacency_matrix){
  common_neigbours <- (adjacency_matrix %*% adjacency_matrix) * adjacency_matrix
  return(common_neigbours)
}

#Computes the number of neigburs of each node that are neigbours between them. For each node the number of neighbours that are neighbours between them will be 
#determined by the rowsums of the matrix generated by the previous funcion divided by two. 

compute_n_i <- function(adjacency_matrix){
  n_i <- rowSums(compute_n_common_neigbours(adjacency_matrix))/2
  n_i[is.na(n_i)] <- 0
  #print(n_i)
  return(n_i)
}

#Calcular el máximo de nodos vecinos comunes para cada nodo en base al número de vecinos. Ojo que aquí lo que le tenemos que pasar es el número de vecinos no es
#el número de vecinos conectados.
compute_max_con_com_neig <- function(adj_mat){
  vec <- rowSums(adj_mat)
  max_conn <- (vec*(vec-1))/2
  max_conn[is.na(max_conn)] <- 0
  #print(max_conn)
  return(max_conn)
}

compute_clustering_coeffincient_unweighed <- function(adj_mat){
  cc <- compute_n_i(adj_mat)/compute_max_con_com_neig(adj_mat)
  cc[is.na(cc)] <- 0
  return(cc)
}


#En el caso de las redes de coexpresión positivas generalizar n_i es facil. Aplicaremos la misma fórmula. La dificualtad estriba en generalizar Pii. ni tiene que
#ser menor que pii y ni tiene que ser igual a pii para una red completamente conectada.

compute_max_con_com_neig_weigh <- function(adj_mat){
  max_conn <- rowSums(adj_mat)^2 - rowSums(adj_mat^2)
  max_conn <- max_conn/2
  return(max_conn)
}

compute_clustering_coeffincient_weighed <- function(adj_mat){
  cc <- compute_n_i(adj_mat)/compute_max_con_com_neig_weigh(adj_mat)
  cc[is.na(cc)] <- 0
  return(cc)
}

#Compute standardized clustering coefficient.

compute_stadardized_clust_con_weigh <- function(adj_mat){
  cc <- compute_clustering_coeffincient_weighed(adj_mat)
  z_cc <- (cc - mean(cc))/sqrt(var(cc))
  return(z_cc)
}


plot_cor_cc_con <- function(adj_mat,col,bool,main_title){
  st_con <- compute_standardized_connectivity_fadj(adj_mat)
  #print(head(st_con))
  st_cc <- compute_stadardized_clust_con_weigh(adj_mat)
  #print(head(st_cc))
  print(bool)
  cor_case <- round(cor(st_con[bool],st_cc[bool],method="spearman"),2)
  print(cor_case)
  cor_control <- round(cor(st_con[!bool],st_cc[!bool],method="spearman"),2)
  print(cor_control)
  plot(st_con,st_cc,col = col,ylab="Z.C",xlab="Z.K",main=paste(main_title,": ","Case: ",cor_case," ","Control: ",cor_control),xlim=c(-4,4),ylim=c(-4,4))
  abline(lm(st_cc[bool]~st_con[bool]),col="red")
  abline(lm(st_cc[!bool]~st_con[!bool]),col="black")
}

#Given a pallete and a vector transforms it to colors.

transform_to_color <- function(vec,pal_name){
  unique_vals <- unique(vec)
  vec_cols <- c()
  n <- length(unique(vec))
  library("RColorBrewer")
  brewer <- brewer.pal(n, pal_name)
  for(i in 1:length(vec)){
    vec_cols <- c(vec_cols,brewer[which(unique_vals %in% vec[i])])
  }
  return(vec_cols)
}

transform_to_color_cust <- function(vec,your_pal_name){
  unique_vals <- unique(vec)
  vec_cols <- c()
  n <- length(unique(vec))
  for(i in 1:length(vec)){
    vec_cols <- c(vec_cols,your_pal_name[which(unique_vals %in% vec[i])])
  }
  return(vec_cols)
}


plot_standardize_connectivity <- function(adj_mat,col){
  st_con <- compute_standardized_connectivity_fadj(adj_mat)
  plot(st_con,col=col)
}

sample_network_analysis <- function(eset,group,case_cont){
  par(mfrow=c(2,2))
  for(i in 1:length(unique(group))){
    print(i)
    group_bool <- group %in% unique(group)[i]
    adj_temp <- compute_adjacency(eset[,group_bool])
    case_cont_temp <- case_cont[group_bool]
    case_cont_bool <- ifelse(case_cont_temp == "T",TRUE,FALSE)
    col = ifelse(case_cont_temp == "T","red","black")
    #print(col)
    plot_cor_cc_con(adj_temp,col,case_cont_bool,main_title = unique(group)[i])
  }
}

#Z.K Connectivities
#Z.C Clustering

#################
##NEW FUNCTIONS##
#################

#Plots a dendrogram with 

plot_Dendro <- function(eset,label){
  library("ClassDiscovery")
  dist <- 1 - compute_adjacency(eset)
  hc <- hclust(as.dist(dist),method = "average") #hc$labels mantiene el orden inicial de las etiquetas
  plotColoredClusters(hc,labs = hc$labels ,col = pData(eset)[,paste("pColor_",label,sep="")],cex =  0.5)
  legend <- unique(pData(eset)[,c(paste("pColor_",label,sep=""),paste("pCh_",label,sep=""))])
  print(legend)
  legend("topright",legend = legend[,2],col=legend[,1],lty=1)
  return(hc)
}

#Computing adjacency for each pair of samples. Technically aij is a signed weighted adjancency matrix. β = 2 results in an adjacency measure that is close to the correlation when the correlation is large (e.g. larger than 0.6, which is often the case among samples in microarray data).

compute_adjacency <- function(eset){
  edata <- exprs(eset)
  power <- 2
  adjacency <- (((cor(exprs(eset),method = "pearson")+1)/2)^power)
  diag(adjacency) <- 0
  return(adjacency)
}

compute_connectivity <- function(adjacency){
  connectivity <- rowSums(adjacency)
  return(connectivity)
}

#K is connectivity

compute_K <- function(eset){
  #WGCNA computation
  #library("WGCNA")
  #n_concepts <- networkConcepts(exprs(eset),networkType = "signed",power = 2) #La adyacecia de las muestras se tiene que
  #Calcular con los parametros anteriores.
  #k_wgcna <- n_concepts$Connectivity 
  adj <- compute_adjacency(eset)
  k_my <- compute_connectivity(adj)
  #print(head(k_wgcna))
  print(head(k_my))
  #print(table(k_wgcna == k_my))
  return(k_my)
}

compute_Ki <- function(eset){
  k <- compute_K(eset)
  k_max <- max(k)
  k_i <- k/k_max
  return(k_i)
}

compute_ZK <- function(eset){
  k_i <- compute_Ki(eset)
  Z_K <- (k_i - mean(k_i))/(sqrt(var(k_i)))
  return(Z_K)
}

plot_ZK <- function(eset,order = "",label,legend_on = TRUE){
  if(length(order) > 1){
    eset <- eset[,order]
  }
  ZK <- compute_ZK(eset)
  plot(ZK, col = pData(eset)[,paste("pColor_",label,sep="")])
  text(1:201,ZK,labels=colnames(eset), cex= 0.4, pos=3)
  abline(h = -2,col = "red")
  legend <- unique(pData(eset)[,c(paste("pColor_",label,sep=""),paste("pCh_",label,sep=""))])
  #print(legend)
  if(legend_on){
    legend("topright",legend = legend[,2],col=legend[,1],lty=1,cex = 0.5)
  }	
}

unravel_vec <- function(vector){
  vec_unique <- unique(vector)
  df <- data.frame(matrix(NA,ncol=length(vec_unique),nrow=length(vector)))
  for(i in 1:length(vec_unique)){
    bool_temp <- vector == vec_unique[i]
    df[,i] <- bool_temp
  }
  colnames(df) <- vec_unique
  return(df)
}

#Evaluar niveles de expression quantiles ver si los el maximo y la mediana son mayoures a un determinado nivel.

transfor_to_log <- function(eset){
  print(quantile(exprs(eset)))
  if(max(exprs(eset)) > 22){
    if(min(exprs(eset))<0){
      print("There are values under 0...")
      print(min(exprs(eset)))
      exprs(eset)[exprs(eset) < 1] <- 1
      print("Min value after offset adition...")
      print(min(exprs(eset)))
      
    }
    exprs(eset) <- log2(exprs(eset))
    print(min(exprs(eset)))
  }
  return(eset)
}


ZK_by_group <- function(eset,group,label,directory,pre_post = c("Pre_Norm")){
  
  eset <- transfor_to_log(eset)
  
  my_color_pallete <- c("yellow4","wheat4","violetred4","turquoise4","tomato4","thistle4","tan4","steelblue4","springgreen4","snow4","slategrey","steelblue","red2","purple3","plum3","palegreen3","darkorange1","black","deeppink","lightgoldenrod")
  print(unique(pData(eset)[,paste("pColor_",label,sep="")]))
  print(length(unique(pData(eset)[,c(paste("pColor_",label,sep=""),paste("pCh_",label,sep=""))])))
  if(length(unique(pData(eset)[,paste("pColor_",label,sep="")])) == 1){
    print("Yes")
    pData(eset)[,paste("pColor_",label,sep="")] <- transform_to_color_cust(pData(eset)[,paste("pColor_",label,sep="")],my_color_pallete)
  }
  bool_df <- unravel_vec(group)
  list_removal_by_group <- list()
  for(i in 1:ncol(bool_df)){
    outdir <- file.path(directory,"Results",pre_post)
    outdir_2 <- file.path(outdir,colnames(bool_df)[i])
    if(!dir.exists(outdir_2)){
      dir.create(outdir_2)
    }
    pdf(file=file.path(outdir_2,paste(colnames(bool_df)[i],"_ZK_Plot.pdf",sep="")))
    eset_temp <- eset[,bool_df[,i]]
    plot_ZK(eset_temp,label=label)
    dev.off()
    ZK_temp <- compute_ZK(eset_temp)
    remove <- colnames(eset_temp)[ZK_temp < -2]
    save(file=file.path(outdir_2,paste(colnames(bool_df)[i],"_To_Remove.Rda",sep="")),remove)
    list_removal_by_group[[i]] <- remove
  }
  names(list_removal_by_group) <- colnames(bool_df)
  
  return(list_removal_by_group)
}


array_quality_metrics_analysis <- function(eset,directory,tagg,do.logtransform = TRUE,pre_post = c("Pre_Norm")){
  if(do.logtransform){
    eset <- transfor_to_log(eset)
  }
  library("arrayQualityMetrics")
  outdir <- file.path(directory,"Results",pre_post)
  outdir_2 <- file.path(outdir,tagg)
  dir.create(outdir_2)
  dat <- arrayQualityMetrics(expressionset = eset, outdir = outdir_2, force = TRUE, do.logtransform = FALSE)
  save(file = file.path(outdir_2,"Raw_array_quality_metrics_Results.Rda"),dat)
}

list_esets_from_groups <- function(eset,var){
  groups <- unique(pData(eset)[,var])
  print(groups)
  list_out <- list()
  for(i in 1:length(groups)){
    print(i)
    print(groups[i])
    eset_temp <- eset[,pData(eset)[,var] == groups[i]]
    list_out[[i]] <- eset_temp
  }
  names(list_out) <- groups
  return(list_out)
}

arrayQualityMetrics_All_Groups <- function(eset,directory,var,do.logtransform = TRUE,pre_post ="Pre_Norm"){
  library("arrayQualityMetrics")
  eset_group_list <- list_esets_from_groups(eset,var)
  group_names <- names(eset_group_list)
  mapply(eset_group_list,FUN = array_quality_metrics_analysis,directory = directory,tagg=group_names,pre_post = pre_post,do.logtransform = do.logtransform)
}

#FUNCTION WITH MISTAKES !!!!


get_list_of_folders <- function(directory = getwd(),type_an = "Group",pre_post ="Pre_Norm"){
  dir_2 <- file.path(directory,"Results",pre_post)
  folders <- dir(dir_2)
  if(type_an == "Group"){
    folders <- folders[!grepl("Whole",folders)]
    folders <- folders[!grepl("Legacy",folders)]
  }
  else{
    folders <- c("Whole_Dataset")
  }
  print(folders)
  list_out_dfs <- list()
  rownames(list)
  name_for_rows <- c()
  for(i in 1:length(folders)){
    results_AQM_temp <- get(load(file=file.path(dir_2,folders[i],"Raw_array_quality_metrics_Results.Rda")))
    df_temp <- data.frame(matrix(0,nrow = length(results_AQM_temp$arrayTable$sampleNames),ncol = 9))
    colnames(df_temp) <- c("pCh_QC_Heat","pCh_QC_Heat_Th","pCh_QC_Box","pCh_QC_Box_Th","pCh_QC_MA","pCh_QC_MA_Th","pCh_QC_ZK","pCh_QC_ZK_Th","pCh_QC_Groups")
    rownames(df_temp) <- results_AQM_temp$arrayTable$sampleNames
    name_for_rows <- c(name_for_rows,results_AQM_temp$arrayTable$sampleNames)
    out_heat <- names(results_AQM_temp$modules$heatmap@outliers@which)
    out_heat_th <- results_AQM_temp$modules$heatmap@outliers@threshold
    df_temp[out_heat,1] <- 1
    df_temp[,2] <- rep(out_heat_th,nrow(df_temp))
    out_box <- names(results_AQM_temp$modules$boxplot@outliers@which)
    out_box_th <- results_AQM_temp$modules$boxplot@outliers@threshold
    df_temp[out_box,3] <- 1
    df_temp[,4] <- rep(out_box_th,nrow(df_temp))
    out_MA <- names(results_AQM_temp$modules$maplot@outliers@which)
    out_MA_th <- results_AQM_temp$modules$maplot@outliers@threshold
    df_temp[out_MA,5] <- 1
    df_temp[,6] <- rep(out_box_th,nrow(df_temp))
    print(file.path(dir_2,folders[i],paste(folders[i],"_To_Remove.Rda",sep="")))
    results_ZK_temp <- get(load(file=file.path(dir_2,folders[i],paste(folders[i],"_To_Remove.Rda",sep=""))))
    df_temp[results_ZK_temp,7] <- 1
    df_temp[,8] <- rep(-2,nrow(df_temp))
    df_temp[,9] <- rep(folders[i],nrow(df_temp)) 
    list_out_dfs[[i]] <- df_temp
    
  }
  names(list_out_dfs) <- folders
  df_out <- do.call("rbind",list_out_dfs)
  rownames(df_out) <- name_for_rows
  #return(list_out_dfs)
  return(df_out)
}

get_samples_to_remove <- function(analysis_results,min_numner_of_failed_tests = 3){
  return(rowSums(analysis_results[,c(1,3,5,7)]) >= 3)
}


#Extracting data from qq (Good Function)

get_outliers_info <- function(data_aqm){
  df_temp <- data.frame(matrix(0,nrow = length(data_aqm$arrayTable$sampleNames),ncol = 6))
  colnames(df_temp) <- c("pCh_QC_Heat","pCh_QC_Heat_Th","pCh_QC_Box","pCh_QC_Box_Th","pCh_QC_MA","pCh_QC_MA_Th")
  rownames(df_temp) <- data_aqm$arrayTable$sampleNames
  out_heat <- names(data_aqm$modules$heatmap@outliers@which)
  out_heat_th <- data_aqm$modules$heatmap@outliers@threshold
  df_temp[out_heat,1] <- 1
  df_temp[,2] <- rep(out_heat_th,nrow(df_temp))
  out_box <- names(data_aqm$modules$boxplot@outliers@which)
  out_box_th <- data_aqm$modules$boxplot@outliers@threshold
  df_temp[out_box,3] <- 1
  df_temp[,4] <- rep(out_box_th,nrow(df_temp))
  out_MA <- data_aqm$modules$maplot@outliers@which
  out_MA_th <- data_aqm$modules$maplot@outliers@threshold
  df_temp[out_MA,5] <- 1
  df_temp[,6] <- rep(out_MA_th,nrow(df_temp))
  return(df_temp)
}


#Creating the first quality measures: The ones that I used until now. Based on the mean IAC.

auto_quality <- function(eset,tagg,directory,pre_post = "Pre_Norm"){
  eset <- transfor_to_log(eset)
  library("lumi")
  library("dendextend")
  library("affycoretools")
  print(file.path(directory,"Results",pre_post))
  if(!dir.exists(file.path(directory,"Results",pre_post))){
    dir.create(file.path(directory,"Results",pre_post))
  }
  dir.create(file.path(directory,"Results",pre_post,"Legacy"))
  pdf(file=file.path(directory,"Results",pre_post,"Legacy",paste(tagg,"Quality_Measures.pdf",sep="")))
  lumi::density(eset,legend = FALSE) #Detecta de manera automática si está en base 2
  lumi::boxplot(eset) #Detecta de manera automática si está en base 2.
  group <- pData(eset)$pNm_Status
  affycoretools::plotPCA(eset,groups = factor(group,levels = c("1","2")),groupnames = c("Control","Disease"),x.coord = 0, y.coord = 0)
  IAC=cor(exprs(eset),use="p")
  hist(IAC[upper.tri(IAC)],sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)),breaks=100)
  cluster1=hclust(as.dist(1-IAC),method="average")
  dend <- as.dendrogram(cluster1)
  dend <- set(dend, "labels_cex", 0.4)
  groupCodes <- pData(eset)$pCh_Status
  colorCodes <- c(T="blue", NT="red")
  labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
  plot(dend)
  meanIAC=apply(IAC,2,mean)
  sdCorr=sd(meanIAC)
  numbersd=(meanIAC-mean(meanIAC))/sdCorr
  plot(numbersd)
  abline(h=-2)
  dev.off()
}

############################
#Functions to mount dataset#
############################

#Loading datasets

load_dataset_list <- function(studies,dir_root,norm_type = c("rma","frma")){
  files_to_load <- paste(dir_root,studies,"/Data/",studies,"_Norm.Rda",sep = "")
  print(files_to_load)
  list_studies <- list()
  for(i in 1:length(files_to_load)){
    if(norm_type == "frma"){
      list_studies[[i]] <- get(load(file = files_to_load[i]))[[2]]
    }else if(norm_type == "rma"){
      list_studies[[i]] <- get(load(file = files_to_load[i]))[[1]]
    }
  }
  names(list_studies) <- studies
  return(list_studies)
}

#Generating data.frame with pheno data.

generate_pData <- function(list_load){
  for(i in 1:length(list_load)){
    if(i == 1){
      pData_All <- pData(list_load[[i]])
    }else{
      pData_All <- rbind.fill(pData_All,pData(list_load[[i]]))
    }
  }
  return(pData_All)
}

#Merging all data.

merge_edata_based_on_shared <- function(list_load){
  list_of_probes <- lapply(list_load,function(x) rownames(x))
  Intersected_Probes <- Reduce("intersect",list_of_probes)
  for(i in 1:length(list_load)){
    if(i == 1){
      eData_All <- exprs(list_load[[i]])[Intersected_Probes,]
    }else{
      eData_All <- cbind(eData_All,exprs(list_load[[i]])[Intersected_Probes,])
    }
  }
  return(eData_All)
}

#Correcting batch effects.

correct_batch <- function(eData_All,pData_All){
  require(sva)
  Batch_Data <- pData_All$pCh_Batch
  Data_Mod <- model.matrix(~as.factor(pData_All$pCh_Status))
  eData_All_Combat = ComBat(dat=eData_All, batch=Batch_Data, mod=Data_Mod, par.prior=TRUE, prior.plots=FALSE)
  return(eData_All_Combat)
}

#Filtering out probes that do not target genes.

filt_by_gene_tar_probes <- function(eData_All_Combat){
  library("hgu133a.db")
  x <- hgu133aSYMBOL
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  eData_All_Combat <- eData_All_Combat[names(xx),]
  return(eData_All_Combat)
}

#Collapsing probes that target the same gene.

collapse_rows <- function(entrez_probe_map,e_set,Method = "MaxMean"){
  if(!(paste("package:","WGCNA",sep="") %in% search())){
    require("WGCNA")
  }
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  xx_filt <- xx[rownames(e_set)]
  if(class(e_set)[1] == "matrix"){
    exprs_eset_collapsed <- collapseRows(e_set,as.character(xx_filt),rownames(e_set),method=Method)
  }else if(class(e_set)[1] == "data.frame"){
    exprs_eset_collapsed <- collapseRows(exprs(e_set),as.character(xx_filt),rownames(e_set),method=Method)
  }
  exprs_eset_collapsed <- exprs_eset_collapsed[[1]][!is.na(rownames(exprs_eset_collapsed[[1]])),]
  return(exprs_eset_collapsed)
}

#Cleaning covariates.

Clean_pheno <- function(pData_All,norm_type = c("rma","frma")){
  #Survival block.
  try({
    pData_All$pCh_DFS_E[pData_All$pCh_DFS_E == "NA"] <- NA
  })
  try({
    pData_All$pCh_DFS_T <- as.numeric( pData_All$pCh_DFS_T)
  })
  try({
    pData_All$pCh_Death_E[pData_All$pCh_Death_E == "NA"] <- NA
  })
  try({
    pData_All$pCh_Death_T <- as.numeric( pData_All$pCh_Death_T)
  })
  #Receptors block
  try({
    pData_All$pCh_ER_Status[pData_All$pCh_ER_Status == "ER+"] <- "1"
  })
  try({
    pData_All$pCh_ER_Status[pData_All$pCh_ER_Status == "I"] <- NA
  })
  try({
    pData_All$pCh_ER_Status[pData_All$pCh_ER_Status == "NA"] <- NA
  })
  try({
    pData_All$pCh_PR_Status[pData_All$pCh_PR_Status == "NA"] <- NA
    
  })
  try({
    pData_All$pCh_HER2_Status[pData_All$pCh_HER2_Status == "I"] <- NA
  })
  try({
    pData_All$pCh_HER2_Status[pData_All$pCh_HER2_Status == "NA"] <- NA
  })
  try({
    pData_All$pCh_HER2_Status[is.na(pData_All$pCh_HER2_Status)] <- pData_All$pCh_ERBB2_Status[is.na(pData_All$pCh_HER2_Status)]
  })
  try({
    pData_All$pCh_HER2_Status[pData_All$pCh_HER2_Status == "NA"] <- NA
  })
  #Age
  try({
    pData_All$pCh_Age <- as.numeric(pData_All$pCh_Age)
  })
  #Grade
  try({
    pData_All$pCh_Grade[pData_All$pCh_Grade == "G1"] <- "1"
    
  })
  try({
    pData_All$pCh_Grade[pData_All$pCh_Grade == "G2"] <- "2"
    
  })
  try({
    pData_All$pCh_Grade[pData_All$pCh_Grade == "G3"] <- "3"
    
  })
  try({
    pData_All$pCh_Grade[pData_All$pCh_Grade == "4=Indeterminate"] <- NA
    
  })
  try({
    pData_All$pCh_Grade[pData_All$pCh_Grade == "NA"] <- NA
    
  })
  #Node
  try({
    pData_All$pCh_Node[pData_All$pCh_Node == "N2"] <- "1"
    
  })
  try({
    pData_All$pCh_Node[pData_All$pCh_Node == "N3"] <- "1"
    
  })
  try({
    pData_All$pCh_Node[pData_All$pCh_Node == "NA"] <- NA
  })
  if(norm_type == "rma"){
    try({
      pData_All$pam50[pData_All$pCh_Status == "NT"] <- "Healthy_breast"
      
    })
  }
  if(norm_type == "frma"){
    try({
      pData_All$pam50_frma[pData_All$pCh_Status == "NT"] <- "Healthy_breast"
    })
  }
  rownames(pData_All) <- pData_All$pCh_Sample_Name
  return(pData_All)
}


#Removing duplicated samples based on checksums.

remove_dup_check <- function(pData_All,eData_All_Combat){
  bool_filt <- !duplicated(pData_All$pCh_md5_to_df)
  eData_All_Combat <- eData_All_Combat[,bool_filt]
  pData_All <- pData_All[bool_filt,]
  return(list(eData_All_Combat,pData_All))
}

#Carry out cox proportional hazard models using one gene at a time as predictor.

#Flattening normal tissue vectors.

flattenign_vectors <- function(normal_tissue_data){
  df_out <- normal_tissue_data
  for(i in 1:ncol(normal_tissue_data)){
    df_out[,i] <- fitted(lm(normal_tissue_data[,i] ~ 0 + ., data = data.frame(normal_tissue_data)[,-i]))
  }
  return(df_out)
}

#wssplot for clustering

wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

#Generating disease vectors.

disease_vector_analysis <- function(eData_All_Combat,pData_All){
  bool_filt <- pData_All$pCh_Status == "NT"
  eData_All_Combat_NT <- eData_All_Combat[,bool_filt]
  bool_filt <- pData_All$pCh_Status == "T"
  eData_All_Combat_T <- eData_All_Combat[,bool_filt]
  flattened_vectors <- flattenign_vectors(eData_All_Combat_NT)
  #Performing dimensionality reduction.
  svd_results <- svd(flattened_vectors)
  g_w_out <- goodnes_w(svd_results$d,nrow(eData_All_Combat_NT))
  print("Done")
  plot(1:(ncol(eData_All_Combat_NT)-1),g_w_out,ylim  = c(1,50))
  number_of_components <- readline(prompt="Enter the selected number of components: ")
  n_esp <- Create_Normal_State_Model(svd_results,as.numeric(number_of_components))
  #n_esp <- n_esp[,1:as.numeric(number_of_components)]
  disease_vectors <- generate_the_disease_vector(eData_All_Combat,n_esp)
  return(disease_vectors)
}

#Generating the disease vectors. 

generate_the_disease_vector <- function(disease_matrix,normal_space){
  output_df <- disease_matrix
  for(i in 1:ncol(disease_matrix)){
    print(i)
    output_df[,i] <- resid(lm(disease_matrix[,i] ~ 0 + ., data = data.frame(normal_space)))
  }
  return(output_df)
}


#Compute the goodness of fit measure for W for for the flat matrix.

goodnes_w <- function(vector_of_sv,n_val){
  output_vector <- c()
  for(i in 1:(length(vector_of_sv)-1)){
    left_eq_num <- vector_of_sv[i]^2
    #print(left_eq_num)
    left_eq_den <- sum((vector_of_sv[(i+1):length(vector_of_sv)])^2)
    #print(left_eq_den)
    left_side <- left_eq_num/left_eq_den
    #print(left_side)
    right_side_num <- (n_val-i-1)*(length(vector_of_sv)-i)
    #print(right_side_num)
    right_side_den <- (n_val + length(vector_of_sv)-(2*i))
    #print(right_side_den)
    right_side <- right_side_num/right_side_den
    #print(right_side)
    W_l <- left_side*right_side
    output_vector <- c(output_vector,W_l)
  }
  return(output_vector)
}

#Create the normal state model.

Create_Normal_State_Model <- function(svd_res,number_of_selected_comp){
  matrix_b <- matrix(0,nrow = length(svd_res$d),ncol = length(svd_res$d))
  vector_diag <- svd_res$d
  vector_diag[(number_of_selected_comp+1):length(vector_diag)] <- 0
  print(vector_diag)
  diag(matrix_b) <- vector_diag
  print(matrix_b)
  normal_space_a <- svd_res$u %*% matrix_b 
  print("This is the dimension of the first matrix multiplication")
  print(dim(normal_space_a))
  print("This is the dimension of the second matrix multiplication")
  normal_space_b <- normal_space_a %*% t(svd_res$v)
  print(dim(normal_space_b))
  return(normal_space_b)
}

#Computing proportional hazard models for each gene.


cox_all_genes <- function(eData,pData){
  require(dplyr)
  require(survival)
  require(survminer)
  list_out <- list()
  for(i in 1:nrow(eData)){
    print(i)
    temp <- summary(coxph(Surv(pData$pCh_DFS_T,as.numeric(pData$pCh_DFS_E))~eData[i,]))$coefficients[1,]
    list_out[[i]] <- temp
  }
  df_out <- data.frame(do.call("rbind",list_out))
  rownames(df_out) <-rownames(eData)
  return(df_out)
}


#Plotting WGCNA functions (fit to power law)

plot_wgcna_power_Tables <- function(x,name_study,thr_fit = 0.85){
  plot(x$fitIndices[,1], -sign(x$fitIndices[,3])*x$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence in ", name_study),ylim = c(0,1))
  text(x$fitIndices[,1], -sign(x$fitIndices[,3])*x$fitIndices[,2],
       labels=powers,cex=1,col="red")
  abline(h=thr_fit,col="red")
}

#Plotting WGCNA functions (Connectivity)


plot_wgcna_power_Tables_con <- function(x,name_study,thr_fit = 0.85){
  plot(x$fitIndices[,1], x$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence in ", name_study))
  text(x$fitIndices[,1], x$fitIndices[,5],
       labels=powers,cex=1,col="red")
}

