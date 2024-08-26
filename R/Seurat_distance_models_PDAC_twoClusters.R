# Author: Anna Lyubetskaya. Date: 20-04-22

# This script compares normalized gene expression values in each spot to a select spot distance metric (e.g., distance between a spot and tumor center)
# The comparison is performed in the following ways:
# - Correlation between distance and each gene expression profile across spots
# - Linear regression fit: GeneX expression ~ Distance + Cell counts
# Differences between tissue sections can be accounted by including section identity as a model covariate or by normalizing the distance metric to within the section


## ENVIRONMENT ----


library(Seurat)
library(foreach)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


source("code/utils/utils_ggplot.R")
source("code/utils/utils_signatures.R")
source("code/utils/utils_lin_regress.R")
source("code/utils/utils_stats.R")

source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


# Path to processed Seurat data
sample_name <- "a_name"
sample_exclude <- NULL  

# Distance variable to use for correlations
dist_var <- "DistanceBySpot"  # Pathology.Distance.Tissue.filled; Pathology.Distance.Epithelium.filled.invertTRUE
# Normalize pathology distances by sample
normalize_dist_by_sample <- FALSE

# Variables to include in the ElNet model as confounders
confounder_vars_num <- NULL #" # c("CellCounts")
confounder_vars_cat <- c("user.Sample_Name", "stromal_subcluster_name")  # e.g., NULL or c("user.Sample_Name", "user.Tissue")

# Filters to select genes
sct_threshold <- 0.5
spot_threshold <- 5

# distance cutoff 
distance_thres <- 10 # fullres pixel thres 2000 # NULL

# cluster filter
keep_cluster <- c("Fibroblasts", "myCAF", "iCAF", "ISG.RS") # NULL
keep_cluster_by <- "stromal_subcluster_name"
reference_cluster <- c("PDACBasal_Hypoxia", "PDACClassical", "Ducts") # The order of the reference matters when fetching the overlay results

# List of genes to highlight in the scatter plot
gene_highlight_list <- NULL


## PATHS ----


# Path to a signature file
sig_path <- "dir/to/siganturefile.txt"

# Location of pre-processed data
input_path <- paste0("/dir/to/integrated/file/", sample_name, ".rds")
distance_file_path <- "/dir/to/distancefile/output/from/previous/script/"
collagen_overlay_path <- "/dir/to/collagen/overlay/results/"

# Output folder
output_path_init <- "/output/path/"
output_path <- paste0(output_path_init, confounder_vars_num, "Reference_", paste0(reference_cluster,collapse = "_"), "_Target_", paste0(keep_cluster,collapse = "_"), "_", dist_var, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# full cluster
cluster_tab <- read.delim("/dir/to/fullcohort_clustering.txt")
cluster_name <- read.delim("/dir/to/dictionary_full_cohort.txt")
cluster_dict <- hash::hash(keys = cluster_name$fullcohort_cluster, values = cluster_name$fullcohort_cluster_name)
cluster_tab$fullcohort_cluster_name <- hash::values(cluster_dict, keys = cluster_tab$fullcohort_cluster)

# epi cluster
cluster_tab_epi <- read.delim("/dir/to/epithelium_cohort_clustering.txt")
cluster_name_epi <- read.delim("/dir/to/dictionary_epi.txt")
cluster_dict_epi <- hash::hash(keys = cluster_name_epi$epi_subcluster, values = cluster_name_epi$epi_subcluster_name)
cluster_tab_epi$epi_subcluster_name <- hash::values(cluster_dict_epi, keys = cluster_tab_epi$epi_subcluster)

# stroma cluster
cluster_tab_stroma <- read.delim("/dir/to/stromal_cohort_clustering.txt")
cluster_name_stroma <- read.delim("/dir/to/dictionary_stromal.txt")
cluster_dict_stroma <- hash::hash(keys = cluster_name_stroma$stromal_subcluster, values = cluster_name_stroma$stromal_subcluster_name)
cluster_tab_stroma$stromal_subcluster_name <- hash::values(cluster_dict_stroma, keys = cluster_tab_stroma$stromal_subcluster)

# merge
cluster_tab <- cluster_tab %>% 
  dplyr::select(c("Coordinate","Pathology.Group","fullcohort_cluster","fullcohort_cluster_name")) %>% 
  dplyr::left_join(cluster_tab_epi %>% dplyr::select(c("Coordinate","epi_subcluster","epi_subcluster_name")), by = "Coordinate") %>% 
  dplyr::left_join(cluster_tab_stroma %>% dplyr::select(c("Coordinate","stromal_subcluster","stromal_subcluster_name")), by = "Coordinate") #%>% 
#dplyr::left_join(data_seurat@meta.data %>% tibble::rownames_to_column("Coordinate") %>% dplyr::select(c("Coordinate", "Sample_Name"="user.Sample_Name")), by = "Coordinate")

data_seurat@meta.data <- data_seurat@meta.data %>% tibble::rownames_to_column("Coordinate") %>% dplyr::left_join(cluster_tab, by = "Coordinate") %>% tibble::column_to_rownames("Coordinate")

# Find samples in the merged object that don't satisfy the gene selection criteria
sample_exclude <- c(sample_exclude, names(table(data_seurat@meta.data$user.Sample_Name))[which(table(data_seurat@meta.data$user.Sample_Name) <= spot_threshold)])

# Remove samples if necessary
if(length(sample_exclude)!=0){
  barcode_list <- rownames(data_seurat@meta.data[which(!data_seurat@meta.data[["user.Sample_Name"]] %in% sample_exclude),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}
if(!is.null(keep_cluster)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[[keep_cluster_by]] %in% keep_cluster),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# List of samples and number of spots in each
table(data_seurat@meta.data$user.Sample_Name)

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, NULL, sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Invert signatures to get an annotated gene list
sig_gene_df <- invert_list_my(signature_list)


## WRANGLE DATA ----


# Read in Distance 
for (ref in reference_cluster){
  file_dir <- list.files(distance_file_path, paste0("table_", ref, ".txt"), full.names = T)
  distance_tmp <- read.delim(file_dir)
  distance_tmp <- distance_tmp %>% dplyr::select(c(!!rlang::sym(dist_var), "Coordinate"))
  colnames(distance_tmp)[1] <- paste0(ref, "_", dist_var)
  data_seurat@meta.data <- data_seurat@meta.data %>% tibble::rownames_to_column("Coordinate") %>% dplyr::left_join(distance_tmp, by = "Coordinate") %>% tibble::column_to_rownames("Coordinate")
}

# do this step before abundant gene identification to keep same gene set.
if(!is.null(distance_thres)){
  dist_cols <- data_seurat@meta.data[,stringr::str_detect(colnames(data_seurat@meta.data), dist_var)]
  barcode_list <- rownames(dist_cols)[rowSums(dist_cols <= distance_thres) > 0] # this way of filtering results in Inf distance for certain samples, since subclusters could be missing there.
  #barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data$Distance <= distance_thres),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Abundant gene list
gene_list <- unique(seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold, 
                                                    assay="SCT", slot="data"))

# Exclude MT/Ribo genes for ease of interpretation
gene_list <- setdiff(gene_list, gene_list[grep("^MT-|^RP[SL]|^GM\\d+", gene_list)])
# Remove NA value
gene_list <- gene_list[which(!is.na(gene_list))]

# Remove all distance cols
data_seurat@meta.data <- data_seurat@meta.data[,!(colnames(data_seurat@meta.data) %in% colnames(dist_cols))]
coef_gene_df_cumu <- data.frame()

for(ref in reference_cluster){
  print(ref)
  output_path_ref <- paste0(output_path, ref, "/")
  
  # Create output folders
  dir.create(output_path_ref, showWarnings = FALSE)
  
  data_seurat_tmp <- data_seurat
  
  file_dir <- list.files(distance_file_path, paste0("table_", ref, ".txt"), full.names = T)
  distance_tmp <- read.delim(file_dir)
  distance_tmp <- distance_tmp %>% dplyr::select(c(!!rlang::sym(dist_var), "Coordinate"))
  data_seurat_tmp@meta.data <- data_seurat_tmp@meta.data %>% tibble::rownames_to_column("Coordinate") %>% dplyr::left_join(distance_tmp, by = "Coordinate") %>% tibble::column_to_rownames("Coordinate")
  
  # Remove spots too far away
  if(!is.null(distance_thres)){
    dist_cols <- data_seurat_tmp@meta.data[,stringr::str_detect(colnames(data_seurat_tmp@meta.data), dist_var)]
    #barcode_list <- rownames(dist_cols)[rowSums(dist_cols <= distance_thres) > 0] # this way of filtering results in Inf distance for certain samples, since subclusters could be missing there.
    barcode_list <- rownames(data_seurat_tmp@meta.data[which(data_seurat_tmp@meta.data[dist_var] <= distance_thres),])
    data_seurat_tmp <- subset(data_seurat_tmp, cells=barcode_list)
  }
  
  # Extract a wide matrix of SCT normalized expression values
  data_wide_df <- tibble::as_tibble(t(as.matrix(Seurat::GetAssayData(data_seurat_tmp, assay="SCT", slot="data"))), rownames="Coordinate")
  
  # Subset data to only relevant genes
  data_wide_df <- data_wide_df %>%
    dplyr::select(dplyr::all_of(intersect(c("Coordinate", gene_list), colnames(data_wide_df))))
  
  # Check that metadata and expression data have the same order our of paranoia
  table(rownames(data_seurat_tmp@meta.data) == data_wide_df$Coordinate)
  
  # Extract the relevant meta data
  meta_data <- data_seurat_tmp@meta.data[c(dist_var, confounder_vars_num, confounder_vars_cat)] %>%
    dplyr::filter(!is.na(!!rlang::sym(dist_var)) & !is.infinite(!!rlang::sym(dist_var))) %>%
    tibble::rownames_to_column("Coordinate")
  
  # Re-normalize pathology distance by sample (% of max)
  if(normalize_dist_by_sample == TRUE){
    dist_norm_df <- data_seurat_tmp@meta.data[c(dist_var, "user.Sample_Name")] %>%
      dplyr::group_by(user.Sample_Name) %>%
      dplyr::mutate(MaxDist = max(!!rlang::sym(dist_var)),
                    DISTANCE = !!rlang::sym(dist_var) / MaxDist * 100)
    
    meta_data[dist_var] <- dist_norm_df$DISTANCE
  }
  
  # Add meta data to the expression data
  data_wide_df <- data_wide_df %>%
    dplyr::inner_join(meta_data, by="Coordinate")
  
  
  ## FIND DISTANCE CORRELATES ----
  
  
  # Calculate correlation between the distance and each of the gene expression vectors
  cor_vec <- sapply(gene_list, function(x) cor(data_wide_df[[x]], meta_data[[dist_var]]))
  names(cor_vec) <- gene_list
  
  # Create a correlation tibble
  cor_df <- tibble::tibble(Symbol = gene_list,
                           R2 = cor_vec) %>%
    dplyr::mutate(R2 = round(R2, 3)) %>%
    dplyr::arrange(R2) %>%
    dplyr::left_join(sig_gene_df, by="Symbol")
  
  # Write correlations to file
  filename <- paste0(output_path_ref, "Distance2", ref, "_correlation_result.txt")
  readr::write_delim(cor_df, filename, delim="\t")
  
  # Create a correlation histogram
  filename <- paste0(output_path_ref, "Distance2", ref, "_correlation_hist")
  create_hist_plot_my(cor_df, x_label="R2", fill_label="InSignature", intercept=c(-0.25, 0, 0.25), 
                      binwidth=0.01, filename=filename, 
                      labels=c("R2", "Gene number", paste0("Correlation with", dist_var)))
  
  # Cleanup before running more analysis
  rm(data_seurat_tmp)
  gc()
  
  
  ## LM OF EXPRESSION ON DISTANCE ----
  
  
  ## Prep data for LM
  
  # Define the outcome and explicitly define row names
  data_wide_df <- data_wide_df %>% 
    dplyr::select(-Coordinate)
  
  # Factorize factor fields
  if(!is.null(confounder_vars_cat) && length(confounder_vars_cat) > 0){
    data_wide_df[confounder_vars_cat] <- lapply(data_wide_df[confounder_vars_cat], as.factor) 
  }
  
  
  # Save gene LMs to a file as a list of models
  filename <- paste0(output_path_ref, "Distance2", ref, "Vgene_LM_result.rds")
  gene_model_list <- list()
  
  # Go through all genes and perform individual fits: Gene Expression ~ Distance to Center
  if(!file.exists(filename)){
    
    # Parallelize the process for speed
    cl <- parallel::makeCluster(parallel::detectCores() / 2)
    doParallel::registerDoParallel(cl)
    
    gene_model_list <- foreach(gene = gene_list, .combine='c') %dopar% {
      
      # Define the outcome and explicitly define row names
      data_wide_loc_df <- data_wide_df %>%
        dplyr::select(dplyr::all_of(c(dist_var, gene, confounder_vars_cat, confounder_vars_num))) %>% 
        dplyr::rename("outcome" = !!rlang::sym(gene))
      
      # Run LM
      model_res <- run_all_models_my(data_wide_loc_df, c("norm"), output_loc=NULL, formula=NULL)
      list(gene = summary(model_res$norm)$coefficients)
      
    }
    
    parallel::stopCluster(cl)
    
    # Add gene names to list names
    names(gene_model_list) <- gene_list
    
    # Save a list of models to file
    saveRDS(gene_model_list, filename)
  } else{
    gene_model_list <- readRDS(filename)
  }
  
  # Analyze the result of the linear model
  coef_gene_df <- t(sapply(names(gene_model_list), function(g) gene_model_list[[g]][dist_var, c(1,4)])) %>%
    tibble::as_tibble(rownames="Symbol") %>%
    dplyr::mutate(Score = -log10(`Pr(>|t|)`),
                  Labels = ifelse(Symbol %in% gene_highlight_list,
                                  Symbol, ""),
                  ref_cluster = ref)
  
  # Establish the estimate threshold as mean +/- 3SD
  estimate_threshold <- mean(abs(coef_gene_df$Estimate)) + sd(abs(coef_gene_df$Estimate)) * 3
  
  # Call significant coefficients
  coef_gene_df <- coef_gene_df %>%
    dplyr::mutate(IsSignificant = `Pr(>|t|)` <= 0.05 / length(gene_list) & abs(Estimate) >= estimate_threshold) %>%
    dplyr::arrange(desc(IsSignificant), desc(Score))
  
  # Visualize the result of linear modeling of gene expression ~ distance to center
  filename <- paste0(output_path_ref, "Distance2", ref, "Vgene_LM_result_scatter")
  create_scatter_plot_my(coef_gene_df, x_label="Estimate", y_label="Score", 
                         fill_label="IsSignificant", shape=21, size=2, dot_labels="Labels", 
                         filename=filename, labels=NULL, do_fit=NULL, stroke=0)
  coef_gene_df$IsSignificant <- factor(coef_gene_df$IsSignificant, levels = c( "TRUE", "FALSE"))
  volcano_plot <- ggplot(coef_gene_df, aes(Estimate, Score, color = IsSignificant))+
    geom_point()+theme_bw()+
    geom_hline(yintercept=-log10(0.05), color = "red")+
    geom_vline(xintercept = c(-estimate_threshold,estimate_threshold),color ="red")+
    scale_color_manual(values=c("red", "black"))+ggrepel::geom_text_repel(data = subset(coef_gene_df, IsSignificant == "TRUE"),
                                                                          aes(label = Symbol),
                                                                          size = 3)+ggtitle(paste0("Gene Expr ~ Distance to ", ref,", adjusted for ", paste0(c(confounder_vars_cat, confounder_vars_num), collapse = "+")))
  png(paste0(output_path_ref, "Volcano_LM_Distance2", ref,".png"), width = 9, height = 7.3, res = 300, units = "in")
  print(volcano_plot)
  dev.off()
  
  coef_gene_df <- coef_gene_df %>% dplyr::left_join(sig_gene_df %>% dplyr::select(c("Symbol", "Sig_Name", "Num_Sigs")), by="Symbol")
  
  filename <- paste0(output_path_ref, "Distance2", ref, "Vgene_LM_result_scatter.txt")
  readr::write_delim(coef_gene_df, filename, delim="\t")
  
  write.table(coef_gene_df$Symbol[coef_gene_df$IsSignificant == T & coef_gene_df$Estimate >0], paste0(output_path_ref, "gene_posAssociated_Distance.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
  write.table(coef_gene_df$Symbol[coef_gene_df$IsSignificant == T & coef_gene_df$Estimate <0], paste0(output_path_ref, "gene_negAssociated_Distance.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
  
  coef_gene_df_cumu <- rbind(coef_gene_df_cumu, coef_gene_df)
}
filename <- paste0(output_path, "Distance2", paste0(reference_cluster, collapse = "_"), "_cumu_LM_result_scatter.txt")
readr::write_delim(coef_gene_df_cumu, filename, delim="\t")


## Visualize pair-wise comparison as scatter plot
possible_combn <- combn(reference_cluster,2)
for(i in 1:ncol(possible_combn)){
  cluster1 <- stringr::str_replace(possible_combn[1,i], "/","_")
  cluster2 <- stringr::str_replace(possible_combn[2,i], "/","_")
  
  collagen_overlay_files <- list.files(dir(collagen_overlay_path, pattern = paste0(cluster1, "v", cluster2), full.names = T), full.names = T)
  collagen_overlay_file <- collagen_overlay_files[stringr::str_detect(collagen_overlay_files, "gene_LM_result_scatter_[:print:]*.txt$")]
  collagen_overlay <- readr::read_delim(collagen_overlay_file)
  
  coef_gene_df_cumu_tmp <- coef_gene_df_cumu %>% dplyr::filter(ref_cluster %in% c(cluster1, cluster2)) %>% 
    dplyr::select(c(Symbol, Estimate, ref_cluster, IsSignificant))
  coef_gene_df_cumu_tmp <- reshape2::dcast(reshape2::melt(coef_gene_df_cumu_tmp, id.vars = c("Symbol","ref_cluster")), Symbol~variable+ref_cluster)
  coef_gene_df_cumu_tmp[,stringr::str_detect(colnames(coef_gene_df_cumu_tmp),"Estimate")] <- apply(coef_gene_df_cumu_tmp[,stringr::str_detect(colnames(coef_gene_df_cumu_tmp),"Estimate")], 2, as.numeric)
  
  coef_gene_df_cumu_tmp$diff <- coef_gene_df_cumu_tmp[[paste0("Estimate_", cluster2)]] - coef_gene_df_cumu_tmp[[paste0("Estimate_", cluster1)]]
  diff_mean <- mean(coef_gene_df_cumu_tmp$diff)
  diff_sd <- sd(coef_gene_df_cumu_tmp$diff)
  coef_gene_df_cumu_tmp$diff_z = (coef_gene_df_cumu_tmp$diff - diff_mean)/diff_sd 
  
  group3sd <- character(length = nrow(coef_gene_df_cumu_tmp))
  for(rowi in 1:nrow(coef_gene_df_cumu_tmp)){
    diff_tmp <- coef_gene_df_cumu_tmp$diff[rowi]
    if(diff_tmp > (diff_mean+3*diff_sd)){
      group3sd[rowi] <- paste0(cluster1, "_outlier")
    }else if(diff_tmp < (diff_mean-3*diff_sd)){
      group3sd[rowi] <- paste0(cluster2, "_outlier")
    }else{
      group3sd[rowi] <- "fit"
    }
  }
  coef_gene_df_cumu_tmp$IsOutlier <- group3sd
  
  # Force more stringent cutoff to call significant association with collagen abs(Estimate) >= 0.5 & IsSignificant == T
  collagen_overlay <- collagen_overlay %>% 
    dplyr::mutate(IsSignificant = ifelse(abs(Estimate) >= 0.3 & IsSignificant == T, T, F))
  coef_gene_df_cumu_tmp <- coef_gene_df_cumu_tmp %>% dplyr::left_join(collagen_overlay %>% 
                                               dplyr::select(c(Symbol, Estimate, Pval = `Pr(>|t|)`, Collagen_Sig = IsSignificant)), by = "Symbol") %>% 
    dplyr::mutate(Pfdr = p.adjust(Pval, method="fdr"))
  
  ## collagen association categorized by significance
  ## Genes whose association with collagen with abs(Estimate) >= 0.3 & Befferoni adjusted p-value < 0.05 colored
  ## All others grey
  ## Names of genes with significant association with collagen/outlier in distance models are labeled
  coef_gene_df_cumu_tmp$Collagen <- NA
  coef_gene_df_cumu_tmp$Collagen[coef_gene_df_cumu_tmp$Collagen_Sig == F] <- "NS"
  coef_gene_df_cumu_tmp$Collagen[coef_gene_df_cumu_tmp$Collagen_Sig == T & coef_gene_df_cumu_tmp$Estimate > 0] <- "Positive"
  coef_gene_df_cumu_tmp$Collagen[coef_gene_df_cumu_tmp$Collagen_Sig == T & coef_gene_df_cumu_tmp$Estimate < 0] <- "Negative"
  
  outlier_tab <- coef_gene_df_cumu_tmp %>% 
    dplyr::mutate(xplusy=abs(!!rlang::sym(paste0("Estimate_", cluster1)))+abs(!!rlang::sym(paste0("Estimate_", cluster2)))) %>%
    dplyr::slice_max(order_by = xplusy, n = 30)
  outlier_tab <- outlier_tab %>% 
    dplyr::filter(!!rlang::sym(paste0("IsSignificant_", cluster1)) ==1|!!rlang::sym(paste0("IsSignificant_",cluster2)) == 1)
  
  tmp_plot <- ggplot(coef_gene_df_cumu_tmp %>% dplyr::arrange(Collagen), aes(!!rlang::sym(paste0("Estimate_", cluster1)), !!rlang::sym(paste0("Estimate_", cluster2)), color = Collagen, size = Collagen))+
    geom_point()+
    scale_color_manual(breaks = c("NS", "Coll_NS_Dist_outlier", "Negative", "Positive"), values = c("grey", "black", "royalblue", "#c10001"))+
    scale_size_manual(breaks = c("NS", "Coll_NS_Dist_outlier", "Negative", "Positive"), values = c(0.5, 1.2, 1.2, 1.2))+
    ggrepel::geom_text_repel(data = coef_gene_df_cumu_tmp %>% 
                               #dplyr::filter(IsOutlier!="fit"),
                               dplyr::filter(Symbol %in% outlier_tab$Symbol, coef_gene_df_cumu_tmp$Symbol[coef_gene_df_cumu_tmp$Collagen != "NS"]), 
                             aes(label = Symbol),size = 3, fontface = "italic")+
    geom_abline(aes(slope=1,intercept=0)) + 
    geom_abline(aes(slope=1,intercept=0 + 3*diff_sd),linetype=3) +
    geom_abline(aes(slope=1,intercept=0 - 3*diff_sd),linetype=3) +
    theme_bw()+labs(color = "Association with Collagen", x = paste0("Estimated association with distance to ", cluster1), y = paste0("Estimated association with distance to ", cluster2))
  pdf(paste0(output_path, "path_LM_result_overlay_collagen_", cluster2, "v", cluster1,".pdf"), width = 8, height = 6)
  plot(tmp_plot)
  dev.off()
  
  filename <- paste0(output_path, "ScatterPlot_", cluster2, "v", cluster1, "_with_Collagen_association.txt")
  readr::write_delim(coef_gene_df_cumu_tmp, filename, delim="\t")
}
