# Author: Anna Lyubetskaya. Date: 20-06-10


read_10X_H5_my <- function(input_path){
  ## Ingest a file to a Seurat object
  
  # Load GEX data
  data_seurat <- Seurat::CreateSeuratObject(counts = Seurat::Read10X_h5(input_path))
  
  # Make sure all gene names are in upper case - important for feature subsetting and matching
  # It's better to do it as early as possible
  rownames(data_seurat@assays$RNA@counts) <- toupper(rownames(data_seurat@assays$RNA@counts))
  rownames(data_seurat@assays$RNA@data) <- toupper(rownames(data_seurat@assays$RNA@data))
  rownames(data_seurat@assays$RNA@meta.features) <- toupper(rownames(data_seurat@assays$RNA@meta.features))
  
  return(data_seurat)
}


read_10X_spatial_folder_my <- function(folder_10X, filename="filtered_feature_bc_matrix.h5"){
  ## Read a 10X Spatial folder 
  
  data_seurat <- Seurat::Load10X_Spatial(folder_10X, filename=filename)
  
  # Make sure all gene names are in upper case - important for feature subsetting and matching
  # It's better to do it as early as possible
  rownames(data_seurat@assays$Spatial@counts) <- toupper(rownames(data_seurat@assays$Spatial@counts))
  rownames(data_seurat@assays$Spatial@data) <- toupper(rownames(data_seurat@assays$Spatial@data))
  rownames(data_seurat@assays$Spatial@meta.features) <- toupper(rownames(data_seurat@assays$Spatial@meta.features))
  
  #if(filename == "raw_feature_bc_matrix.h5" && ncol(data_seurat) != nrow(data_seurat@images$slice1@coordinates)){
  #  cells <- intersect(rownames(data_seurat@images$slice1@coordinates), colnames(data_seurat))
  #  data_seurat@images$slice1@coordinates <- data_seurat@images$slice1@coordinates[cells,]
  #}
  
  return(data_seurat) 
}


read_10X_qc_metrics_my <- function(folder_10X){
  ## Return default 10X QC metrics
  
  # Define metrics file full path
  qc_files <- paste0(folder_10X, "/metrics_summary.csv")
  
  qc_df <- NULL
  if(file.exists(qc_files)){
    qc_df <- readr::read_delim(qc_files, delim=",")
  }
  
  return(qc_df)
}


read_rds_list_simple_my <- function(file_list, feature_filter=100, sig_clean=TRUE){
  ## Ingest a set of RDS Seurat objects
  
  # Go through list of files
  seurat_list <- list()
  for(rds_file in file_list){

    # Ingest a Seurat object
    data_seurat <- readRDS(rds_file)
    
    # Identify sample name
    sample_name <- data_seurat@misc$user.Sample_Name
    # In case this object hasn't been given a sample name yet
    if(is.null(sample_name)){
      sample_name <- gsub("^.+/|\\.[^\\.]+$", "", rds_file) 
    }
    
    # Only select barcodes with enough features
    if(!is.null(feature_filter)){
      # Select barcodes by their number of features
      if("nFeature_SCT" %in% colnames(data_seurat@meta.data)){
        barcodes <- rownames(data_seurat@meta.data)[which(data_seurat@meta.data$nFeature_Spatial >= feature_filter)]
      } else{
        barcodes <- rownames(data_seurat@meta.data)
      }
      
      data_seurat <- subset(data_seurat, cells=barcodes)
    }
    
    # Remove any old signature information
    if(sig_clean == TRUE){
      data_seurat@meta.data <- data_seurat@meta.data[which(!grepl("sig.", colnames(data_seurat@meta.data)))]
    }
    
    seurat_list[sample_name] <- data_seurat
  }
  
  return(seurat_list)
}


read_probe_offtargets <- function(probe_file){
  ## Read and format offtarget probe data
  
  probe_df <- readr::read_delim(probe_file, delim=",") %>%
    dplyr::mutate(probe_id = toupper(gsub("^[^\\|]+\\||\\|[^\\|]+$", "", probe_id)),
                  off_target_genes = toupper(gsub("ENSG\\d+\\|", "", off_target_genes))) %>%
    unique() %>%
    dplyr::group_by(probe_id) %>%
    dplyr::summarize(off_target_genes = paste0(sort(unique(strsplit(off_target_genes, ";")[[1]])), collapse=";"),
                     off_target_count = length(strsplit(off_target_genes, ";")[[1]]),
                     off_target_genes_short = gsub(";*A[A-Z]\\d+\\.\\d;*", "", off_target_genes),
                     gene_name = paste0(c(unique(probe_id), ":", strsplit(off_target_genes_short, ";")[[1]][1], ":", off_target_count), collapse="")) %>%
    dplyr::ungroup()
  
  
  return(probe_df)
}


restore_filtered_probes <- function(data_seurat_raw, data_seurat, probe_file, do_rename=FALSE){
  ## Gather both raw and filtered 10x H5 object, restore missing probes from raw to filtered object, and rename the probes based on their offtargets
  
  # Read in offtarget probe data
  probe_df <- read_probe_offtargets(probe_file) %>%
    dplyr::mutate(probe_id = toupper(probe_id))
  
  # Subset the Seurat object
  data_seurat <- subset(data_seurat_raw, cells=colnames(data_seurat), features=unique(c(rownames(data_seurat), probe_df$probe_id)))
  
  if(do_rename == TRUE){
    # Update gene names for off-targets
    gene1_list <- rownames(data_seurat_raw@assays$Spatial@counts)
    #gene2_list <- rownames(data_seurat_raw@assays$Spatial@data)
    #gene3_list <- rownames(data_seurat_raw@assays$Spatial@meta.features)
    #table(gene1_list == gene2_list)
    #table(gene1_list == gene3_list)
    gene_df <- tibble::tibble(Ref = gene1_list) %>%
      dplyr::left_join(probe_df, by=c("Ref" = "probe_id")) %>%
      dplyr::mutate(gene_name = ifelse(is.na(gene_name), Ref, gene_name)) %>%
      unique()
    
    data_seurat_raw@assays$Spatial@counts@Dimnames[[1]] <- gene_df$gene_name
    data_seurat_raw@assays$Spatial@data@Dimnames[[1]] <- gene_df$gene_name
    data_seurat_raw@assays$Spatial@meta.features <- data.frame(matrix(ncol = 0, nrow = nrow(gene_df)), row.names = gene_df$gene_name)
    
    
    gene1_list <- rownames(data_seurat@assays$Spatial@counts)
    #gene2_list <- rownames(data_seurat@assays$Spatial@data)
    #gene3_list <- rownames(data_seurat@assays$Spatial@meta.features)
    #table(gene1_list == gene2_list)
    #table(gene1_list == gene3_list)
    gene_df <- tibble::tibble(Ref = gene1_list) %>%
      dplyr::left_join(probe_df, by=c("Ref" = "probe_id")) %>%
      dplyr::mutate(gene_name = ifelse(is.na(gene_name), Ref, gene_name))
    
    data_seurat@assays$Spatial@counts@Dimnames[[1]] <- gene_df$gene_name
    data_seurat@assays$Spatial@data@Dimnames[[1]] <- gene_df$gene_name
    data_seurat@assays$Spatial@meta.features <- data.frame(matrix(ncol = 0, nrow = nrow(gene_df)), row.names = gene_df$gene_name)
  }
  
  
  return(data_seurat)  
}
