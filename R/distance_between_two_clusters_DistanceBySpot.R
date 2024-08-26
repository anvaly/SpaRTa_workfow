# Measure distance between two signature-defined features


## SETUP ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_signatures.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_heatmap.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_image.R")
source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_signatures.R")


## PARAMETERS ----


# The cohort of interest regex ID
# Seurat RDS files are tagged as follows

# Following lines specify filters for subclusters at each compartment level. It aims to exclude subclusters included in targeted full cohort clusters but not of focus, like nerve in stromal compartment. 
# These are all subsets of Cluster1 and Cluster2
target_clusters_col <- "epi_subcluster_name"
target_clusters <- c("Ducts", "PanINs", "PDAC_Classical", "PDAC_Proliferative", "PDAC_Intermediate", "PDAC_Basal", "PDAC_FibHigh") # c("Ducts") c("PDACBasal/Hypoxia")
combine_target_clusters <- FALSE # TRUE means consider all target_clusters subclusters as a whole when caculating distance, FALSE means iterate over all subclusters
target_shortname <- "Tumor"

reference_clusters_col <- "stromal_subcluster_name"
reference_clusters <- c("Fibroblasts", "myCAF", "iCAF", "ISG.RS", "Stellate1", "Stellate2", "Endothelial")
reference_shortname <- "Fibs"
#DE_cluster <- c("Fibroblasts", "myCAF", "iCAF", "ISG.RS")

reference_name <- ifelse(length(reference_clusters)==1, stringr::str_replace(reference_clusters, "/","_"), reference_shortname)
target_name <- ifelse(length(target_clusters)==1, stringr::str_replace(target_clusters, "/","_"), target_shortname)

cohort_name <- paste0("PDAC_naive_", reference_name, "2", target_name, "_combine_", combine_target_clusters, "_UnitSpot")

# Real target and reference cluster, names from full cohort clusters, not subclusters of compartment.
cluster1 <- c("Muscle", "Stroma")
cluster2 <- c("Ducts", "Tumor_Cl_Int", "Tumor_Bas")

# distance cutoff 
distance_thres <- 10 # How many spot lengths

# color of reference clusters
color_tmp <- setNames(c("#E5E5AA","#B5B572","#5E5E39","#49392E", "#184928","#A0D5B5","#358E5B"), c("Fibroblasts", "myCAF", "iCAF", "ISG.RS", "Stellate1", "Stellate2", "Endothelial"))


## PATHS ----


# Input folder
input_path <- "/dir/to/integrated.rds"

# Output information
output_path_init <- "/output/dir/"
output_path <- paste0(output_path_init, cohort_name, "/")

# Create output folder
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)



## INGEST DATA ----


# Ingest a merged RDS Seurat objects
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
  dplyr::left_join(cluster_tab_stroma %>% dplyr::select(c("Coordinate","stromal_subcluster","stromal_subcluster_name")), by = "Coordinate") %>% 
  dplyr::left_join(data_seurat@meta.data %>% tibble::rownames_to_column("Coordinate") %>% dplyr::select(c("Coordinate", "Sample_Name"="user.Sample_Name")), by = "Coordinate")

#cluster_meta <- data_seurat@meta.data %>% tibble::rownames_to_column("Coordinate") %>% dplyr::left_join(cluster_tab, by = "Coordinate")

# Calculate global sample stats
sample_stat_df <- data_seurat@meta.data %>%
  tibble::rownames_to_column("Coordinate") %>% 
  dplyr::group_by(user.Sample_Name) %>%
  dplyr::summarise(spot_num = dplyr::n_distinct(Coordinate)) %>%
  dplyr::ungroup()


data_summary_df_cumu <- data.frame()
for(target_clusters_tmp in target_clusters){
  # keep only relevant spots
  cluster_tab_tmp <- cluster_tab[cluster_tab[[target_clusters_col]] %in% target_clusters_tmp | cluster_tab[[reference_clusters_col]] %in% reference_clusters,]
  #table(cluster_tab[cluster_tab[[target_clusters_col]] %in% target_clusters_tmp,]$fullcohort_cluster_name, cluster_tab[cluster_tab[[target_clusters_col]] %in% target_clusters_tmp,][[target_clusters_col]])
  target_clusters_tmp <- stringr::str_replace(target_clusters_tmp, "/", "_")
  # Count the number of spots selected by both signatures and their intersection
  spot_list1 <- cluster_tab_tmp %>%
    dplyr::filter(fullcohort_cluster_name %in% cluster1) %>%
    dplyr::pull(Coordinate)
  
  spot_list2 <- cluster_tab_tmp %>%
    dplyr::filter(fullcohort_cluster_name %in% cluster2) %>%
    dplyr::pull(Coordinate)
  
  spot_list3 <- intersect(spot_list1, spot_list2)
  
  print(c(length(spot_list1), length(spot_list2), length(spot_list3)))
  
  
  ## CALCULATE DISTANCE BETWEEN SPOTS ----
  
  
  spot_count_list <- list()
  spot_dist_list <- list()
  
  for(sample in unique(data_seurat@meta.data$user.Sample_Name)){
    
    # Select spots that are positive one or the other signature in this sample
    spots_select1 <- cluster_tab_tmp %>%
      dplyr::filter(fullcohort_cluster_name %in% cluster1 & Sample_Name == sample) %>%
      dplyr::pull(Coordinate) %>%
      unique
    
    # Select spots that are positive one or the other signature in this sample
    spots_select2 <- cluster_tab_tmp %>%
      dplyr::filter(fullcohort_cluster_name %in% cluster2 & Sample_Name == sample) %>%
      dplyr::pull(Coordinate) %>%
      unique
    
    # Number of spots corresponding to each of the signatures in the given sample
    spot_count_list[[sample]] <- c(sample, length(spots_select1), length(spots_select2))
    
    # Find spot coordinates for the 1st signature
    image1_df <- tibble::as_tibble(data_seurat@images[[sample]]@coordinates[c("imagerow", "imagecol")], rownames="Coordinate") %>%
      dplyr::filter(Coordinate %in% spots_select1) %>%
      dplyr::rename(X = imagecol, Y = imagerow)
    
    # Find spot coordinates for the 2nd signature
    image2_df <- tibble::as_tibble(data_seurat@images[[sample]]@coordinates[c("imagerow", "imagecol")], rownames="ClosestCoordinate") %>%
      dplyr::filter(ClosestCoordinate %in% spots_select2) %>%
      dplyr::rename(X = imagecol, Y = imagerow)
    
    # Find the closest barcode from one tibble in another
    spot_dist_list[[sample]] <- nearest_neighbor_df_my(df = image1_df, ref_df = image2_df, value="ClosestCoordinate")
  }
  
  # Bind together spot number data
  spot_count_df <- t(dplyr::bind_rows(spot_count_list)) %>%
    tibble::as_tibble() %>%
    dplyr::rename(Sample_Name = V1, Cluster1_Spots = V2, Cluster2_Spots = V3) %>%
    dplyr::inner_join(sample_stat_df, by=c("Sample_Name" = "user.Sample_Name")) %>%
    dplyr::mutate(Cluster1_Perc = round(as.numeric(Cluster1_Spots) / spot_num * 100, 1),
                  Cluster2_Perc = round(as.numeric(Cluster2_Spots) / spot_num * 100, 1))
  
  # Bind together spot distance data
  spot_dist_df <- dplyr::bind_rows(spot_dist_list) %>%
    dplyr::mutate(Sample_Name = stringr::str_remove(Coordinate, ":[:print:]*$"),
                  Distance_log2 = round(log2(Distance + 1), 2))
  
  # Extract spot radius
  radius_df <- data.frame(spot_radius = sapply(data_seurat@images, function(x) {x@scale.factors$spot_diameter_fullres/2})) %>% tibble::rownames_to_column("Sample_Name")
  radius_df$C2C_dist <- radius_df$spot_radius/27.5 * 100
  spot_dist_df <- spot_dist_df %>% dplyr::left_join(radius_df, by="Sample_Name") %>% dplyr::mutate(DistanceBySpot = Distance/C2C_dist)
  
  
  ## VISUALIZE DATA ----
  
  
  # Percent spots corresponding to the pair of signatures in each sample
  filename <- paste0(output_path, target_clusters_tmp, "_bar_cluster_spot_percent")
  create_bar_plot_my(spot_count_df %>%
                       dplyr::select(Sample_Name, Cluster1_Perc, Cluster2_Perc) %>%
                       df_wide2long_my(key="Cluster", val="SpotPercent"),
                     x_label="Sample_Name", y_label="SpotPercent", fill_label="Cluster", 
                     position="stack", filename=filename, labels=c("Sample_Name", "% spots", ""), reorder_x=TRUE)
  
  # Distances between two sets of signature positive spots in each sample
  filename <- paste0(output_path, target_clusters_tmp, "_box_sig_spot_distance")
  create_box_plot_my(spot_dist_df, x_label="Sample_Name", y_label="DistanceBySpot", fill_label="Sample_Name", 
                     filename=filename, labels=c("Sample Name", "Distance, spot", paste(paste0(cluster1, collapse = "+"), paste0(cluster2, collapse = "+"))), reorder_x=TRUE)
  
  
  # Spot status
  if("Status" %in% colnames(data_seurat@meta.data)){
    data_seurat@meta.data <- data_seurat@meta.data[,-which(colnames(data_seurat@meta.data) == "Status"), drop =F]
  }
  data_seurat@meta.data[which(rownames(data_seurat@meta.data) %in% spot_dist_df$Coordinate[spot_dist_df$DistanceBySpot <= distance_thres]), "Status"] <- reference_name
  data_seurat@meta.data[which(rownames(data_seurat@meta.data) %in% spot_dist_df$Coordinate[spot_dist_df$DistanceBySpot > distance_thres]), "Status"] <- paste0("Far_", reference_name)
  data_seurat@meta.data[which(rownames(data_seurat@meta.data) %in% spot_dist_df$ClosestCoordinate), "Status"] <- paste0("Proximal_", target_clusters_tmp)
  data_seurat@meta.data[which(rownames(data_seurat@meta.data) %in% setdiff(spot_list2, spot_dist_df$ClosestCoordinate)), "Status"] <- paste0("Other_", target_clusters_tmp)
  data_seurat@meta.data[which(is.na(data_seurat@meta.data$Status)), "Status"] <- "Exclude"
  
  
  # Spatial plot of selected spots
  cols <- c("blue","darkgrey", "red", "green", "lightgrey")
  names(cols) <- c(reference_name, paste0("Far_", reference_name),paste0("Proximal_", target_clusters_tmp), paste0("Other_", target_clusters_tmp), "Exclude" = "lightgrey")
  p <- spatial_dim_plot_my(data_seurat, group.by="Status",cols=cols, num_col = 7)
  
  # Write combo plot to file
  filename <- paste0(output_path, target_clusters_tmp, "_spatial_spots")
  write_plot2file_my(p, filename, num_row=7, num_col=7)
  
  
  ## SAVE DATA ----
  
  
  # Signature scores, calls, and distances joined
  data_summary_df <- cluster_tab_tmp %>%
    dplyr::inner_join(spot_dist_df, by=c("Coordinate", "Sample_Name")) %>% 
    dplyr::left_join(cluster_tab_tmp %>% dplyr::select("Coordinate", "ClosestCoordinateType"="epi_subcluster_name"), by = c("ClosestCoordinate" = "Coordinate"))
  
  data_summary_df[[reference_clusters_col]] <- factor(data_summary_df[[reference_clusters_col]], levels = reference_clusters)
  SpotFreq_plt <- data_summary_df %>%
    dplyr::filter(DistanceBySpot <= distance_thres) %>% 
    dplyr::mutate(DistanceBins = cut(DistanceBySpot, breaks = seq(min(DistanceBySpot)%/%1 *1, distance_thres, by = 1))) %>% 
    dplyr::group_by(DistanceBins, stromal_subcluster_name) %>% 
    dplyr::summarise(SpotCnt = dplyr::n()) %>% 
    dplyr::group_by(DistanceBins) %>% 
    dplyr::mutate(SpotFreq = SpotCnt/sum(SpotCnt)) %>% 
    ggplot(aes(x=DistanceBins, y= SpotFreq, fill = stromal_subcluster_name))+
    scale_fill_manual(values=color_tmp)+
    geom_bar(stat="identity")+theme_classic()+Seurat::RotatedAxis()+xlab(paste0("Number of spots from ", target_clusters_tmp))+ylab("Proportion of Stromal spots")
  png(paste0(output_path, target_clusters_tmp, "_SpotFreq_by_Distance.png"), width = 6, height = 5, res = 300, units = "in")
  print(SpotFreq_plt)
  dev.off()
  
  # Write signature and distance data to a file
  filename <- paste0(output_path, "table_", target_clusters_tmp, ".txt")
  readr::write_delim(data_summary_df, filename, delim="\t")
  
  data_summary_df$reference_cluster <- target_clusters_tmp
  data_summary_df_cumu <- rbind(data_summary_df_cumu, data_summary_df)
  
  # Add distance information to the original object
  meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
    dplyr::inner_join(data_summary_df %>%
                        dplyr::select("Coordinate", "ClosestCoordinate", "ClosestCoordinateType","Distance", "DistanceBySpot"), by="Coordinate")
  
  if(F){
    data_seurat_reference <- data_seurat
    # Update meta.data information
    data_seurat_reference@meta.data <- meta_df %>%
      tibble::column_to_rownames("Coordinate")
    
    
    # Spatial plot of distances
    #p <- spatial_feature_plot_my(data_seurat_reference, feature="DistanceBySpot", min.cutoff = NA, max.cutoff = NA, num_col = 7)?? 
    p <- Seurat::SpatialFeaturePlot(object = data_seurat_reference, 
                                    features = "DistanceBySpot",
                                    alpha = 1, 
                                    ncol = 7, 
                                    pt.size.factor = find_pt_size_factor_my(data_seurat), 
                                    min.cutoff = NA, 
                                    max.cutoff = NA, 
                                    crop = T, 
                                    slot = "data",
                                    images = NULL, 
                                    stroke = 0)
    # Write combo plot to file
    filename <- paste0(output_path, target_clusters_tmp, "_spatial_distance.png")
    png(filename, width = 12, height = 12, res = 300, units = "in")
    print(p)
    dev.off()
    
    
    # Subset and save object
    #data_seurat_reference <- subset(data_seurat_reference, cells=spot_dist_df$Coordinate)
    #filename <- paste0(output_path,"PDAC108_path14_fibs2", target_clusters_tmp, "_dist.rds")
    #saveRDS(data_seurat_reference, filename)
    #remove(data_seurat_reference)
    gc()
  }
  
  
}


# Write signature and distance data to a file
filename <- paste0(output_path, "Summary_dist_table.txt")
readr::write_delim(data_summary_df_cumu, filename, delim="\t")

# Visualize the proximal spots
data_summary_df_cumu <- data_summary_df_cumu %>%
  dplyr::filter(DistanceBySpot <= distance_thres) %>% 
  dplyr::mutate(DistanceBins = cut(DistanceBySpot, breaks = seq(min(DistanceBySpot)%/%1 *1, distance_thres, by = 1)))
for(distbin in sort(unique(data_summary_df_cumu$DistanceBins))){
  data_tmp <- data_summary_df_cumu[data_summary_df_cumu$DistanceBins == distbin,]
  bar_plt <- data_tmp %>%
    #dplyr::filter(DistanceBySpot <= distance_thres) %>% 
    #dplyr::mutate(DistanceBins = cut(DistanceBySpot, breaks = seq(min(DistanceBySpot)%/%1 *1, distance_thres, by = 1))) %>% 
    dplyr::group_by(reference_cluster, stromal_subcluster_name) %>% 
    dplyr::summarise(SpotCnt = dplyr::n()) %>% 
    dplyr::group_by(reference_cluster) %>% 
    dplyr::mutate(SpotFreq = SpotCnt/sum(SpotCnt)) %>% 
    dplyr::mutate(reference_cluster = factor(reference_cluster, levels=unique(data_summary_df_cumu$reference_cluster)),
                  stromal_subcluster_name = factor(stromal_subcluster_name, levels = names(color_tmp))) %>% 
    ggplot(aes(x=reference_cluster, y= SpotFreq, fill = stromal_subcluster_name))+
    scale_fill_manual(values=color_tmp)+
    geom_bar(stat="identity")+theme_classic()+Seurat::RotatedAxis()+xlab(paste0("Reference cluster"))+ylab(paste0("Proportion of stromal spots ", distbin, " away from reference"))
  
  png(paste0(output_path, "SpotFreq_by_Reference_", distbin, ".png"), width = 6, height = 5, res = 300, units = "in")
  print(bar_plt)
  dev.off()
}
for(ref in unique(data_summary_df_cumu$reference_cluster)){
  line_plt <- data_summary_df_cumu %>%
    dplyr::filter(DistanceBySpot <= distance_thres & reference_cluster == ref) %>% 
    dplyr::mutate(DistanceBins = cut(DistanceBySpot, breaks = seq(min(DistanceBySpot)%/%1 *1, distance_thres, by = 1))) %>% 
    dplyr::group_by(DistanceBins, stromal_subcluster_name) %>% 
    dplyr::summarise(SpotCnt = dplyr::n()) %>% 
    dplyr::group_by(DistanceBins) %>% 
    dplyr::mutate(SpotFreq = SpotCnt/sum(SpotCnt)) %>% 
    dplyr::mutate(stromal_subcluster_name = factor(stromal_subcluster_name, levels = names(color_tmp))) %>% 
    ggplot(aes(x=DistanceBins, y= SpotFreq, color = stromal_subcluster_name))+
    geom_line(aes(group=stromal_subcluster_name))+
    geom_point()+
    scale_color_manual(values=color_tmp)+theme_classic()+Seurat::RotatedAxis()+xlab(paste0("Number of spots from ", ref))+ylab("Proportion of Stromal spots")
  png(paste0(output_path, ref, "_line_SpotFreq_by_Distance.png"), width = 6, height = 5, res = 300, units = "in")
  print(line_plt)
  dev.off()
}

proportion_tab <- data_summary_df_cumu %>%
  dplyr::filter(DistanceBySpot <= distance_thres) %>% 
  dplyr::mutate(DistanceBins = cut(DistanceBySpot, breaks = seq(min(DistanceBySpot)%/%1 *1, distance_thres, by = 1))) %>% 
  dplyr::mutate(reference_cluster = factor(reference_cluster, levels=unique(data_summary_df_cumu$reference_cluster)),
                stromal_subcluster_name = factor(stromal_subcluster_name, levels = names(color_tmp))) %>% 
  dplyr::group_by(DistanceBins, stromal_subcluster_name, reference_cluster) %>% 
  dplyr::summarise(SpotCnt = dplyr::n()) %>% 
  dplyr::group_by(DistanceBins, reference_cluster) %>% 
  dplyr::mutate(SpotFreq = SpotCnt/sum(SpotCnt)) %>% 
  dplyr::mutate(stromal_subcluster_name = factor(stromal_subcluster_name, levels = names(color_tmp)))
filename <- paste0(output_path, "Summary_proportion_table.txt")
readr::write_delim(proportion_tab, filename, delim="\t")

bar_plt <- data_summary_df_cumu %>%
  dplyr::filter(DistanceBySpot <= distance_thres) %>% 
  dplyr::mutate(DistanceBins = cut(DistanceBySpot, breaks = seq(min(DistanceBySpot)%/%1 *1, distance_thres, by = 1))) %>% 
  dplyr::mutate(reference_cluster = factor(reference_cluster, levels=unique(data_summary_df_cumu$reference_cluster)),
                stromal_subcluster_name = factor(stromal_subcluster_name, levels = names(color_tmp))) %>% 
  dplyr::group_by(DistanceBins, stromal_subcluster_name, reference_cluster) %>% 
  dplyr::summarise(SpotCnt = dplyr::n()) %>% 
  dplyr::group_by(DistanceBins, reference_cluster) %>% 
  dplyr::mutate(SpotFreq = SpotCnt/sum(SpotCnt)) %>% 
  dplyr::mutate(stromal_subcluster_name = factor(stromal_subcluster_name, levels = names(color_tmp))) %>% 
  ggplot(aes(x=DistanceBins, y= SpotFreq, fill = stromal_subcluster_name))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=color_tmp)+theme_classic()+Seurat::RotatedAxis()+xlab(paste0("Number of spots from tumor spots"))+ylab("Proportion of Stromal spots")+
  facet_wrap(~reference_cluster)
line_plt <- data_summary_df_cumu %>%
  dplyr::filter(DistanceBySpot <= distance_thres) %>% 
  dplyr::mutate(DistanceBins = cut(DistanceBySpot, breaks = seq(min(DistanceBySpot)%/%1 *1, distance_thres, by = 1))) %>% 
  dplyr::mutate(reference_cluster = factor(reference_cluster, levels=unique(data_summary_df_cumu$reference_cluster)),
                stromal_subcluster_name = factor(stromal_subcluster_name, levels = names(color_tmp))) %>% 
  dplyr::group_by(DistanceBins, stromal_subcluster_name, reference_cluster) %>% 
  dplyr::summarise(SpotCnt = dplyr::n()) %>% 
  dplyr::group_by(DistanceBins, reference_cluster) %>% 
  dplyr::mutate(SpotFreq = SpotCnt/sum(SpotCnt)) %>% 
  dplyr::mutate(stromal_subcluster_name = factor(stromal_subcluster_name, levels = names(color_tmp))) %>% 
  ggplot(aes(x=DistanceBins, y= SpotFreq, color = stromal_subcluster_name))+
  geom_line(aes(group=stromal_subcluster_name))+
  geom_point()+
  scale_color_manual(values=color_tmp)+theme_classic()+Seurat::RotatedAxis()+xlab(paste0("Number of spots from tumor spots"))+ylab("Proportion of Stromal spots")+
  facet_wrap(~reference_cluster)
area_plt <- data_summary_df_cumu %>%
  dplyr::filter(DistanceBySpot <= distance_thres) %>% 
  dplyr::mutate(DistanceBins = cut(DistanceBySpot, breaks = seq(min(DistanceBySpot)%/%1 *1, distance_thres, by = 1))) %>% 
  dplyr::mutate(reference_cluster = factor(reference_cluster, levels=unique(data_summary_df_cumu$reference_cluster)),
                stromal_subcluster_name = factor(stromal_subcluster_name, levels = names(color_tmp))) %>% 
  dplyr::group_by(DistanceBins, stromal_subcluster_name, reference_cluster) %>% 
  dplyr::summarise(SpotCnt = dplyr::n()) %>% 
  dplyr::group_by(DistanceBins, reference_cluster) %>% 
  dplyr::mutate(SpotFreq = SpotCnt/sum(SpotCnt)) %>% 
  dplyr::mutate(stromal_subcluster_name = factor(stromal_subcluster_name, levels = names(color_tmp))) %>% 
  ggplot(aes(group =stromal_subcluster_name))+
  geom_area(aes(x=DistanceBins,y=SpotFreq, fill = stromal_subcluster_name))+
  scale_fill_manual(values=color_tmp)+theme_classic()+Seurat::RotatedAxis()+xlab(paste0("Number of spots from tumor spots"))+ylab("Proportion of Stromal spots")+
  facet_wrap(~reference_cluster)

png(paste0(output_path, "Summary_bar_SpotFreq_by_Distance.png"), width = 10, height = 8, res = 300, units = "in")
print(bar_plt)
dev.off()

png(paste0(output_path, "Summary_line_SpotFreq_by_Distance.png"), width = 10, height = 8, res = 300, units = "in")
print(line_plt)
dev.off()

png(paste0(output_path, "Summary_area_SpotFreq_by_Distance.png"), width = 10, height = 8, res = 300, units = "in")
print(area_plt)
dev.off()
