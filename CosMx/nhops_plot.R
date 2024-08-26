library(tidyverse)
library(monocle3)

home_path <- "/out/cosmx_TAP/"
int_path <- "/out/cosmx_TAP/intermediate/" 

tumor <- readRDS(file = paste0(int_path, "all_slides_tumor_integrated.RDS"))

pdf(paste0(home_path,"man_plots/F3N_tumor_neighbors_average_heatmaps.pdf"), h=3.5, w=9)
for (i in c(1:9)) {
    comp_mat <- read.csv(paste0(int_path, "TME_analysis/tme_subtype/composition_full_qc_20_", as.character(i), "hop" ,".csv"))
    meta <- read.csv(paste0(int_path, "TME_analysis/tme_subtype/composition_meta_", as.character(i), "hop" , ".csv"))
    ##all_edges <- readRDS(file = paste0(int_path, "all_edges.RDS"))
    
    comp_sub <- merge(comp_mat, meta, by = c("Subgraph.ID", "FOV"))

    comp_sub$slide <- stringr::str_split_fixed(comp_sub$fov, "_", 2)[,1]
    comp_sub$slide <- gsub('slide', "S", comp_sub$slide)
    comp_sub$center_cell_id <- paste0(comp_sub$center_cell, "_", comp_sub$slide)
    row.names(comp_sub) <- comp_sub$center_cell_id
    comp_sub$subgraph <- paste0(comp_sub$fov, "_", comp_sub$Subgraph.ID) 
    
    comp_sub_long <- tidyr::separate_rows(comp_sub, cell_members, sep = ",")
    comp_sub_long$cell_members <- gsub("\\[", "", comp_sub_long$cell_members)
    comp_sub_long$cell_members <- gsub("\\]", "", comp_sub_long$cell_members)
    comp_sub_long$cell_members <- gsub(" ", "", comp_sub_long$cell_members)
    comp_sub_long$cell_members <- gsub("'", "", comp_sub_long$cell_members)
    
    comp_sub_long$cell_members <- paste0(comp_sub_long$cell_members, "_", comp_sub_long$slide)
    comp_sub$center_cell_type <- 
        pData(tumor)[comp_sub$center_cell_id,]$comb_labels_clusts
    
    for_plot <- comp_sub[comp_sub$center_cell_type %in%  c("Ductal (6)",
                                                           "Tumor Classical (1,2,3,5,7,8,11,14)",
                                                           "Tumor Classical / Proliferative (9)",
                                                           "Tumor Mixed (10)",
                                                           "Tumor Basal (4,15)",
                                                           "Tumor Mixed / Fibro adjacent (12)"),]
    
    for_plot_count <- for_plot %>%
        group_by(center_cell_type) %>% 
        summarize(Acinar = mean(Acinar),
                  B.cell.Plasma = mean(B.cell.Plasma),
                  Basal = mean(Basal),
                  Classical = mean(Classical),
                  Dendritic = mean(Dendritic),
                  Ductal = mean(Ductal),
                  Endocrine = mean(Endocrine),
                  Endothelial = mean(Endothelial),
                  Excluded = mean(Excluded),
                  Fib.high = mean(Fib.high),
                  Fibroblast = mean(Fibroblast),
                  Low.quality = mean(Low.quality),
                  Macrophage = mean(Macrophage),
                  Mast = mean(Mast),
                  Mixed = mean(Mixed),
                  Other.Tumor = mean(Other.Tumor),
                  Proliferative = mean(Proliferative),
                  Stellate = mean(Stellate),
                  T.cell = mean(T.cell))
    
    for_plot_count <- as.data.frame(for_plot_count)
    row.names(for_plot_count) <- for_plot_count$center_cell_type
    for_plot_count$center_cell_type <- NULL
    
    for_plot_count <- for_plot_count[c("Ductal (6)",
                                       "Tumor Classical (1,2,3,5,7,8,11,14)",
                                       "Tumor Classical / Proliferative (9)",
                                       "Tumor Mixed (10)",
                                       "Tumor Basal (4,15)",
                                       "Tumor Mixed / Fibro adjacent (12)"),
                                     c("Ductal",
                                       "Classical",
                                       "Proliferative",
                                       "Mixed",
                                       "Basal",
                                       "Fib.high",
                                       "Acinar",
                                       "B.cell.Plasma",
                                       "Dendritic",
                                       "Endocrine",
                                       "Endothelial",
                                       "Fibroblast",
                                       "Macrophage",
                                       "Mast",
                                       "Stellate",
                                       "T.cell")]
    
    
    ##png(paste0(home_path, "man_plots/F3S2G_tumor_neighbors_average_heatmap_", as.character(i), "hop" ,".png"), width = 9, height = 3.5, res = 300, units = "in")
    print(pheatmap::pheatmap(for_plot_count, cluster_rows = F, cluster_cols = F, 
                             color = RColorBrewer::brewer.pal(9, "Blues"), 
                             display_numbers = round(for_plot_count * 100, digits = 1), 
                             number_color = "darkgrey"))
    ##dev.off()
}
dev.off()

pdf(paste0(home_path,"man_plots/F3N_tumor_neighbors_heatmaps.pdf"), h=3.5, w=9)
for (i in c(1:9)) {
    comp_mat <- read.csv(paste0(int_path, "TME_analysis/tme_subtype/composition_full_qc_20_", as.character(i), "hop" ,".csv"))
    meta <- read.csv(paste0(int_path, "TME_analysis/tme_subtype/composition_meta_", as.character(i), "hop" , ".csv"))
    ##all_edges <- readRDS(file = paste0(int_path, "all_edges.RDS"))
    
    comp_sub <- merge(comp_mat, meta, by = c("Subgraph.ID", "FOV"))

    comp_sub$slide <- stringr::str_split_fixed(comp_sub$fov, "_", 2)[,1]
    comp_sub$slide <- gsub('slide', "S", comp_sub$slide)
    comp_sub$center_cell_id <- paste0(comp_sub$center_cell, "_", comp_sub$slide)
    row.names(comp_sub) <- comp_sub$center_cell_id
    comp_sub$subgraph <- paste0(comp_sub$fov, "_", comp_sub$Subgraph.ID)
    
    pData(tumor)$comb_labels_clusts <- as.character(pData(tumor)$comb_labels_clusts)
    
    comp_sub$center_cell_type <- pData(tumor)[comp_sub$center_cell_id,]$comb_labels_clusts
    for_plot <- comp_sub[comp_sub$center_cell_type %in%  c("Ductal (6)",
                                                           "Tumor Classical (1,2,3,5,7,8,11,14)",
                                                           "Tumor Classical / Proliferative (9)",
                                                           "Tumor Mixed (10)",
                                                           "Tumor Basal (4,15)",
                                                           "Tumor Mixed / Fibro adjacent (12)"),] 
    
    for_plot_count <- for_plot %>% group_by(center_cell_type) %>% 
        summarize(Acinar = sum(Acinar > 0)/n(),
                  B.cell.Plasma = sum(B.cell.Plasma > 0)/n(),
                  Basal = sum(Basal > 0)/n(),
                  Classical = sum(Classical > 0)/n(),
                  Dendritic = sum(Dendritic > 0)/n(),
                  Ductal = sum(Ductal > 0)/n(),
                  Endocrine = sum(Endocrine > 0)/n(),
                  Endothelial = sum(Endothelial > 0)/n(),
                  Excluded = sum(Excluded > 0)/n(),
                  Fib.high = sum(Fib.high > 0)/n(),
                  Fibroblast = sum(Fibroblast > 0)/n(),
                  Low.quality = sum(Low.quality > 0)/n(),
                  Macrophage = sum(Macrophage > 0)/n(),
                  Mast = sum(Mast > 0)/n(),
                  Mixed = sum(Mixed > 0)/n(),
                  Other.Tumor = sum(Other.Tumor > 0)/n(),
                  Proliferative = sum(Proliferative > 0)/n(),
                  Stellate = sum(Stellate > 0)/n(),
                  T.cell = sum(T.cell > 0)/n())
    
    for_plot_count <- as.data.frame(for_plot_count)
    row.names(for_plot_count) <- for_plot_count$center_cell_type
    for_plot_count$center_cell_type <- NULL
    
    for_plot_count <- for_plot_count[c("Ductal (6)",
                                       "Tumor Classical (1,2,3,5,7,8,11,14)",
                                       "Tumor Classical / Proliferative (9)",
                                       "Tumor Mixed (10)",
                                       "Tumor Basal (4,15)",
                                       "Tumor Mixed / Fibro adjacent (12)"),
                                     c("Ductal",
                                       "Classical",
                                       "Proliferative",
                                       "Mixed",
                                       "Basal",
                                       "Fib.high",
                                       "Acinar",
                                       "B.cell.Plasma",
                                       "Dendritic",
                                       "Endocrine",
                                       "Endothelial",
                                       "Fibroblast",
                                       "Macrophage",
                                       "Mast",
                                       "Stellate",
                                       "T.cell")]
    
    colmins <- apply(for_plot_count, 2, min)
    colmax <- apply(for_plot_count, 2, max)
    colrange <- colmax - colmins
    for_plot_count_sc <- t((t(for_plot_count) - colmins)/colrange)
    
    ##png(paste0(home_path, "man_plots/F3N_tumor_neighbors_heatmap_", as.character(i), "hop" ,".png"), width = 9, height = 3.5, res = 300, units = "in")
    print(pheatmap::pheatmap(for_plot_count_sc,
                             cluster_rows = F,
                             cluster_cols = F,
                             main = paste0('Tumor neighbors heatmap ', as.character(i), "hop"),
                             color = RColorBrewer::brewer.pal(9, "Blues"), 
                             display_numbers = round(for_plot_count * 100, digits = 1), 
                             number_color = "darkgrey"))
    ##dev.off()
}
dev.off()
