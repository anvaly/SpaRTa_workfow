suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(monocle3)
  library(clusterProfiler)
  library(piano)
  library(dplyr)
  library(SeuratWrappers)
  library(igraph)
  library(data.table)
})

stash_data <- "/path/to/data/cosmx/"
stash_base <- "/path/to/data/"
home_path <- "/out/cosmx_TAP/"
int_path <- "/out/cosmx_TAP/intermediate/"

set.seed(2023)

patient_ROI_labels <- list("S1_P1" = "C2_D6",
                           "S1_P2" = "C2_D1",
                           "S2_P1" = "C2_D10",
                           "S2_P2" = "C2_D8",
                           "S3_P1" = "C3_D4_ROI4",
                           "S3_P2" = "C3_D4_ROI2",
                           "S4_P1" = "C3_D12_ROI2",
                           "S4_P2" = "C3_D12_ROI3",
                           "S5_P1" = "C3_D9_ROI3",
                           "S5b2_P2" = "C3_D9_ROI1a",
                           "S5b3_P2" = "C3_D9_ROI1b",
                           "S6_P1" = "C3_D5",
                           "S6_P2" = "C3_D6")

patient_labels <- list("S1_P1" = "C2_D6",
                       "S1_P2" = "C2_D1",
                       "S2_P1" = "C2_D10",
                       "S2_P2" = "C2_D8",
                       "S3_P1" = "C3_D4",
                       "S3_P2" = "C3_D4",
                       "S4_P1" = "C3_D12",
                       "S4_P2" = "C3_D12",
                       "S5_P1" = "C3_D9",
                       "S5b2_P2" = "C3_D9",
                       "S5b3_P2" = "C3_D9",
                       "S6_P1" = "C3_D5",
                       "S6_P2" = "C3_D6")

all_cols <- list('Acinar 1' = "#056DB5",
                 'Acinar 2' = "#056DB5",
                 'Acinar' = "#056DB5",
                 'Acinar (16,21)' = "#056DB5",
                 'Endocrine' = "#153C65",
                 'Endocrine (23)' = "#153C65",
                 'Ductal' = "#4CB0E1",
                 'Ductal (12)' = "#4CB0E1",
                 'Tumor' = "#9A2626",
                 'Tumor (4,6,8,9,10,15,17,19,26)' = "#9A2626",
                 'Tumor 1' = "#9A2626",
                 'Tumor 2' = "#9A2626",
                 'Tumor 3' = "#9A2626",
                 'Tumor 4' = "#9A2626",
                 'Tumor 5' = "#9A2626",
                 'Tumor 6' = "#9A2626",
                 'Tumor 7' = "#9A2626",
                 'Tumor 8' = "#9A2626",
                 'Tumor 9' = "#9A2626",
                 'Fibroblast' = "#E5E5AA",
                 'Fibroblast (1,3,20)' = "#E5E5AA",
                 'Fibroblast 1' = "#E5E5AA",
                 'Fibroblast 2' = "#E5E5AA",
                 'Fibroblast 3' = "#E5E5AA",
                 'Stellate' = "#A0D5B5",
                 'Stellate (18)' = "#A0D5B5",
                 'Endothelial' = "#358E5B",
                 'Endothelial (14,25)' = "#358E5B",
                 'Endothelial 1' = "#358E5B",
                 'Endothelial 2' = "#358E5B",
                 'B cell/Plasma' = "#FFCE07",
                 'B cell/Plasma (13,24)' = "#FFCE07",
                 'B cell/Plasma 1' = "#FFCE07",
                 'B cell/Plasma 2' = "#FFCE07",
                 'T cell' = "#D66100",
                 'T cell (5)' = "#D66100",
                 'Macrophage' = "#FF9F2C",
                 'Macrophage (2,11,28)' = "#FF9F2C",
                 'Macrophage 1' = "#FF9F2C",
                 'Macrophage 2' = "#FF9F2C",
                 'Macrophage 3' = "#FF9F2C",
                 'Dendritic' = "#9B5E2C",
                 'Dendritic (27)' = "#9B5E2C",
                 'Mast' = "#B1571B",
                 'Mast (22)' = "#B1571B",
                 'Low quality' = "lightgray",
                 'Low quality (7)' = "lightgray")

clust_cols <- list('16' = "#056DB5",
                   '21' = "#056DB5",
                   '23' = "#153C65",
                   '12' = "#4CB0E1",
                   '4' = "#9A2626",
                   '6' = "#9A2626",
                   '8' = "#9A2626",
                   '10' = "#9A2626",
                   '15' = "#9A2626",
                   '17' = "#9A2626",
                   '19' = "#9A2626",
                   '9' = "#9A2626",
                   '26' = "#9A2626",
                   '1' = "#E5E5AA",
                   '3' = "#E5E5AA",
                   '20' = "#E5E5AA",
                   '18' = "#A0D5B5",
                   '14' = "#358E5B",
                   '25' = "#358E5B",
                   '13' = "#FFCE07",
                   '24' = "#FFCE07",
                   '5' = "#D66100",
                   '2' = "#FF9F2C",
                   '11' = "#FF9F2C",
                   '28' = "#FF9F2C",
                   '27' = "#9B5E2C",
                   '22' = "#B1571B",
                   '7' = "lightgray")

tumor_clust_cols <- list('6' = "#4CB0E1",
                         '1' = "#F16666",
                         '2' = "#F16666",
                         '3' = "#F16666",
                         '5' = "#F16666",
                         '7' = "#F16666",
                         '8' = "#F16666",
                         '14' = "#F16666",
                         '9' = "#9A2626",
                         '11' = "#F16666",
                         '10' = "#CA3FAE",
                         '12' = "#2C9FAF",
                         '4' = "#681D59",
                         '15' = "#681D59",
                         '13' = "lightgrey")

tumor_cols <- list('Ductal (6)' = "#4CB0E1",
                   'Tumor Classical (1,2,3,5,7,8,11,14)' = "#F16666",
                   'Tumor Classical / Proliferative (9)' = "#9A2626",
                   'Tumor Mixed (10)' = '#CA3FAE',
                   'Tumor Mixed / Fibro adjacent (12)' = "#2C9FAF",
                   'Tumor Basal (4,15)' = "#681D59",
                   'Tumor Other (13)' = "lightgrey")

stroma_clust_cols <- list('15' = "#E5E5AA",
                          '1' = "#E5E5AA",
                          '2' = "#E5E5AA",
                          '5' = "#E5E5AA",
                          '14' = "#E5E5AA",
                          '13' = "#E5E5AA",
                          '10' = "#E5E5AA",
                          '7' = "#E5E5AA",
                          '9' = "#E5E5AA",
                          '4' = "#E5E5AA",
                          '6' = "#184928",
                          '8' = "#A0D5B5",
                          '3' = "#358E5B",
                          '12' = "#358E5B",
                          '11' = "#49392E",
                          '16' = "#939393")

stroma_cols <- list('Fibroblast (1,2,4,5,7,9,10,13,14,15)' = "#E5E5AA",
                    'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = "#E5E5AA",
                    'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = "#E5E5AA",
                    'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = "#E5E5AA",
                    'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = "#E5E5AA",
                    'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = "#E5E5AA",
                    'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = "#E5E5AA",
                    'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = "#E5E5AA",
                    'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = "#E5E5AA",
                    'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = "#E5E5AA",
                    'Stellate 1 (6)' = "#184928",
                    'Stellate 2 (8)' = "#A0D5B5",
                    'Endothelial (3,12)' = "#358E5B",
                    'Endothelial (3,12)' = "#358E5B",
                    'Fibroblasts/IFNa pathway/ resistance (11)' = "#49392E",
                    'Nerve (16)' = "#939393")

imm_clust_cols <- list('7' = "#FEDC4B",
                       '9' = "#FEDC4B",
                       '17' = "#FEDC4B",
                       '16' = "#FFCE07",
                       '23' = "gray",
                       '1' = "#FF9F2C",
                       '2' = "#FF9F2C",
                       '5' = "#FF9F2C",
                       '8' = "#FF9F2C",
                       '10' = "#FF9F2C",
                       '11' = "#FF9F2C",
                       '13' = "#FF9F2C",
                       '14' = "#FF9F2C",
                       '18' = "#FF9F2C",
                       '20' = "#FF9F2C",
                       '21' = "#FF9F2C",
                       '22' = "#FF9F2C",
                       '3' = "#D66100",
                       '4' = "#D66100",
                       '6' = "#D66100",
                       '15' = "#D66100",
                       '19' = "#D66100",
                       '12' = "#E0B400")

imm_cols <- list('Plasma (7,9,17)' = "#FEDC4B",
                 'Plasma (7,9,17)' = "#FEDC4B",
                 'Plasma (7,9,17)' = "#FEDC4B",
                 'B cell (16)' = "#FFCE07",
                 '? (23)' = "gray",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'Macrophage (1,2,5,8,10,11,13,14,18,20,21,22)' = "#FF9F2C",
                 'T cell (3,4,6,15,19)' = "#D66100",
                 'T cell (3,4,6,15,19)' = "#D66100",
                 'T cell (3,4,6,15,19)' = "#D66100",
                 'T cell (3,4,6,15,19)' = "#D66100",
                 'T cell (3,4,6,15,19)' = "#D66100",
                 'Mast (12)' = "#E0B400")


type_sigs <- c("PDAC.P19.Acinar",
               "PDAC.P19.Endocrine",
               "PDAC.P19.Ductal_1",
               "PDAC.collisson.classical",
               "PDAC.moffitt.basal",
               "PDAC.P19.Fibroblast",
               "PDAC.P19.Stellate",
               "PDAC.P19.Endothelial",
               "PDAC.P19.Macrophage",
               "PDAC.P19.Tcell",
               "PDAC.P19.Bcell",
               "PDAC.U.Nervous",
               "PDAC.P19.Ductal_2",
               "Seurat.Proliferation")


slide1 <- readRDS(paste0(int_path, "slide1.RDS"))
slide2 <- readRDS(paste0(int_path, "slide2.RDS"))
slide3 <- readRDS(paste0(int_path, "slide3.RDS"))
slide4 <- readRDS(paste0(int_path, "slide4.RDS"))
slide5 <- readRDS(paste0(int_path, "slide5.RDS"))
slide5@images$P2 <- NULL
slide5_P1 <- slide5[,slide5@meta.data$CenterX_global_px > 0]
slide5b2 <- readRDS(paste0(int_path, "slide5b2.RDS"))
slide5b2@images$P1 <- NULL
slide5b2_P2 <- slide5b2[,slide5b2@meta.data$CenterX_global_px < 0]
slide5b3 <- readRDS(paste0(int_path, "slide5b3.RDS"))
slide5b3@images$P1 <- NULL
slide5b3_P2 <- slide5b3[,slide5b3@meta.data$CenterX_global_px < 0]
slide6 <- readRDS(paste0(int_path, "slide6.RDS"))

slide_list <- list(slide1, slide2, slide3, slide4, slide5_P1, slide5b2_P2, slide5b3_P2, slide6)
names(slide_list) <- c(1,2,3,4,5,'5b2', '5b3', 6)







pData(all_slides)$tumor_comb_labels <- NA
pData(tumor)$comb_labels_clusts <- as.character(pData(tumor)$comb_labels_clusts)
pData(all_slides)$tumor_comb_labels <- pData(tumor)[row.names(pData(all_slides)),]$comb_labels_clusts
pData(all_slides)$tumor_comb_labels <-  as.character(pData(all_slides)$tumor_comb_labels)

pData(all_slides)$tme_subtype <- ifelse(is.na(pData(all_slides)$tumor_comb_labels),
                                         as.character(pData(all_slides)$comb_labels),
                                         as.character(pData(all_slides)$tumor_comb_labels))

pData(all_slides)$tme_subtype[pData(all_slides)$tme_subtype == "Tumor Basal (4,15)"] <- "Basal"
pData(all_slides)$tme_subtype[pData(all_slides)$tme_subtype %in% c("Tumor Classical (1,2,3,5,7,8,11,14)")] <- "Classical"
pData(all_slides)$tme_subtype[pData(all_slides)$tme_subtype %in% c("Tumor Classical / Proliferative (9)")] <- "Proliferative"
pData(all_slides)$tme_subtype[pData(all_slides)$tme_subtype %in% c("Tumor Mixed (10)")] <- "Mixed"
pData(all_slides)$tme_subtype[pData(all_slides)$tme_subtype %in% c("Tumor Mixed / Fibro adjacent (12)")] <- "Fib-high"
pData(all_slides)$tme_subtype[pData(all_slides)$tme_subtype == "Tumor Other (13)"] <- "Other Tumor"
pData(all_slides)$tme_subtype[pData(all_slides)$tme_subtype == "Ductal (6)"] <- "Ductal"


## One set with tumor subtype
for(slide_num in c(1,2,3,4,5,'5b2', '5b3', 6)) {
  slide_list[[slide_num]]@meta.data[['tme_subtype']] <- pData(all_slides)[paste0(row.names(slide_list[[slide_num]]@meta.data), 
                                                               "_S", slide_num),][['tme_subtype']]
  temp <- as.data.frame(slide_list[[slide_num]]@images[[paste0("Slide", slide_num)]]@boundaries$centroids@coords)
  temp$CELL_ID <- slide_list[[slide_num]]@images[[paste0("Slide", slide_num)]]@boundaries$centroids@cells
  names(temp) <- c("X", "Y", "CELL_ID")
  temp <- temp[,c("CELL_ID", "X", "Y")]
  temp$X <- round(temp$X + abs(min(temp$X)) + 1, 0)
  temp$Y <- round(temp$Y + abs(min(temp$Y)) + 1, 0)
  temp$fov <- stringr::str_split_fixed(temp$CELL_ID, "_", 2)[,2]
  for(fov in unique(temp$fov)) {
    write.csv(temp[temp$fov == fov, c(1:3)], quote = FALSE, row.names = FALSE, 
              file = paste0(home_path, "data/intermediate/TME_analysis/tme_subtype/slide", 
                            slide_num, "_FOV",fov, ".cell_data.csv"))
  }
  temp <- data.frame(CELL_ID = row.names(slide_list[[slide_num]]@meta.data),
                     CELL_TYPE = slide_list[[slide_num]]@meta.data$tme_subtype)
  temp$CELL_TYPE[is.na(temp$CELL_TYPE)] <- "Excluded"
  temp$fov <- stringr::str_split_fixed(temp$CELL_ID, "_", 2)[,2]
  for(fov in unique(temp$fov)) {
    write.csv(temp[temp$fov == fov, c(1:2)], quote = FALSE, row.names = FALSE, 
              file = paste0(home_path, "/data/intermediate/TME_analysis/tme_subtype/slide", 
                            slide_num, "_FOV",fov, ".cell_types.csv"))
  }
}


gen_edge_list <- function(field, seur_obj, max_dist = 150) {
  temp <- as.data.table(seur_obj@images[[field]]$centroids@coords)
  names(temp) <- c("sdimx", "sdimy")
  temp$cell_ID <- seur_obj@images[[field]]$centroids@cells
  dtest <- Giotto:::create_delaunayNetwork_geometry(temp)
  
  dtest$delaunay_network_DT$dist <- sqrt((dtest$delaunay_network_DT$sdimx_begin - dtest$delaunay_network_DT$sdimx_end)^2 + 
                                           (dtest$delaunay_network_DT$sdimy_begin - dtest$delaunay_network_DT$sdimy_end)^2)
  
  edges <- dtest$delaunay_network_DT[dtest$delaunay_network_DT$dist < max_dist,]
  return(edges)
}

gen_mean_neighbor_mat <- function(edge_list, seur_obj) {
  gr <- graph_from_edgelist(as.matrix(edge_list[,c("from", "to")]), 
                            directed = FALSE)
  gr_adj <- igraph::as_adjacency_matrix(gr)
  expr_fov <- seur_obj@assays$Nanostring@counts[,row.names(gr_adj)]
  neighbor_exp <- expr_fov %*% gr_adj
  neighbor_mean <- t(t(neighbor_exp)/rowSums(gr_adj))
  return(neighbor_mean)
}

slide1 <- readRDS(paste0(int_path, "slide1.RDS"))
slide2 <- readRDS(paste0(int_path, "slide2.RDS"))
slide3 <- readRDS(paste0(int_path, "slide3.RDS"))
slide4 <- readRDS(paste0(int_path, "slide4.RDS"))
slide5 <- readRDS(paste0(int_path, "slide5.RDS"))
slide5@images$P2 <- NULL
slide5_P1 <- slide5[,slide5@meta.data$CenterX_global_px > 0]
slide5b2 <- readRDS(paste0(int_path, "slide5b2.RDS"))
slide5b2@images$P1 <- NULL
slide5b2_P2 <- slide5b2[,slide5b2@meta.data$CenterX_global_px < 0]
slide5b3 <- readRDS(paste0(int_path, "slide5b3.RDS"))
slide5b3@images$P1 <- NULL
slide5b3_P2 <- slide5b3[,slide5b3@meta.data$CenterX_global_px < 0]
slide6 <- readRDS(paste0(int_path, "slide6.RDS"))
rm(slide5, slide5b2, slide5b3)
gc()

## Generate edges
edges1 <- gen_edge_list("Slide1", slide1)  
edges2 <- gen_edge_list("Slide2", slide2)  
edges3 <- gen_edge_list("Slide3", slide3)  
edges4 <- gen_edge_list("Slide4", slide4)  
edges5 <- gen_edge_list("Slide5", slide5_P1)
edges5b2 <- gen_edge_list("Slide5b2", slide5b2_P2)
edges5b3 <- gen_edge_list("Slide5b3", slide5b3_P2)
edges6 <- gen_edge_list("Slide6", slide6)  

gr <- graph_from_edgelist(as.matrix(edges1[,c("from", "to")]), directed = FALSE)
gr_adj <- igraph::as_adjacency_matrix(gr)
neighbor_count1 <- rowSums(gr_adj)
names(neighbor_count1) <- paste0(names(neighbor_count1), "_S1")

gr <- graph_from_edgelist(as.matrix(edges2[,c("from", "to")]), directed = FALSE)
gr_adj <- igraph::as_adjacency_matrix(gr)
neighbor_count2 <- rowSums(gr_adj)
names(neighbor_count2) <- paste0(names(neighbor_count2), "_S2")

gr <- graph_from_edgelist(as.matrix(edges3[,c("from", "to")]), directed = FALSE)
gr_adj <- igraph::as_adjacency_matrix(gr)
neighbor_count3<- rowSums(gr_adj)
names(neighbor_count3) <- paste0(names(neighbor_count3), "_S3")

gr <- graph_from_edgelist(as.matrix(edges4[,c("from", "to")]), directed = FALSE)
gr_adj <- igraph::as_adjacency_matrix(gr)
neighbor_count4 <- rowSums(gr_adj)
names(neighbor_count4) <- paste0(names(neighbor_count4), "_S4")

gr <- graph_from_edgelist(as.matrix(edges5[,c("from", "to")]), directed = FALSE)
gr_adj <- igraph::as_adjacency_matrix(gr)
neighbor_count5 <- rowSums(gr_adj)
names(neighbor_count5) <- paste0(names(neighbor_count5), "_S5")

gr <- graph_from_edgelist(as.matrix(edges5b2[,c("from", "to")]), directed = FALSE)
gr_adj <- igraph::as_adjacency_matrix(gr)
neighbor_count5b2 <- rowSums(gr_adj)
names(neighbor_count5b2) <- paste0(names(neighbor_count5b2), "_S5b2")

gr <- graph_from_edgelist(as.matrix(edges5b3[,c("from", "to")]), directed = FALSE)
gr_adj <- igraph::as_adjacency_matrix(gr)
neighbor_count5b3 <- rowSums(gr_adj)
names(neighbor_count5b3) <- paste0(names(neighbor_count5b3), "_S5b3")

gr <- graph_from_edgelist(as.matrix(edges6[,c("from", "to")]), directed = FALSE)
gr_adj <- igraph::as_adjacency_matrix(gr)
neighbor_count6 <- rowSums(gr_adj)
names(neighbor_count6) <- paste0(names(neighbor_count6), "_S6")

neighbor_count <- c(neighbor_count1, neighbor_count2, neighbor_count3, 
                    neighbor_count4, neighbor_count5, neighbor_count5b2, 
                    neighbor_count5b3, neighbor_count6)
pData(all_slides)$neighbor_count <- neighbor_count[row.names(pData(all_slides))]
pData(all_slides)$neighbor_count[is.na(pData(all_slides)$neighbor_count)] <- 0

neighbor_mean1 <- gen_mean_neighbor_mat(edges1, slide1)
colnames(neighbor_mean1) <- paste0(colnames(neighbor_mean1), "_S1")
neighbor_mean2 <- gen_mean_neighbor_mat(edges2, slide2)
colnames(neighbor_mean2) <- paste0(colnames(neighbor_mean2), "_S2")
neighbor_mean3 <- gen_mean_neighbor_mat(edges3, slide3)
colnames(neighbor_mean3) <- paste0(colnames(neighbor_mean3), "_S3")
neighbor_mean4 <- gen_mean_neighbor_mat(edges4, slide4)
colnames(neighbor_mean4) <- paste0(colnames(neighbor_mean4), "_S4")
neighbor_mean5 <- gen_mean_neighbor_mat(edges5, slide5_P1)
colnames(neighbor_mean5) <- paste0(colnames(neighbor_mean5), "_S5")
neighbor_mean5b2 <- gen_mean_neighbor_mat(edges5b2, slide5b2_P2)
colnames(neighbor_mean5b2) <- paste0(colnames(neighbor_mean5b2), "_S5b2")
neighbor_mean5b3 <- gen_mean_neighbor_mat(edges5b3, slide5b3_P2)
colnames(neighbor_mean5b3) <- paste0(colnames(neighbor_mean5b3), "_S5b3")
neighbor_mean6 <- gen_mean_neighbor_mat(edges6, slide6)
colnames(neighbor_mean6) <- paste0(colnames(neighbor_mean6), "_S6")

full_neighbor_mat <- do.call(cbind, list(neighbor_mean1, neighbor_mean2, 
                                         neighbor_mean3, neighbor_mean4, 
                                         neighbor_mean5, neighbor_mean5b2, 
                                         neighbor_mean5b3, neighbor_mean6))

zeros <- as(matrix(0, ncol = sum(pData(all_slides)$neighbor_count == 0), 
                   nrow = nrow(full_neighbor_mat)),'dgCMatrix')
row.names(zeros) <- row.names(full_neighbor_mat)
colnames(zeros) <- row.names(pData(all_slides)[pData(all_slides)$neighbor_count == 0,])

full_neighbor_mat <- cbind(full_neighbor_mat, zeros)

edges1$from <- paste0(edges1$from, "_S1")
edges1$to <- paste0(edges1$to, "_S1")
edges2$from <- paste0(edges2$from, "_S2")
edges2$to <- paste0(edges2$to, "_S2")
edges3$from <- paste0(edges3$from, "_S3")
edges3$to <- paste0(edges3$to, "_S3")
edges4$from <- paste0(edges4$from, "_S4")
edges4$to <- paste0(edges4$to, "_S4")
edges5$from <- paste0(edges5$from, "_S5")
edges5$to <- paste0(edges5$to, "_S5")
edges5b2$from <- paste0(edges5b2$from, "_S5b2")
edges5b2$to <- paste0(edges5b2$to, "_S5b2")
edges5b3$from <- paste0(edges5b3$from, "_S5b3")
edges5b3$to <- paste0(edges5b3$to, "_S5b3")
edges6$from <- paste0(edges6$from, "_S6")
edges6$to <- paste0(edges6$to, "_S6")

all_edges <- do.call(rbind, list(edges1, edges2, edges3, edges4, edges5, 
                                 edges5b2, edges5b3, edges6))

saveRDS(all_edges, file = paste0(int_path, "all_edges.RDS"))

rm(neighbor_mean1, neighbor_mean2, neighbor_mean3, 
   neighbor_mean4, neighbor_mean5,neighbor_mean5b2, neighbor_mean5b3, 
   neighbor_mean6,  edges1, edges2, edges3, edges4, edges5, edges5b2, edges5b3, 
   edges6, gr, gr_adj,neighbor_count1, neighbor_count2, neighbor_count3,
   neighbor_count4, neighbor_count5, neighbor_count5b2, neighbor_count5b3, 
   neighbor_count6, zeros, slide2, slide4, slide5b2_P2, slide5b3_P2, slide6,
   full_neighbor_mat, slide5_P1)
gc()
sessionInfo()
