options(width = 100)

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


all_sig <- read.table("~/Code/P02567_TBIO-3021_10X_pilot/data/import/Signatures/signatures_pdac_20230720.txt",
                      header = TRUE, sep = "\t")
all_sig <- tidyr::separate_rows(all_sig, Gene_list, sep = ",")
all_sig_gene <- all_sig %>% group_by(Gene_list) %>% 
  summarize(sig_list = paste(Signature_name, collapse = ","))

all_sig_clean <- all_sig[grepl("PDAC", all_sig$Signature_name),]

all_sig_gene_clean <- all_sig_clean %>% group_by(Gene_list) %>% 
  summarize(sig_list = paste(Signature_name, collapse = ","))

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

add_tissue_segs <- function(seur_obj, slide_num) {
  fov_met <- seur_obj@meta.data[seur_obj@meta.data$CenterX_global_px > 0,]
  fov.crop <- Crop(seur_obj[[paste0("Slide", slide_num)]], 
                   y = c(min(fov_met$CenterX_global_px),
                         max(fov_met$CenterX_global_px)), 
                   x = c(min(fov_met$CenterY_global_px), 
                         max(fov_met$CenterY_global_px)))
  seur_obj[[paste0("P1")]] <- fov.crop
  DefaultBoundary(seur_obj[[paste0("P1")]]) <- "segmentation"
  
  fov_met <- seur_obj@meta.data[seur_obj@meta.data$CenterX_global_px < 0,]
  fov.crop <- Crop(seur_obj[[paste0("Slide", slide_num)]], 
                   y = c(min(fov_met$CenterX_global_px),
                         max(fov_met$CenterX_global_px)), 
                   x = c(min(fov_met$CenterY_global_px), 
                         max(fov_met$CenterY_global_px)))
  seur_obj[[paste0("P2")]] <- fov.crop
  DefaultBoundary(seur_obj[[paste0("P2")]]) <- "segmentation"
  return(seur_obj)
}

all_meta <- read.csv(paste0(stash_base, "/processed_metadata.csv"))
all_meta$cell_number <- stringr::str_split_fixed(all_meta$cell_ID, "_", 4)[,4]

slide1 <- LoadNanostring(data.dir = paste0(stash_data, "/Run5558_Slide1"), 
                         fov = "Slide1")

slide1_meta <- all_meta[all_meta$slide_ID_numeric == 1,]
row.names(slide1_meta) <- paste0(slide1_meta$cell_number, "_", slide1_meta$fov)

meta <- read.csv(paste0(stash_data, 
                        "/Run5558_Slide1/Run5558_Slide1_metadata_file.csv"))
row.names(meta) <- paste0(meta$cell_ID, "_", meta$fov)

slide1@meta.data <- cbind(slide1@meta.data, meta[row.names(slide1@meta.data),])

slide1_meta <- slide1_meta[,!names(slide1_meta) %in% names(slide1@meta.data)]
slide1@meta.data <- cbind(slide1@meta.data, slide1_meta[row.names(slide1@meta.data),])

#Add tissue specific segmentations (internal structure thing for future)
slide1 <- add_tissue_segs(slide1, 1)

slide1$QC <- ifelse(slide1$nCount_Nanostring > 30 & 
                      slide1$nFeature_Nanostring > 10,
                    "Pass", "Fail")
#Save slide
saveRDS(slide1, file = paste0(int_path, "slide1.RDS"))

cds_s1 <- SeuratWrappers::as.cell_data_set(slide1[,slide1$QC == "Pass"], assay = NULL)
fData(cds_s1) <- data.frame(row.names = row.names(exprs(cds_s1)), 
                            gene_short_name = row.names(exprs(cds_s1)))
pData(cds_s1)$slide <- "S1"
pData(cds_s1)$position <- ifelse(pData(cds_s1)$CenterX_global_px > 0, 
                                 "P1", "P2")
pData(cds_s1)$tissue <- paste0(pData(cds_s1)$slide, "_", 
                               pData(cds_s1)$position)

rm(slide1)
gc()

slide2 <- LoadNanostring(data.dir = paste0(stash_data, "/Run5560_Slide2"), 
                         fov = "Slide2")

slide2_meta <- all_meta[all_meta$slide_ID_numeric == 2,]
row.names(slide2_meta) <- paste0(slide2_meta$cell_number, "_", slide2_meta$fov)

meta <- read.csv(paste0(stash_data, 
                        "/Run5560_Slide2/Run5560_Slide2_metadata_file.csv"))
row.names(meta) <- paste0(meta$cell_ID, "_", meta$fov)

slide2@meta.data <- cbind(slide2@meta.data, 
                          meta[row.names(slide2@meta.data),])

slide2_meta <- slide2_meta[,!names(slide2_meta) %in% names(slide2@meta.data)]
slide2@meta.data <- cbind(slide2@meta.data, 
                          slide2_meta[row.names(slide2@meta.data),])

#Add tissue specific segmentations (internal structure thing for future)
slide2 <- add_tissue_segs(slide2, 2)

slide2$QC <- ifelse(slide2$nCount_Nanostring > 30 & 
                      slide2$nFeature_Nanostring > 10,
                    "Pass", "Fail")

#Save slide
saveRDS(slide2, file = paste0(int_path, "slide2.RDS"))

cds_s2 <- SeuratWrappers::as.cell_data_set(slide2[,slide2$QC == "Pass"], assay = NULL)
fData(cds_s2) <- data.frame(row.names = row.names(exprs(cds_s2)), 
                            gene_short_name = row.names(exprs(cds_s2)))
pData(cds_s2)$slide <- "S2"
pData(cds_s2)$position <- ifelse(pData(cds_s2)$CenterX_global_px > 0, 
                                 "P1", "P2")
pData(cds_s2)$tissue <- paste0(pData(cds_s2)$slide, "_", 
                               pData(cds_s2)$position)

rm(slide2)
gc()

slide3 <- LoadNanostring(data.dir = paste0(stash_data, "/Run5583_Slide3"), 
                         fov = "Slide3")

slide3_meta <- all_meta[all_meta$slide_ID_numeric == 3,]
row.names(slide3_meta) <- paste0(slide3_meta$cell_number, "_", slide3_meta$fov)

meta <- read.csv(paste0(stash_data, "/Run5583_Slide3/Run5583_Slide3_metadata_file.csv"))
row.names(meta) <- paste0(meta$cell_ID, "_", meta$fov)

slide3@meta.data <- cbind(slide3@meta.data, meta[row.names(slide3@meta.data),])

slide3_meta <- slide3_meta[,!names(slide3_meta) %in% names(slide3@meta.data)]
slide3@meta.data <- cbind(slide3@meta.data, slide3_meta[row.names(slide3@meta.data),])

#Add tissue specific segmentations (internal structure thing for future)
slide3 <- add_tissue_segs(slide3, 3)

slide3$QC <- ifelse(slide3$nCount_Nanostring > 30 & 
                      slide3$nFeature_Nanostring > 10,
                    "Pass", "Fail")

#Save slide
saveRDS(slide3, file = paste0(int_path, "slide3.RDS"))

cds_S3 <- SeuratWrappers::as.cell_data_set(slide3[,slide3$QC == "Pass"], assay = NULL)
fData(cds_S3) <- data.frame(row.names = row.names(exprs(cds_S3)), 
                            gene_short_name = row.names(exprs(cds_S3)))

pData(cds_S3)$slide <- "S3"
pData(cds_S3)$position <- ifelse(pData(cds_S3)$CenterX_global_px > 0, 
                                 "P1", "P2")
pData(cds_S3)$tissue <- paste0(pData(cds_S3)$slide, "_", 
                               pData(cds_S3)$position)

rm(slide3)
gc()

slide4 <- LoadNanostring(data.dir = paste0(stash_data, "/Run5583_Slide4"), 
                         fov = "Slide4")

slide4_meta <- all_meta[all_meta$slide_ID_numeric == 4,]
row.names(slide4_meta) <- paste0(slide4_meta$cell_number, "_", slide4_meta$fov)

meta <- read.csv(paste0(stash_data, 
                        "/Run5583_Slide4/Run5583_Slide4_metadata_file.csv"))
row.names(meta) <- paste0(meta$cell_ID, "_", meta$fov)

slide4@meta.data <- cbind(slide4@meta.data, meta[row.names(slide4@meta.data),])

slide4_meta <- slide4_meta[,!names(slide4_meta) %in% names(slide4@meta.data)]
slide4@meta.data <- cbind(slide4@meta.data, slide4_meta[row.names(slide4@meta.data),])

#Add tissue specific segmentations (internal structure thing for future)
slide4 <- add_tissue_segs(slide4, 4)

slide4$QC <- ifelse(slide4$nCount_Nanostring > 30 & 
                      slide4$nFeature_Nanostring > 10,
                    "Pass", "Fail")

#Save slide
saveRDS(slide4, file = paste0(int_path, "slide4.RDS"))

cds_S4 <- SeuratWrappers::as.cell_data_set(slide4[,slide4$QC == "Pass"], assay = NULL)
fData(cds_S4) <- data.frame(row.names = row.names(exprs(cds_S4)), 
                            gene_short_name = row.names(exprs(cds_S4)))

pData(cds_S4)$slide <- "S4"
pData(cds_S4)$position <- ifelse(pData(cds_S4)$CenterX_global_px > 0, 
                                 "P1", "P2")
pData(cds_S4)$tissue <- paste0(pData(cds_S4)$slide, "_", 
                               pData(cds_S4)$position)

rm(slide4)
gc()

slide5 <- LoadNanostring(data.dir = paste0(stash_data, "/Run5592_Slide5"), 
                         fov = "Slide5")

slide5_meta <- all_meta[all_meta$slide_ID_numeric == 5,]
row.names(slide5_meta) <- paste0(slide5_meta$cell_number, "_", slide5_meta$fov)

meta <- read.csv(paste0(stash_data, 
                        "/Run5592_Slide5/Run5592_Slide5_metadata_file.csv"))
row.names(meta) <- paste0(meta$cell_ID, "_", meta$fov)

slide5@meta.data <- cbind(slide5@meta.data, meta[row.names(slide5@meta.data),])

slide5_meta <- slide5_meta[,!names(slide5_meta) %in% names(slide5@meta.data)]
slide5@meta.data <- cbind(slide5@meta.data, slide5_meta[row.names(slide5@meta.data),])

#Add tissue specific segmentations (internal structure thing for future)
slide5 <- add_tissue_segs(slide5, 5)

slide5$QC <- ifelse(slide5$nCount_Nanostring > 30 & 
                      slide5$nFeature_Nanostring > 10,
                    "Pass", "Fail")

#Save slide
saveRDS(slide5, file = paste0(int_path, "slide5.RDS"))

cds_S5 <- SeuratWrappers::as.cell_data_set(slide5[,slide5$QC == "Pass"], assay = NULL)
fData(cds_S5) <- data.frame(row.names = row.names(exprs(cds_S5)), 
                            gene_short_name = row.names(exprs(cds_S5)))

pData(cds_S5)$slide <- "S5"
pData(cds_S5)$position <- ifelse(pData(cds_S5)$CenterX_global_px > 0, 
                                 "P1", "P2")
pData(cds_S5)$tissue <- paste0(pData(cds_S5)$slide, "_", 
                               pData(cds_S5)$position)

rm(slide5)
gc()

slide5b2 <- LoadNanostring(data.dir = paste0(stash_data, "/Run5746_Slide5_b2"),
                           fov = "Slide5b2")

slide5b2_meta <- all_meta[all_meta$slide_ID_numeric == 6,]
row.names(slide5b2_meta) <- paste0(slide5b2_meta$cell_number, "_", slide5b2_meta$fov)

meta <- read.csv(paste0(stash_data, 
                        "/Run5746_Slide5_b2/Run5746_Slide5_b2_metadata_file.csv"))
row.names(meta) <- paste0(meta$cell_ID, "_", meta$fov)

slide5b2@meta.data <- cbind(slide5b2@meta.data, meta[row.names(slide5b2@meta.data),])

slide5b2_meta <- slide5b2_meta[,!names(slide5b2_meta) %in% names(slide5b2@meta.data)]
slide5b2@meta.data <- cbind(slide5b2@meta.data, slide5b2_meta[row.names(slide5b2@meta.data),])

#Add tissue specific segmentations (internal structure thing for future)
slide5b2 <- add_tissue_segs(slide5b2, "5b2")

slide5b2$QC <- ifelse(slide5b2$nCount_Nanostring > 30 & 
                        slide5b2$nFeature_Nanostring > 10,
                      "Pass", "Fail")

#Save slide
saveRDS(slide5b2, file = paste0(int_path, "slide5b2.RDS"))

cds_S5b2 <- SeuratWrappers::as.cell_data_set(slide5b2[,slide5b2$QC == "Pass"], assay = NULL)
fData(cds_S5b2) <- data.frame(row.names = row.names(exprs(cds_S5b2)), 
                              gene_short_name = row.names(exprs(cds_S5b2)))

pData(cds_S5b2)$slide <- "S5b2"
pData(cds_S5b2)$position <- ifelse(pData(cds_S5b2)$CenterX_global_px > 0, 
                                   "P1", "P2")
pData(cds_S5b2)$tissue <- paste0(pData(cds_S5b2)$slide, "_", 
                                 pData(cds_S5b2)$position)

rm(slide5b2)
gc()

slide5b3 <- LoadNanostring(data.dir = paste0(stash_data, "/Run5783_Slide5_b3"),
                           fov = "Slide5b3")

slide5b3_meta <- all_meta[all_meta$slide_ID_numeric == 7,]
row.names(slide5b3_meta) <- paste0(slide5b3_meta$cell_number, "_", slide5b3_meta$fov)

meta <- read.csv(paste0(stash_data, 
                        "/Run5783_Slide5_b3/Run5783_Slide5_b3_metadata_file.csv"))
row.names(meta) <- paste0(meta$cell_ID, "_", meta$fov)

slide5b3@meta.data <- cbind(slide5b3@meta.data, meta[row.names(slide5b3@meta.data),])

slide5b3_meta <- slide5b3_meta[,!names(slide5b3_meta) %in% names(slide5b3@meta.data)]
slide5b3@meta.data <- cbind(slide5b3@meta.data, slide5b3_meta[row.names(slide5b3@meta.data),])

#Add tissue specific segmentations (internal structure thing for future)
fov_met <- slide5b3@meta.data[slide5b3@meta.data$CenterX_global_px < 0,]
fov.crop <- Crop(slide5b3[[paste0("Slide5b3")]], 
                 y = c(min(fov_met$CenterX_global_px),
                       max(fov_met$CenterX_global_px)), 
                 x = c(min(fov_met$CenterY_global_px), 
                       max(fov_met$CenterY_global_px)))
slide5b3[[paste0("P2")]] <- fov.crop
DefaultBoundary(slide5b3[[paste0("P2")]]) <- "segmentation"

slide5b3$QC <- ifelse(slide5b3$nCount_Nanostring > 30 &
                        slide5b3$nFeature_Nanostring > 10,
                      "Pass", "Fail")

#Save slide
saveRDS(slide5b3, file = paste0(int_path, "slide5b3.RDS"))

cds_S5b3 <- SeuratWrappers::as.cell_data_set(slide5b3[,slide5b3$QC == "Pass"], assay = NULL)
fData(cds_S5b3) <- data.frame(row.names = row.names(exprs(cds_S5b3)), 
                              gene_short_name = row.names(exprs(cds_S5b3)))

pData(cds_S5b3)$slide <- "S5b3"
pData(cds_S5b3)$position <- ifelse(pData(cds_S5b3)$CenterX_global_px > 0, 
                                   "P1", "P2")
pData(cds_S5b3)$tissue <- paste0(pData(cds_S5b3)$slide, "_", 
                                 pData(cds_S5b3)$position)

rm(slide5b3)
gc()

slide6 <- LoadNanostring(data.dir = paste0(stash_data, "/Run5588_Slide6"), 
                         fov = "Slide6")

slide6_meta <- all_meta[all_meta$slide_ID_numeric == 8,]
row.names(slide6_meta) <- paste0(slide6_meta$cell_number, "_", slide6_meta$fov)

meta <- read.csv(paste0(stash_data, 
                        "/Run5588_Slide6/Run5588_Slide6_metadata_file.csv"))
row.names(meta) <- paste0(meta$cell_ID, "_", meta$fov)

slide6@meta.data <- cbind(slide6@meta.data, meta[row.names(slide6@meta.data),])

slide6_meta <- slide6_meta[,!names(slide6_meta) %in% names(slide6@meta.data)]
slide6@meta.data <- cbind(slide6@meta.data, slide6_meta[row.names(slide6@meta.data),])

#Add tissue specific segmentations (internal structure thing for future)
slide6 <- add_tissue_segs(slide6, 6)

slide6$QC <- ifelse(slide6$nCount_Nanostring > 30 & 
                      slide6$nFeature_Nanostring > 10,
                    "Pass", "Fail")

#Save slide
saveRDS(slide6, file = paste0(int_path, "slide6.RDS"))

cds_S6 <- SeuratWrappers::as.cell_data_set(slide6[,slide6$QC == "Pass"], assay = NULL)
fData(cds_S6) <- data.frame(row.names = row.names(exprs(cds_S6)), 
                            gene_short_name = row.names(exprs(cds_S6)))

pData(cds_S6)$slide <- "S6"
pData(cds_S6)$position <- ifelse(pData(cds_S6)$CenterX_global_px > 0, 
                                 "P1", "P2")
pData(cds_S6)$tissue <- paste0(pData(cds_S6)$slide, "_", 
                               pData(cds_S6)$position)

rm(slide6)
gc()

all_slides <- combine_cds(list(S1 = cds_s1,
                               S2 = cds_s2,
                               S3 = cds_S3,
                               S4 = cds_S4,
                               S5 = cds_S5,
                               S5b2 = cds_S5b2,
                               S5b3 = cds_S5b3,
                               S6 = cds_S6), sample_col_name = "slide")

# Exclude S5_P2 (failed), S5b2_P1 (repeat of S5_P1)
all_slides <- all_slides[,!pData(all_slides)$tissue %in% c("S5_P2", "S5b2_P1")]

saveRDS(all_slides, file = paste0(int_path, "all_slides.RDS"))

rm(all_slides)
gc()

## takes forever
all_slides <- readRDS(file = paste0(int_path, "all_slides.RDS"))
all_slides <- estimate_size_factors(all_slides)
all_slides <- preprocess_cds(all_slides)
all_slides <- align_cds(all_slides, alignment_group = "tissue") 
all_slides <- reduce_dimension(all_slides, preprocess_method = "Aligned")
all_slides <- cluster_cells(all_slides, reduction_method = "Aligned")
pData(all_slides)$aligned_clusts <- all_slides@clusters$Aligned$clusters
saveRDS(all_slides, file = paste0(int_path, "all_slides_int.RDS"))

all_slides <- readRDS(file = paste0(int_path, "all_slides_int.RDS"))
all_slides@clusters$UMAP <- all_slides@clusters$Aligned
all_seur <- as.Seurat(all_slides, assay = NULL)
all_seur <- Seurat::SCTransform(all_seur, assay="originalexp", variable.features.n=988,
                                return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")

##takes forever
all_seur <- Seurat::SetIdent(all_seur, value="aligned_clusts")
markers_matrix <- Seurat::FindAllMarkers(object = all_seur, assay="originalexp",
                                         min.pct=0.1, logfc.threshold=0.25,
                                         test.use = "wilcox")

markers_df <- markers_matrix %>%
  tibble::as_tibble(rownames = "Symbol") %>%
  dplyr::mutate(Symbol = gsub("\\.\\d+$", "", Symbol)) %>%
  dplyr::mutate(pct_1 = round(pct.1 * 100),
                pct_2 = round(pct.2 * 100),
                avg_logFC = round(avg_log2FC, 3),
                p_val_adj_neg_log10 = round(-log10(p_val_adj), 1),
                direction = ifelse(avg_logFC >= 0, "UP", "DN")) %>%
  dplyr::select(c(Symbol, cluster, direction, avg_logFC, p_val_adj_neg_log10, pct_1, pct_2))

markers_df <- markers_df %>%
  dplyr::mutate(IsSignificant = abs(avg_logFC) >= 0.5 & 
                  p_val_adj_neg_log10 >= 2 &
                  pct_1 >= 0.2*100)

# Filter markers down to signficant only
markers_filt_df <- markers_df %>%
  dplyr::filter(IsSignificant == TRUE) %>%
  dplyr::select(-IsSignificant)


markers_filt_df <- merge(markers_filt_df, all_sig_gene, 
                         by.x = "Symbol", by.y = "Gene_list", all.x = T)
write.table(markers_filt_df, 
            file = paste0(int_path, "20240109_cosmx_derived_signatures_all_clusters_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

pData(all_slides)$man_labels <- pData(all_slides)$aligned_clusts
pData(all_slides)$comb_labels <- pData(all_slides)$aligned_clusts
pData(all_slides)$comb_labels_clusts <- pData(all_slides)$aligned_clusts
levels(pData(all_slides)$comb_labels_clusts) <- list('Acinar (16,21)' = '16',
                                                     'Acinar (16,21)' = '21',
                                                     'Endocrine (23)' = '23',
                                                     'Ductal (12)' = '12',
                                                     'Tumor (4,6,8,9,10,15,17,19,26)' = '4',
                                                     'Tumor (4,6,8,9,10,15,17,19,26)' = '6',
                                                     'Tumor (4,6,8,9,10,15,17,19,26)' = '8',
                                                     'Tumor (4,6,8,9,10,15,17,19,26)' = '10',
                                                     'Tumor (4,6,8,9,10,15,17,19,26)' = '15',
                                                     'Tumor (4,6,8,9,10,15,17,19,26)' = '17',
                                                     'Tumor (4,6,8,9,10,15,17,19,26)' = '19',
                                                     'Tumor (4,6,8,9,10,15,17,19,26)' = '9',
                                                     'Tumor (4,6,8,9,10,15,17,19,26)' = '26',
                                                     'Fibroblast (1,3,20)' = '1',
                                                     'Fibroblast (1,3,20)' = '3',
                                                     'Fibroblast (1,3,20)' = '20',
                                                     'Stellate (18)' = '18',
                                                     'Endothelial (14,25)' = '14',
                                                     'Endothelial (14,25)' = '25',
                                                     'B cell/Plasma (13,24)' = '13',
                                                     'B cell/Plasma (13,24)' = '24',
                                                     'T cell (5)' = '5',
                                                     'Macrophage (2,11,28)' = '2',
                                                     'Macrophage (2,11,28)' = '11',
                                                     'Macrophage (2,11,28)' = '28',
                                                     'Dendritic (27)' = '27',
                                                     'Mast (22)' = '22',
                                                     'Low quality (7)' = '7')

levels(pData(all_slides)$comb_labels) <- list('Acinar' = '16',
                                              'Acinar' = '21',
                                              'Endocrine' = '23',
                                              'Ductal' = '12',
                                              'Tumor' = '4',
                                              'Tumor' = '6',
                                              'Tumor' = '8',
                                              'Tumor' = '10',
                                              'Tumor' = '15',
                                              'Tumor' = '17',
                                              'Tumor' = '19',
                                              'Tumor' = '9',
                                              'Tumor' = '26',
                                              'Fibroblast' = '1',
                                              'Fibroblast' = '3',
                                              'Fibroblast' = '20',
                                              'Stellate' = '18',
                                              'Endothelial' = '14',
                                              'Endothelial' = '25',
                                              'B cell/Plasma' = '13',
                                              'B cell/Plasma' = '24',
                                              'T cell' = '5',
                                              'Macrophage' = '2',
                                              'Macrophage' = '11',
                                              'Macrophage' = '28',
                                              'Dendritic' = '27',
                                              'Mast' = '22',
                                              'Low quality' = '7')

levels(pData(all_slides)$man_labels) <- list('Acinar 1' = '16',
                                             'Acinar 2' = '21',
                                             'Endocrine' = '23',
                                             'Ductal' = '12',
                                             'Tumor 1' = '4',
                                             'Tumor 2' = '6',
                                             'Tumor 3' = '8',
                                             'Tumor 4' = '10',
                                             'Tumor 5' = '15',
                                             'Tumor 6' = '17',
                                             'Tumor 7' = '19',
                                             'Tumor 8' = '9',
                                             'Tumor 9' = '26',
                                             'Fibroblast 1' = '1',
                                             'Fibroblast 2' = '3',
                                             'Fibroblast 3' = '20',
                                             'Stellate' = '18',
                                             'Endothelial 1' = '14',
                                             'Endothelial 2' = '25',
                                             'B cell/Plasma 1' = '13',
                                             'B cell/Plasma 2' = '24',
                                             'T cell' = '5',
                                             'Macrophage 1' = '2',
                                             'Macrophage 2' = '11',
                                             'Macrophage 3' = '28',
                                             'Dendritic' = '27',
                                             'Mast' = '22',
                                             'Low quality' = '7')

saveRDS(all_slides, file = paste0(int_path, "all_slides_int.RDS"))

all_seur <- as.Seurat(all_slides, assay = NULL)
all_seur <- Seurat::SCTransform(all_seur, assay="originalexp", variable.features.n=988,
                                return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")

all_seur <- Seurat::SetIdent(all_seur, value="comb_labels_clusts")
markers_matrix <- Seurat::FindAllMarkers(object = all_seur, assay="originalexp",
                                         min.pct=0.1, logfc.threshold=0.25,
                                         test.use = "wilcox")

markers_df <- markers_matrix %>%
  tibble::as_tibble(rownames = "Symbol") %>%
  dplyr::mutate(Symbol = gsub("\\.\\d+$", "", Symbol)) %>%
  dplyr::mutate(pct_1 = round(pct.1 * 100),
                pct_2 = round(pct.2 * 100),
                avg_logFC = round(avg_log2FC, 3),
                p_val_adj_neg_log10 = round(-log10(p_val_adj), 1),
                direction = ifelse(avg_logFC >= 0, "UP", "DN")) %>%
  dplyr::select(c(Symbol, cluster, direction, avg_logFC, p_val_adj_neg_log10, pct_1, pct_2))

markers_df <- markers_df %>%
  dplyr::mutate(IsSignificant = abs(avg_logFC) >= 0.5 & 
                  p_val_adj_neg_log10 >= 2 &
                  pct_1 >= 0.2*100)

# Filter markers down to signficant only
markers_filt_df <- markers_df %>%
  dplyr::filter(IsSignificant == TRUE) %>%
  dplyr::select(-IsSignificant)


markers_filt_df <- merge(markers_filt_df, all_sig_gene, 
                         by.x = "Symbol", by.y = "Gene_list", all.x = T)
write.table(markers_filt_df, 
            file = paste0(int_path,  "20240109_cosmx_derived_signatures_all_subtype_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

for(sig in type_sigs) {
  sig_genes <- all_sig[all_sig$Signature_name == sig,]$Gene_list
  sig_genes <- sig_genes[sig_genes %in% row.names(fData(all_slides))]
  all_seur <- UCell::AddModuleScore_UCell(all_seur, features = list(sig_genes),
                                          name = sig)
}

##takes forever
sig_table <- all_seur@meta.data[,paste0("signature_1", type_sigs)]
names(sig_table) <- gsub("signature_1", "", names(sig_table))
saveRDS(sig_table, file  = paste0(int_path,  "signature_table.RDS"))

sig_table$cluster <- all_seur$aligned_clusts

sig_summary <- sig_table %>% group_by(cluster) %>% 
  summarize(PDAC.P19.Acinar = mean(PDAC.P19.Acinar),
            PDAC.P19.Endocrine = mean(PDAC.P19.Endocrine),
            PDAC.P19.Ductal_1 = mean(PDAC.P19.Ductal_1),
            PDAC.collisson.classical = mean(PDAC.collisson.classical),
            PDAC.moffitt.basal = mean(PDAC.moffitt.basal),
            PDAC.P19.Fibroblast = mean(PDAC.P19.Fibroblast),
            PDAC.P19.Stellate = mean(PDAC.P19.Stellate),
            PDAC.P19.Endothelial = mean(PDAC.P19.Endothelial),
            PDAC.P19.Macrophage = mean(PDAC.P19.Macrophage),
            PDAC.P19.Tcell = mean(PDAC.P19.Tcell),
            PDAC.P19.Bcell = mean(PDAC.P19.Bcell),
            PDAC.U.Nervous = mean(PDAC.U.Nervous),
            PDAC.P19.Ductal_2 = mean(PDAC.P19.Ductal_2),
            Seurat.Proliferation = mean(Seurat.Proliferation))

write.table(sig_summary, 
            file = paste0(int_path,  "20240109_cosmx_mean_signatures_all_clusters_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

tumor <- all_slides[,all_slides$comb_labels %in% c("Ductal", "Tumor")]

tumor <- estimate_size_factors(tumor)
tumor <- preprocess_cds(tumor)
tumor <- align_cds(tumor, alignment_group = "tissue")
tumor <- reduce_dimension(tumor, preprocess_method = "Aligned")

 set.seed(250)
tumor <- cluster_cells(tumor, reduction_method = "Aligned", resolution = 2e-4)
pData(tumor)$tumor_aligned_clusts <- tumor@clusters$Aligned$clusters

tumor@clusters$UMAP <- tumor@clusters$Aligned
tumor_seur <- as.Seurat(tumor, assay = NULL)
tumor_seur <- NormalizeData(tumor_seur)
tumor_seur$mon_cluster <- monocle3::clusters(tumor)
peng_duct1 <- all_sig[all_sig$Signature_name == "PDAC.P19.Ductal_1",]$Gene_list
peng_duct1 <- peng_duct1[peng_duct1 %in% row.names(fData(tumor))]
tumor_seur <- UCell::AddModuleScore_UCell(tumor_seur, features = list(peng_duct1),
                                           name = "peng_duct1")
pData(tumor)$ucell_peng_ductal1_sig <- tumor_seur$signature_1peng_duct1
tumor_seur$ucell_peng_ductal1_sig <- tumor_seur$signature_1peng_duct1

saveRDS(tumor, file = paste0(int_path, "all_slides_tumor_integrated.RDS"))

tumor <- readRDS(file = paste0(int_path, "all_slides_tumor_integrated.RDS"))

tumor_seur <- as.Seurat(tumor, assay = NULL)
tumor_seur <- Seurat::SCTransform(tumor_seur, assay="originalexp", variable.features.n=988,
                                  return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")

tumor_seur <- Seurat::SetIdent(tumor_seur, value="tumor_aligned_clusts")
markers_matrix <- Seurat::FindAllMarkers(object = tumor_seur, assay="originalexp",
                                         min.pct=0.1, logfc.threshold=0.25,
                                         test.use = "wilcox")

markers_df <- markers_matrix %>%
  tibble::as_tibble(rownames = "Symbol") %>%
  dplyr::mutate(Symbol = gsub("\\.\\d+$", "", Symbol)) %>%
  dplyr::mutate(pct_1 = round(pct.1 * 100),
                pct_2 = round(pct.2 * 100),
                avg_logFC = round(avg_log2FC, 3),
                p_val_adj_neg_log10 = round(-log10(p_val_adj), 1),
                direction = ifelse(avg_logFC >= 0, "UP", "DN")) %>%
  dplyr::select(c(Symbol, cluster, direction, avg_logFC, p_val_adj_neg_log10, pct_1, pct_2))

markers_df <- markers_df %>%
  dplyr::mutate(IsSignificant = abs(avg_logFC) >= 0.5 & 
                  p_val_adj_neg_log10 >= 2 &
                  pct_1 >= 0.2*100)

# Filter markers down to signficant only
markers_filt_df <- markers_df %>%
  dplyr::filter(IsSignificant == TRUE) %>%
  dplyr::select(-IsSignificant)

markers_filt_df <- merge(markers_filt_df, all_sig_gene, 
                         by.x = "Symbol", by.y = "Gene_list", all.x = T)
write.table(markers_filt_df, 
            file = paste0(int_path, "20240109_cosmx_derived_signatures_tumor_clusters_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

PDAC.collisson.classical <- c("AGR2","ATP10B","CAPN8","CEACAM5","CEACAM6",
                              "ELF3","ERBB3","FOXQ1","FXYD3","GPRC5A","GPX2",
                              "LGALS4","MUC13","PLS1","S100P","SDR16C5",
                              "ST6GALNAC1","TFF1","TFF3","TMEM45B","TOX3",
                              "TSPAN8")


PDAC.moffitt.basal <- c("ANXA8L1","AREG","CST6","CTSV","DHRS9","FAM83A",
                        "FGFBP1","GPR87","KRT15","KRT17","KRT6A","KRT6C",
                        "KRT7","LEMD1","LY6D","S100A2","SCEL","SERPINB3",
                        "SERPINB4","SLC2A1","SPRR1B","SPRR3",
                        "TNS4","UCA1","VGLL1")
PDAC.collisson.classical_int <- intersect(row.names(fData(tumor)), 
                                          PDAC.collisson.classical)

PDAC.moffitt.basal_int <- intersect(row.names(fData(tumor)), 
                                    PDAC.moffitt.basal)
tumor_seur <- UCell::AddModuleScore_UCell(tumor_seur, 
                                          features = list(PDAC.collisson.classical_int), 
                                          name = "collisson_classical_sig")
tumor_seur <- UCell::AddModuleScore_UCell(tumor_seur, 
                                          features = list(PDAC.moffitt.basal_int), 
                                          name = "moffitt_basal_sig")
pData(tumor)$ucell_collisson_classical_sig <- tumor_seur$signature_1collisson_classical_sig
pData(tumor)$ucell_moffitt_basal_sig <- tumor_seur$signature_1moffitt_basal_sig

for_plot <- do.call(rbind, list(
  data.frame(tumor_aligned_clusts = pData(tumor)$tumor_aligned_clusts,
             signature_val = pData(tumor)$ucell_collisson_classical_sig,
             signature = "Classical Signature (Ucell)"),
  data.frame(tumor_aligned_clusts = pData(tumor)$tumor_aligned_clusts,
             signature_val = pData(tumor)$ucell_moffitt_basal_sig,
             signature = "Basal Signature (Ucell)"),
  data.frame(tumor_aligned_clusts = pData(tumor)$tumor_aligned_clusts,
             signature_val = pData(tumor)$ucell_peng_ductal1_sig,
             signature = "Ductal Signature (Ucell)")
))

ggplot(for_plot, aes(tumor_aligned_clusts, signature_val, fill = signature)) + 
  geom_boxplot(size = 0.2, outlier.size = 0.01) +
  scale_fill_manual(values=list("Ductal Signature (Ucell)" = "#4CB0E1",
                                "Basal Signature (Ucell)" = "#96257D",
                                "Classical Signature (Ucell)" = "#9A2626")) + 
  theme_bw(base_size = 6) + monocle3:::monocle_theme_opts() + 
  theme(legend.title = element_blank()) + 
  labs(x = "de novo cluster", y = "Ucell signature value")

pData(tumor)$man_labels <- pData(tumor)$tumor_aligned_clusts
pData(tumor)$comb_labels <- pData(tumor)$tumor_aligned_clusts
pData(tumor)$comb_labels_clusts <- pData(tumor)$tumor_aligned_clusts
levels(pData(tumor)$comb_labels_clusts) <- list('Ductal (6)' = '6',
                                                'Tumor Classical (1,2,3,5,7,8,11,14)' = '1',
                                                'Tumor Classical (1,2,3,5,7,8,11,14)' = '2',
                                                'Tumor Classical (1,2,3,5,7,8,11,14)' = '3',
                                                'Tumor Classical (1,2,3,5,7,8,11,14)' = '5',
                                                'Tumor Classical (1,2,3,5,7,8,11,14)' = '7',
                                                'Tumor Classical (1,2,3,5,7,8,11,14)' = '8',
                                                'Tumor Classical (1,2,3,5,7,8,11,14)' = '14',
                                                'Tumor Classical / Proliferative (9)' = '9',
                                                'Tumor Classical (1,2,3,5,7,8,11,14)' = '11',
                                                'Tumor Mixed (10)' = '10',
                                                'Tumor Mixed / Fibro adjacent (12)' = '12',
                                                'Tumor Basal (4,15)' = '4',
                                                'Tumor Basal (4,15)' = '15',
                                                'Tumor Other (13)' = '13')

levels(pData(tumor)$comb_labels) <- list('Ductal' = '6',
                                         'Tumor Classical 1' = '1',
                                         'Tumor Classical 2' = '2',
                                         'Tumor Classical 3' = '3',
                                         'Tumor Classical 4' = '5',
                                         'Tumor Classical 5' = '7',
                                         'Tumor Classical 6' = '8',
                                         'Tumor Classical 7' = '14',
                                         'Tumor Classical 8 / Proliferative' = '9',
                                         'Tumor Classical 9' = '11',
                                         'Tumor Mixed' = '10',
                                         'Tumor Mixed / Fibro adjacent' = '12',
                                         'Tumor Basal 1' = '4',
                                         'Tumor Basal 2' = '15',
                                         'Tumor Other' = '13')



tumor_seur <- as.Seurat(tumor, assay = NULL)
tumor_seur <- Seurat::SCTransform(tumor_seur, assay="originalexp", variable.features.n=988,
                                  return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")

tumor_seur <- Seurat::SetIdent(tumor_seur, value="comb_labels_clusts")
markers_matrix <- Seurat::FindAllMarkers(object = tumor_seur, assay="originalexp",
                                         min.pct=0.1, logfc.threshold=0.25,
                                         test.use = "wilcox")

markers_df <- markers_matrix %>%
  tibble::as_tibble(rownames = "Symbol") %>%
  dplyr::mutate(Symbol = gsub("\\.\\d+$", "", Symbol)) %>%
  dplyr::mutate(pct_1 = round(pct.1 * 100),
                pct_2 = round(pct.2 * 100),
                avg_logFC = round(avg_log2FC, 3),
                p_val_adj_neg_log10 = round(-log10(p_val_adj), 1),
                direction = ifelse(avg_logFC >= 0, "UP", "DN")) %>%
  dplyr::select(c(Symbol, cluster, direction, avg_logFC, p_val_adj_neg_log10, pct_1, pct_2))

markers_df <- markers_df %>%
  dplyr::mutate(IsSignificant = abs(avg_logFC) >= 0.5 & 
                  p_val_adj_neg_log10 >= 2 &
                  pct_1 >= 0.2*100)

# Filter markers down to signficant only
markers_filt_df <- markers_df %>%
  dplyr::filter(IsSignificant == TRUE) %>%
  dplyr::select(-IsSignificant)


markers_filt_df <- merge(markers_filt_df, all_sig_gene, 
                         by.x = "Symbol", by.y = "Gene_list", all.x = T)
write.table(markers_filt_df, 
            file = paste0(int_path, "20240109_cosmx_derived_signatures_tumor_subtype_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

saveRDS(tumor, file = paste0(int_path, "all_slides_tumor_integrated.RDS"))

sig_tumor <- sig_table[row.names(pData(tumor)),]
sig_tumor$cluster <- pData(tumor)$tumor_aligned_clusts
sig_tumor <- sig_tumor[!is.na(sig_tumor$PDAC.P19.Acinar),]
sig_summary <- sig_tumor %>% group_by(cluster) %>% 
  summarize(PDAC.P19.Ductal_1 = mean(PDAC.P19.Ductal_1),
            PDAC.collisson.classical = mean(PDAC.collisson.classical),
            PDAC.moffitt.basal = mean(PDAC.moffitt.basal),
            PDAC.P19.Fibroblast = mean(PDAC.P19.Fibroblast),
            PDAC.P19.Ductal_2 = mean(PDAC.P19.Ductal_2),
            Seurat.Proliferation = mean(Seurat.Proliferation))

write.table(sig_summary, 
            file = paste0(int_path,  "20240109_cosmx_mean_signatures_tumor_clusters_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

fibro <- all_slides[,all_slides$comb_labels %in% c("Fibroblast", "Stellate", "Endothelial")]

fibro <- estimate_size_factors(fibro)
fibro <- preprocess_cds(fibro)
fibro <- reduce_dimension(fibro, preprocess_method = "PCA")
plot_cells(fibro, color_cells_by = "tissue")

fibro <- align_cds(fibro, alignment_group = "tissue")
fibro <- reduce_dimension(fibro, preprocess_method = "Aligned")

set.seed(250)
fibro <- cluster_cells(fibro, reduction_method = "Aligned", resolution = 2e-4)
pData(fibro)$aligned_clusts <- fibro@clusters$Aligned$clusters
fibro@clusters$UMAP <- fibro@clusters$Aligned

saveRDS(fibro, file = paste0(int_path, "all_slides_fibro_integrated.RDS"))

fibro <- readRDS(file = paste0(int_path, "all_slides_fibro_integrated.RDS"))
fibro@clusters$UMAP <- fibro@clusters$Aligned
fibro_seur <- as.Seurat(fibro, assay = NULL)
fibro_seur <- Seurat::SCTransform(fibro_seur, assay="originalexp", variable.features.n=988,
                                  return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")

fibro_seur <- Seurat::SetIdent(fibro_seur, value="aligned_clusts")
markers_matrix <- Seurat::FindAllMarkers(object = fibro_seur, assay="originalexp",
                                         min.pct=0.1, logfc.threshold=0.25,
                                         test.use = "wilcox")

markers_df <- markers_matrix %>%
  tibble::as_tibble(rownames = "Symbol") %>%
  dplyr::mutate(Symbol = gsub("\\.\\d+$", "", Symbol)) %>%
  dplyr::mutate(pct_1 = round(pct.1 * 100),
                pct_2 = round(pct.2 * 100),
                avg_logFC = round(avg_log2FC, 3),
                p_val_adj_neg_log10 = round(-log10(p_val_adj), 1),
                direction = ifelse(avg_logFC >= 0, "UP", "DN")) %>%
  dplyr::select(c(Symbol, cluster, direction, avg_logFC, p_val_adj_neg_log10, pct_1, pct_2))

markers_df <- markers_df %>%
  dplyr::mutate(IsSignificant = abs(avg_logFC) >= 0.5 & 
                  p_val_adj_neg_log10 >= 2 &
                  pct_1 >= 0.2*100)

# Filter markers down to signficant only
markers_filt_df <- markers_df %>%
  dplyr::filter(IsSignificant == TRUE) %>%
  dplyr::select(-IsSignificant)


markers_filt_df <- merge(markers_filt_df, all_sig_gene, 
                         by.x = "Symbol", by.y = "Gene_list", all.x = T)
write.table(markers_filt_df, 
            file = paste0(int_path, "20240109_cosmx_derived_signatures_stroma_cluster_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

iCAF_sig <- c("IL6", "PDGFRA", "CXCL12", "CFD", "DPT", "LMNA", "AGTR1", "HAS1", "CXCL1", "CXCL2", "CCL2", "IL8")
myCAF_sig <- c("ACTA2", "TAGLN", "MMP11", "MYL9", "HOPX", "POSTN", "TPM1", "TPM2")
panCAF_sig <- c("COL1A1", "FAP", "PDPN", "DCN", "VIM")
fibro_seur <- UCell::AddModuleScore_UCell(fibro_seur, features=list(iCAF_sig), name="iCAF")
fibro_seur <- UCell::AddModuleScore_UCell(fibro_seur, features=list(myCAF_sig), name="myCAF")
fibro_seur <- UCell::AddModuleScore_UCell(fibro_seur, features=list(panCAF_sig), name="panCAF")

pData(fibro)$icaf_ucell <- fibro_seur$signature_1iCAF
pData(fibro)$mycaf_ucell <- fibro_seur$signature_1myCAF
pData(fibro)$pancaf_ucell <- fibro_seur$signature_1panCAF

fibro_seur <- UCell::AddModuleScore_UCell(fibro_seur, 
                                          features=list(all_sig[all_sig$Signature_name == "PDAC.moffitt.normalstroma",]$Gene_list), 
                                          name="PDAC.moffitt.normalstroma")
fibro_seur <- UCell::AddModuleScore_UCell(fibro_seur, 
                                          features=list(all_sig[all_sig$Signature_name == "PDAC.moffitt.activatedstroma",]$Gene_list), 
                                          name="PDAC.moffitt.activatedstroma")

pData(fibro)$act_strom_ucell <- fibro_seur$signature_1PDAC.moffitt.activatedstroma
pData(fibro)$norm_strom_ucell <- fibro_seur$signature_1PDAC.moffitt.normalstroma

temp <- fibro[,sample(1:nrow(pData(fibro)), nrow(pData(fibro)), 
                      replace = FALSE)]

p1 <- plot_cells(temp, color_cells_by = "act_strom_ucell", 
                 cell_size = 0.1, alpha = 0.4) + 
  ggtitle("Activated stroma") + theme_classic(base_size = 6) + 
  NoLegend() + scale_color_distiller(palette = "OrRd", direction = 1)

p2 <- plot_cells(temp, color_cells_by = "norm_strom_ucell", 
                 cell_size = 0.1, alpha = 0.4) + 
  ggtitle("Normal stroma") + theme_classic(base_size = 6) + 
  NoLegend() + scale_color_distiller(palette = "OrRd", direction = 1)

(p1 | p2)

PDAC.U.Nervous <- c("CDH19","CRYAB","PMP22","S100B","SCN7A","SOX10","ZEB2")
fibro_seur <- UCell::AddModuleScore_UCell(fibro_seur, 
                                          features=list(all_sig[all_sig$Signature_name == "PDAC.P19.Fibroblast",]$Gene_list), 
                                          name="PDAC.P19.Fibroblast")
fibro_seur <- UCell::AddModuleScore_UCell(fibro_seur, 
                                          features=list(all_sig[all_sig$Signature_name == "PDAC.P19.Stellate",]$Gene_list), 
                                          name="PDAC.P19.Stellate")
fibro_seur <- UCell::AddModuleScore_UCell(fibro_seur, 
                                          features=list(all_sig[all_sig$Signature_name == "PDAC.P19.Endothelial",]$Gene_list), 
                                          name="PDAC.P19.Endothelial")
fibro_seur <- UCell::AddModuleScore_UCell(fibro_seur, 
                                          features=list(PDAC.U.Nervous), 
                                          name="PDAC.U.Nervous")

pData(fibro)$PDAC.P19.Fibroblast_ucell <- fibro_seur$signature_1PDAC.P19.Fibroblast
pData(fibro)$PDAC.P19.Stellate_ucell <- fibro_seur$signature_1PDAC.P19.Stellate
pData(fibro)$PDAC.P19.Endothelial_ucell <- fibro_seur$signature_1PDAC.P19.Endothelial
pData(fibro)$PDAC.U.Nervous_ucell <- fibro_seur$signature_1PDAC.U.Nervous


for_plot <- do.call(rbind, list(
  data.frame(fibro_aligned_clusts = pData(fibro)$aligned_clusts,
             signature_val = pData(fibro)$PDAC.P19.Fibroblast_ucell,
             signature = "Fibroblast Signature (Ucell)"),
  data.frame(fibro_aligned_clusts = pData(fibro)$aligned_clusts,
             signature_val = pData(fibro)$PDAC.P19.Stellate_ucell,
             signature = "Stellate Signature (Ucell)"),
  data.frame(fibro_aligned_clusts = pData(fibro)$aligned_clusts,
             signature_val = pData(fibro)$PDAC.P19.Endothelial_ucell,
             signature = "Endothelial Signature (Ucell)"),
  data.frame(fibro_aligned_clusts = pData(fibro)$aligned_clusts,
             signature_val = pData(fibro)$PDAC.U.Nervous_ucell,
             signature = "Nerves Signature (Ucell)")
))

ggplot(for_plot, aes(fibro_aligned_clusts, signature_val, fill = signature)) + 
  geom_boxplot(size = 0.2, outlier.size = 0.01) +
  scale_fill_manual(values=list("Fibroblast Signature (Ucell)" = "#E5E5AA",
                                "Stellate Signature (Ucell)" = "#A0D5B5",
                                "Endothelial Signature (Ucell)" = "#9B5E2C",
                                "Nerves Signature (Ucell)" = "#FAE573")) + 
  theme_bw(base_size = 6) + monocle3:::monocle_theme_opts() + 
  theme(legend.title = element_blank()) + 
  labs(x = "de novo cluster", y = "Ucell signature value")

table(pData(fibro)$aligned_clusts)

pData(fibro)$man_labels <- pData(fibro)$aligned_clusts
pData(fibro)$comb_labels <- pData(fibro)$aligned_clusts
pData(fibro)$comb_labels_clusts <- pData(fibro)$aligned_clusts
levels(pData(fibro)$comb_labels_clusts) <- list('Fibroblast (1,2,4,5,7,9,10,13,14,15)' = '15',
                                                'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = '1',
                                                'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = '2',
                                                'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = '5',
                                                'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = '14',
                                                'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = '13',
                                                'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = '10',
                                                'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = '7',
                                                'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = '9',
                                                'Fibroblast (1,2,4,5,7,9,10,13,14,15)' = '4',
                                                'Stellate 1 (6)' = '6',
                                                'Stellate 2 (8)' = '8',
                                                'Endothelial (3,12)' = '3',
                                                'Endothelial (3,12)' = '12',
                                                'Fibroblasts/IFNa pathway/ resistance (11)' = '11',
                                                'Nerve (16)' = '16')

levels(pData(fibro)$comb_labels) <- list('Fibroblasts 1' = '15',
                                         'Fibroblasts 2' = '1',
                                         'Fibroblasts 3' = '2',
                                         'Fibroblasts 4' = '5',
                                         'Fibroblasts 5' = '14',
                                         'Fibroblasts 6' = '13',
                                         'Fibroblasts 7' = '10',
                                         'Fibroblasts 8' = '7',
                                         'Fibroblasts 9' = '9',
                                         'Fibroblasts 10' = '4',
                                         'Stellate 1' = '6',
                                         'Stellate 2' = '8',
                                         'Endothelial 1' = '3',
                                         'Endothelial 2' = '12',
                                         'Fibroblasts/IFNa pathway/ resistance' = '11',
                                         'Nerve' = '16')

fibro_seur <- as.Seurat(fibro, assay = NULL)
fibro_seur <- Seurat::SCTransform(fibro_seur, assay="originalexp", variable.features.n=988,
                                return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")

fibro_seur <- Seurat::SetIdent(fibro_seur, value="comb_labels_clusts")
markers_matrix <- Seurat::FindAllMarkers(object = fibro_seur, assay="originalexp",
                                         min.pct=0.1, logfc.threshold=0.25,
                                         test.use = "wilcox")

markers_df <- markers_matrix %>%
  tibble::as_tibble(rownames = "Symbol") %>%
  dplyr::mutate(Symbol = gsub("\\.\\d+$", "", Symbol)) %>%
  dplyr::mutate(pct_1 = round(pct.1 * 100),
                pct_2 = round(pct.2 * 100),
                avg_logFC = round(avg_log2FC, 3),
                p_val_adj_neg_log10 = round(-log10(p_val_adj), 1),
                direction = ifelse(avg_logFC >= 0, "UP", "DN")) %>%
  dplyr::select(c(Symbol, cluster, direction, avg_logFC, p_val_adj_neg_log10, pct_1, pct_2))

markers_df <- markers_df %>%
  dplyr::mutate(IsSignificant = abs(avg_logFC) >= 0.5 & 
                  p_val_adj_neg_log10 >= 2 &
                  pct_1 >= 0.2*100)

# Filter markers down to signficant only
markers_filt_df <- markers_df %>%
  dplyr::filter(IsSignificant == TRUE) %>%
  dplyr::select(-IsSignificant)


markers_filt_df <- merge(markers_filt_df, all_sig_gene, 
                         by.x = "Symbol", by.y = "Gene_list", all.x = T)
write.table(markers_filt_df, 
            file = paste0(int_path,  "20240109_cosmx_derived_signatures_stroma_subtype_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

saveRDS(fibro, file = paste0(int_path, "all_slides_fibro_integrated.RDS"))

sig_fibro <- sig_table[row.names(pData(fibro)),]
sig_fibro$cluster <- pData(fibro)$aligned_clusts
sig_fibro <- sig_fibro[!is.na(sig_fibro$PDAC.P19.Acinar),]
sig_summary <- sig_fibro %>% group_by(cluster) %>% 
  summarize(PDAC.P19.Fibroblast = mean(PDAC.P19.Fibroblast),
            PDAC.P19.Stellate = mean(PDAC.P19.Stellate),
            PDAC.P19.Endothelial = mean(PDAC.P19.Endothelial),
            PDAC.U.Nervous = mean(PDAC.U.Nervous),
            Seurat.Proliferation = mean(Seurat.Proliferation))

write.table(sig_summary, 
            file = paste0(int_path,  "20240109_cosmx_mean_signatures_fibro_clusters_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

imm <- all_slides[,all_slides$comb_labels %in% c("Macrophage", "T cell", "B cell/Plasma", "Mast", "Dendritic")]

imm <- estimate_size_factors(imm)
imm <- preprocess_cds(imm)
imm <- reduce_dimension(imm, preprocess_method = "PCA")

plot_cells(imm, color_cells_by = "tissue", label_cell_groups = F)

imm <- align_cds(imm, alignment_group = "tissue")
imm <- reduce_dimension(imm, preprocess_method = "Aligned")

set.seed(250)
imm <- cluster_cells(imm, reduction_method = "Aligned", resolution = 3e-4)
pData(imm)$aligned_clusts <- imm@clusters$Aligned$clusters
imm@clusters$UMAP <- imm@clusters$Aligned
saveRDS(imm, file = paste0(int_path, "all_slides_immune_integrated.RDS"))

imm <- readRDS(file = paste0(int_path, "all_slides_immune_integrated.RDS"))

imm@clusters$UMAP <- imm@clusters$Aligned
imm_seur <- as.Seurat(imm, assay = NULL)
imm_seur <- Seurat::SCTransform(imm_seur, assay="originalexp", variable.features.n=988,
                                return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")

imm_seur <- Seurat::SetIdent(imm_seur, value="aligned_clusts")
markers_matrix <- Seurat::FindAllMarkers(object = imm_seur, assay="originalexp",
                                         min.pct=0.1, logfc.threshold=0.25,
                                         test.use = "wilcox")

markers_df <- markers_matrix %>%
  tibble::as_tibble(rownames = "Symbol") %>%
  dplyr::mutate(Symbol = gsub("\\.\\d+$", "", Symbol)) %>%
  dplyr::mutate(pct_1 = round(pct.1 * 100),
                pct_2 = round(pct.2 * 100),
                avg_logFC = round(avg_log2FC, 3),
                p_val_adj_neg_log10 = round(-log10(p_val_adj), 1),
                direction = ifelse(avg_logFC >= 0, "UP", "DN")) %>%
  dplyr::select(c(Symbol, cluster, direction, avg_logFC, p_val_adj_neg_log10, pct_1, pct_2))

markers_df <- markers_df %>%
  dplyr::mutate(IsSignificant = abs(avg_logFC) >= 0.5 & 
                  p_val_adj_neg_log10 >= 2 &
                  pct_1 >= 0.2*100)

# Filter markers down to signficant only
markers_filt_df <- markers_df %>%
  dplyr::filter(IsSignificant == TRUE) %>%
  dplyr::select(-IsSignificant)


markers_filt_df <- merge(markers_filt_df, all_sig_gene, 
                         by.x = "Symbol", by.y = "Gene_list", all.x = T)
write.table(markers_filt_df, 
            file = paste0(int_path, "20240109_cosmx_derived_signatures_immune_clusters_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

imm_seur <- UCell::AddModuleScore_UCell(imm_seur, 
                                        features=list(all_sig[all_sig$Signature_name == "PDAC.P19.Macrophage",]$Gene_list), 
                                        name="PDAC.P19.Macrophage")
imm_seur <- UCell::AddModuleScore_UCell(imm_seur, 
                                        features=list(all_sig[all_sig$Signature_name == "PDAC.P19.Tcell",]$Gene_list), 
                                        name="PDAC.P19.Tcell")
imm_seur <- UCell::AddModuleScore_UCell(imm_seur, 
                                        features=list(all_sig[all_sig$Signature_name == "PDAC.P19.Bcell",]$Gene_list), 
                                        name="PDAC.P19.Bcell")

pData(imm)$PDAC.P19.Macrophage_ucell <- imm_seur$signature_1PDAC.P19.Macrophage
pData(imm)$PDAC.P19.Tcell_ucell <- imm_seur$signature_1PDAC.P19.Tcell
pData(imm)$PDAC.P19.Bcell_ucell <- imm_seur$signature_1PDAC.P19.Bcell

for_plot <- do.call(rbind, list(
  data.frame(aligned_clusts = pData(imm)$aligned_clusts,
             signature_val = pData(imm)$PDAC.P19.Macrophage_ucell,
             signature = "Macrophage Signature (Ucell)"),
  data.frame(aligned_clusts = pData(imm)$aligned_clusts,
             signature_val = pData(imm)$PDAC.P19.Tcell_ucell,
             signature = "T cell Signature (Ucell)"),
  data.frame(aligned_clusts = pData(imm)$aligned_clusts,
             signature_val = pData(imm)$PDAC.P19.Bcell_ucell,
             signature = "B cell Signature (Ucell)")
))

ggplot(for_plot, aes(aligned_clusts, signature_val, fill = signature)) + 
  geom_boxplot(size = 0.2, outlier.size = 0.01) +
  scale_fill_manual(values=list("Macrophage Signature (Ucell)" = "#FFE269",
                                "T cell Signature (Ucell)" = "#FF9F2C",
                                "B cell Signature (Ucell)" = "#DF7126")) + 
  theme_bw(base_size = 6) + monocle3:::monocle_theme_opts() + 
  theme(legend.title = element_blank()) + 
  labs(x = "de novo cluster", y = "Ucell signature value")

pData(imm)$man_labels <- pData(imm)$aligned_clusts
pData(imm)$comb_labels <- pData(imm)$aligned_clusts
pData(imm)$comb_labels_clusts <- pData(imm)$aligned_clusts
levels(pData(imm)$comb_labels_clusts) <- list('Plasma (7,9,17)' = '7',
                                              'Plasma (7,9,17)' = '9',
                                              'Plasma (7,9,17)' = '17',
                                              'B cell (16)' = '16',
                                              '? (23)' = '23',
                                              'Dendritic (14)' = '14',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '1',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '2',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '5',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '8',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '10',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '11',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '13',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '18',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '20',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '21',
                                              'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = '22',
                                              'T cell (3,4,6,15,19)' = '3',
                                              'T cell (3,4,6,15,19)' = '4',
                                              'T cell (3,4,6,15,19)' = '6',
                                              'T cell (3,4,6,15,19)' = '15',
                                              'T cell (3,4,6,15,19)' = '19',
                                              'Mast (12)' = '12')

levels(pData(imm)$comb_labels) <- list('Plasma 1' = '7',
                                       'Plasma 2' = '9',
                                       'Plasma 3' = '17',
                                       'B cell' = '16',
                                       '?' = '23',
                                       'Macrophage 1' = '1',
                                       'Macrophage 2' = '2',
                                       'Macrophage 3' = '5',
                                       'Macrophage 4' = '8',
                                       'Macrophage 5' = '10',
                                       'Macrophage 6' = '11',
                                       'Macrophage 7' = '13',
                                       'Macrophage 8' = '18',
                                       'Macrophage 9' = '20',
                                       'Macrophage 10' = '21',
                                       'Macrophage 11' = '22',
                                       'Dendritic' = '14',
                                       'T cell 1' = '3',
                                       'T cell 2' = '4',
                                       'T cell 3' = '6',
                                       'T cell 4' = '15',
                                       'T cell 5 / Prolif' = '19',
                                       'Mast' = '12')

imm_seur <- as.Seurat(imm, assay = NULL)
imm_seur <- Seurat::SCTransform(imm_seur, assay="originalexp", variable.features.n=988,
                                return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")

imm_seur <- Seurat::SetIdent(imm_seur, value="comb_labels_clusts")
markers_matrix <- Seurat::FindAllMarkers(object = imm_seur, assay="originalexp",
                                         min.pct=0.1, logfc.threshold=0.25,
                                         test.use = "wilcox")

markers_df <- markers_matrix %>%
  tibble::as_tibble(rownames = "Symbol") %>%
  dplyr::mutate(Symbol = gsub("\\.\\d+$", "", Symbol)) %>%
  dplyr::mutate(pct_1 = round(pct.1 * 100),
                pct_2 = round(pct.2 * 100),
                avg_logFC = round(avg_log2FC, 3),
                p_val_adj_neg_log10 = round(-log10(p_val_adj), 1),
                direction = ifelse(avg_logFC >= 0, "UP", "DN")) %>%
  dplyr::select(c(Symbol, cluster, direction, avg_logFC, p_val_adj_neg_log10, pct_1, pct_2))

markers_df <- markers_df %>%
  dplyr::mutate(IsSignificant = abs(avg_logFC) >= 0.5 & 
                  p_val_adj_neg_log10 >= 2 &
                  pct_1 >= 0.2*100)

# Filter markers down to signficant only
markers_filt_df <- markers_df %>%
  dplyr::filter(IsSignificant == TRUE) %>%
  dplyr::select(-IsSignificant)


markers_filt_df <- merge(markers_filt_df, all_sig_gene, 
                         by.x = "Symbol", by.y = "Gene_list", all.x = T)
write.table(markers_filt_df, 
            file = paste0(int_path,  "20240109_cosmx_derived_signatures_immune_subtype_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

saveRDS(imm, file = paste0(int_path, "all_slides_immune_integrated.RDS"))

sig_imm <- sig_table[row.names(pData(imm)),]
sig_imm$cluster <- pData(imm)$aligned_clusts
sig_imm <- sig_imm[!is.na(sig_imm$PDAC.P19.Acinar),]
sig_summary <- sig_imm %>% group_by(cluster) %>% 
  summarize(PDAC.P19.Macrophage = mean(PDAC.P19.Macrophage),
            PDAC.P19.Tcell = mean(PDAC.P19.Tcell),
            PDAC.P19.Bcell = mean(PDAC.P19.Bcell),
            Seurat.Proliferation = mean(Seurat.Proliferation))

write.table(sig_summary, 
            file = paste0(int_path,  "20240109_cosmx_mean_signatures_immune_clusters_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

tcell_sub <- imm[,pData(imm)$comb_labels_clusts == "T cell (3,4,6,15,19)"]
tcell_sub <- estimate_size_factors(tcell_sub)
tcell_sub <- preprocess_cds(tcell_sub)
tcell_sub <- align_cds(tcell_sub, alignment_group = "tissue")
tcell_sub <- reduce_dimension(tcell_sub, preprocess_method = "Aligned")
tcell_sub <- cluster_cells(tcell_sub, reduction_method = "Aligned", resolution = 6e-4)#resolution = 1e-3
pData(tcell_sub)$aligned_clusts <- tcell_sub@clusters$Aligned$clusters
tcell_sub@clusters$UMAP <- tcell_sub@clusters$Aligned

saveRDS(tcell_sub, file = paste0(int_path, "all_slides_tcell_integrated.RDS"))

tcell_seur <- as.Seurat(tcell_sub, assay = NULL)
tcell_seur <- Seurat::SCTransform(tcell_seur, assay="originalexp", variable.features.n=988,
                                return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")

tcell_seur <- Seurat::SetIdent(tcell_seur, value="aligned_clusts")
markers_matrix <- Seurat::FindAllMarkers(object = tcell_seur, assay="originalexp",
                                         min.pct=0.1, logfc.threshold=0.25,
                                         test.use = "wilcox")

markers_df <- markers_matrix %>%
  tibble::as_tibble(rownames = "Symbol") %>%
  dplyr::mutate(Symbol = gsub("\\.\\d+$", "", Symbol)) %>%
  dplyr::mutate(pct_1 = round(pct.1 * 100),
                pct_2 = round(pct.2 * 100),
                avg_logFC = round(avg_log2FC, 3),
                p_val_adj_neg_log10 = round(-log10(p_val_adj), 1),
                direction = ifelse(avg_logFC >= 0, "UP", "DN")) %>%
  dplyr::select(c(Symbol, cluster, direction, avg_logFC, p_val_adj_neg_log10, pct_1, pct_2))

markers_df <- markers_df %>%
  dplyr::mutate(IsSignificant = abs(avg_logFC) >= 0.5 & 
                  p_val_adj_neg_log10 >= 2 &
                  pct_1 >= 0.2*100)

# Filter markers down to signficant only
markers_filt_df <- markers_df %>%
  dplyr::filter(IsSignificant == TRUE) %>%
  dplyr::select(-IsSignificant)


markers_filt_df <- merge(markers_filt_df, all_sig_gene, 
                         by.x = "Symbol", by.y = "Gene_list", all.x = T)
write.table(markers_filt_df, 
            file = paste0(int_path, "20240109_cosmx_derived_signatures_tcells_Seurat.txt"), 
            quote = FALSE, row.names = FALSE,sep = "\t") 

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

## One set with minimal type
for(slide_num in c(1,2,3,4,5,'5b2', '5b3', 6)) {
  slide_list[[slide_num]]@meta.data[['tme_type']] <- pData(all_slides)[paste0(row.names(slide_list[[slide_num]]@meta.data), 
                                                               "_S", slide_num),][['comb_labels']]
  temp <- as.data.frame(slide_list[[slide_num]]@images[[paste0("Slide", slide_num)]]@boundaries$centroids@coords)
  temp$CELL_ID <- slide_list[[slide_num]]@images[[paste0("Slide", slide_num)]]@boundaries$centroids@cells
  names(temp) <- c("X", "Y", "CELL_ID")
  temp <- temp[,c("CELL_ID", "X", "Y")]
  temp$X <- round(temp$X + abs(min(temp$X)) + 1, 0)
  temp$Y <- round(temp$Y + abs(min(temp$Y)) + 1, 0)
  temp$fov <- stringr::str_split_fixed(temp$CELL_ID, "_", 2)[,2]
  for(fov in unique(temp$fov)) {
      write.csv(temp[temp$fov == fov, c(1:3)], quote = FALSE, row.names = FALSE, 
                file = paste0(home_path, "data/intermediate/TME_analysis/tme_type/slide", 
                              slide_num, "_FOV",fov, ".cell_data.csv"))
  }
  temp <- data.frame(CELL_ID = row.names(slide_list[[slide_num]]@meta.data),
                     CELL_TYPE = slide_list[[slide_num]]@meta.data$tme_type)
  temp$CELL_TYPE[is.na(temp$CELL_TYPE)] <- "Excluded"
  temp$fov <- stringr::str_split_fixed(temp$CELL_ID, "_", 2)[,2]
  for(fov in unique(temp$fov)) {
    write.csv(temp[temp$fov == fov, c(1:2)], quote = FALSE, row.names = FALSE, 
              file = paste0(home_path, "/data/intermediate/TME_analysis/tme_type/slide", 
                            slide_num, "_FOV",fov, ".cell_types.csv"))
  }
}

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

rm(slide1, slide2, slide3, slide4, slide5, slide5_P1, slide5b2, slide5b2_P2, 
   slide5b3, slide5b3_P2, slide6)
gc()

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
