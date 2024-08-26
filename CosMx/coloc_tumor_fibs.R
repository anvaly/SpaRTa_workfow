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



home_path <- "/out/cosmx_TAP/"
int_path <- "/out/cosmx_TAP/intermediate/" 

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

tumor_cols <- list('Ductal (6)' = "#4CB0E1",
                   'Tumor Classical (1,2,3,5,7,8,11,14)' = "#F16666",
                   'Tumor Classical / Proliferative (9)' = "#9A2626",
                   'Tumor Classical / Mixed (11)' = "#96257D",
                   'Tumor Mixed (10)' = '#CA3FAE',
                   'Tumor Mixed / Fibro adjacent (12)' = "#2C9FAF",
                   'Tumor Basal (4,15)' = "#681D59",
                   'Tumor Other (13)' = "lightgrey")

temp_cols <- list('Acinar' = "#056DB5",
                  'Endocrine' = "#153C65",
                  'Ductal' = "#4CB0E1",
                  'Tumor' = "#9A2626",
                  'Fibroblast' = "#E5E5AA",
                  'Stellate 1' = "#184928",
                  'Stellate 2' = "#A0D5B5",
                  'Endothelial' = "#358E5B",
                  'Immune' = "#D66100",
                  'Low quality' = "lightgray",
                  "Nerve" = "#939393")

tumor <- readRDS(file = paste0(int_path, "all_slides_tumor_integrated.RDS"))

pData(tumor)$slide_fov <- paste0(pData(tumor)$tissue, "_", pData(tumor)$fov)
fib_list <- as.data.frame(table(pData(tumor)[,c("slide_fov", "comb_labels_clusts")]))
fib_list <- fib_list[fib_list$comb_labels_clusts == "Tumor Mixed / Fibro adjacent (12)",]
fib_list <- fib_list[order(fib_list$Freq, decreasing = T),]

temp_cols <- c(temp_cols, all_cols, tumor_cols)

all_slides <- readRDS(file = paste0(int_path, "all_slides_int.RDS"))
all_slides@clusters$UMAP <- all_slides@clusters$Aligned
all_slides <- all_slides[,!(pData(all_slides)$tissue == "S2_P1" & pData(all_slides)$fov %in% c(22, 23))]


pData(tumor)$comb_labels_clusts <- as.character(pData(tumor)$comb_labels_clusts)
pData(all_slides)$tumor_comb_labels <- ""

pData(all_slides)$tumor_comb_labels <- pData(tumor)[row.names(pData(all_slides)),]$comb_labels_clusts

pData(all_slides)$tumor_labels <- ifelse(!is.na(pData(all_slides)$tumor_comb_labels), as.character(pData(all_slides)$tumor_comb_labels), 
                                         as.character(pData(all_slides)$comb_labels_clusts))
pData(all_slides)$tumor_labels[pData(all_slides)$tumor_labels %in% c("B cell/Plasma (13,24)", "Dendritic (27)", "Macrophage (2,11,28)", "Mast (22)", "T cell (5)")] <- "Immune"


colsvec=unlist(temp_cols)
make_fig3_fib_plots <- function(slide_num, fov_num, slide_obj, minval = 0.25) {
    fov_met <- slide_obj@meta.data[slide_obj@meta.data$fov == fov_num,]
    fov.crop <- Crop(slide_obj[[paste0("Slide", slide_num)]], 
                     y = c(min(fov_met$CenterX_global_px),
                           max(fov_met$CenterX_global_px)), 
                     x = c(min(fov_met$CenterY_global_px), 
                           max(fov_met$CenterY_global_px)))
    slide_obj[[paste0("FOV", fov_num)]] <- fov.crop
    DefaultBoundary(slide_obj[[paste0("FOV", fov_num)]]) <- "segmentation"
  
    print(ImageDimPlot(slide_obj,
                       fov = paste0("FOV", fov_num),
                       molecules = c("CEACAM6", "KRT6A", "COL1A1"),
                       nmols = 100000,
                       dark.background = FALSE,
                       alpha = 0.6,
                       border.size = 0.1,
                       border.color = "darkgray",
                       cols = "lightgray",
                       mols.cols = c("#A50104", "#324376", "#0E7C7B")) +
        theme(legend.title = element_blank()) +
        ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
        scale_y_reverse() + scale_x_reverse())
    
    slide_obj@meta.data$tumor_labels <- NA
    slide_obj@meta.data$tumor_labels <- pData(all_slides)[paste0(row.names(slide_obj@meta.data), 
                                                                 "_S", slide_num),]$tumor_labels
    
    print(ImageDimPlot(slide_obj,
                       fov = paste0("FOV", fov_num), #molecules = c("CEACAM6", "KRT17"), nmols = 100000,
                       group.by = "tumor_labels",
                       cols = temp_cols,
                       na.value = "lightgrey",
                       alpha = 0.8,
                       border.size = 0.1,
                       border.color = "darkgray",
                       dark.background = FALSE) +
        theme(legend.title = element_blank()) +
        ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
        scale_y_reverse() + scale_x_reverse())

}


slide1 <- readRDS(paste0(int_path, "slide1.RDS"))
slide_obj=slide1
fov_num=17
slide_num=1
slide2 <- readRDS(paste0(int_path, "slide2.RDS"))
slide6 <- readRDS(paste0(int_path, "slide6.RDS"))



pdf(paste0(home_path, "man_plots/fibhigh_tumors.pdf"), 
    width = 9, height = 11)

make_fig3_fib_plots(1, 18, slide1)
make_fig3_fib_plots(1, 17, slide1)
make_fig3_fib_plots(1, 23, slide1)
make_fig3_fib_plots(1, 13, slide1)
make_fig3_fib_plots(1, 14, slide1)
make_fig3_fib_plots(1, 22, slide1)
make_fig3_fib_plots(2, 12, slide2)
make_fig3_fib_plots(6, 15, slide6)

dev.off()
