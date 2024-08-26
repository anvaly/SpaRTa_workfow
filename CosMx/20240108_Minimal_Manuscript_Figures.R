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

tumor_cols <- list('Ductal (6)' = "#4CB0E1",
                   'Tumor Classical (1,2,3,5,7,8,11,14)' = "#F16666",
                   'Tumor Classical / Proliferative (9)' = "#9A2626",
                   'Tumor Classical / Mixed (11)' = "#96257D",
                   'Tumor Mixed (10)' = '#CA3FAE',
                   'Tumor Mixed / Fibro adjacent (12)' = "#2C9FAF",
                   'Tumor Basal (4,15)' = "#681D59",
                   'Tumor Other (13)' = "lightgrey")

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

imm_cols <- list('Plasma (7,9,17)' = "#FEDC4B",
                 'B cell (16)' = "#FFCE07",
                 '? (23)' = "gray",
                 'Macrophage (1,2,5,8,10,11,13,18,20,21,22)' = "#FF9F2C",
                 'Dendritic (14)' = "#9B5E2C",
                 'T cell (3,4,6,15,19)' = "#D66100",
                 'Mast (12)' = "#E0B400")


all_sig <- read.table(paste0(home_path, "data/signatures.txt"), 
                      header = TRUE, sep = "\t")
all_sig <- tidyr::separate_rows(all_sig, Gene_list, sep = ",")

all_slides <- readRDS(file = paste0(int_path, "all_slides_int.RDS"))
all_slides@clusters$UMAP <- all_slides@clusters$Aligned

# exclude 8 cells in count because they're in two failed FOVs (one with 2 cells (S2_P1 FOV22), one with 6 cells (S2_P1 FOV23))
all_slides <- all_slides[,!(pData(all_slides)$tissue == "S2_P1" & pData(all_slides)$fov %in% c(22, 23))]
# cell count
nrow(pData(all_slides))

temp <- as.data.frame(table(pData(all_slides)$tissue, pData(all_slides)$fov))
# FOV count
sum(temp$Freq > 0)

sum(pData(all_slides)$comb_labels_clusts == "Tumor (4,6,8,9,10,15,17,19,26)")/nrow(pData(all_slides))

slide6 <- readRDS(file = paste0(int_path, "slide6.RDS"))
slide6@meta.data$cell_type_cluster <- NA
slide6@meta.data$cell_type_cluster <- pData(all_slides)[paste0(row.names(slide6@meta.data), "_S6"),]$comb_labels_clusts

png(paste0(home_path, "man_plots/F1B_S6_P1_all_fov_cell_type.png"), width = 16, height = 13, res = 300 ,units = "in")
ImageDimPlot(slide6, fov = "P1", group.by = "cell_type_cluster", cols = all_cols,
             alpha = 0.8, border.size = 0.1, border.color = "darkgray", dark.background = FALSE) + 
  theme(legend.title = element_blank()) +
  ggtitle(paste0("Slide 6, Position 1")) + scale_y_reverse() + scale_x_reverse()
dev.off()

seur_plot <- as.Seurat(all_slides, assay = NULL)

png(paste0(home_path, "man_plots/F1C_all_UMAP_seur_labels_wide.png"), width = 9, 
    height = 5, res = 300, units = "in")
Seurat::DimPlot(object = seur_plot, 
                group.by = "comb_labels_clusts", 
                pt.size = 0.1,
                reduction = "UMAP",
                cols = all_cols, 
                raster = FALSE) + 
  ggplot2::labs(x="Dim 1", y="Dim 2", title="UMAP") + 
  ggplot2::theme(plot.title = element_text(size=6+2),
                 legend.title = element_blank(), 
                 legend.text = element_text(size=6, hjust = 1, angle = 0),
                 legend.key.size = unit(0.25, "line"), 
                 legend.position = "right")
dev.off()

marker_list <- read.table("/path/to/biomarkers.txt",
                          header = T)
seur_plot$comb_labels_clusts2 <- factor(seur_plot$comb_labels_clusts, 
                                        levels = rev(c("Acinar (16,21)", "Endocrine (23)",
                                                       "Ductal (12)", "Tumor (4,6,8,9,10,15,17,19,26)",
                                                       "Fibroblast (1,3,20)", "Stellate (18)", 
                                                       "Endothelial (14,25)", "Macrophage (2,11,28)", 
                                                       "Mast (22)", "B cell/Plasma (13,24)",
                                                       "T cell (5)", "Dendritic (27)",
                                                       "Low quality (7)")))


png(paste0(home_path, "man_plots/F1D_full_umap_dotplot_aligned_markers.png"), width = 13, height = 4, 
    res = 300, units = "in")
DotPlot(object = seur_plot, 
        features = marker_list$Symbol, 
        group.by = "comb_labels_clusts2") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  theme(legend.title = element_text(size = 10), 
              legend.text  = element_text(size = 10),
              legend.key.size = unit(0.7, "lines"))#+
  #scale_color_gradient(low = "lightgrey", high = "blue", oob = scales::squish, limits = c(0, 1))
dev.off()

fibro <- readRDS(file = paste0(int_path, "all_slides_fibro_integrated.RDS"))
fibro_plot <- as.Seurat(fibro, assay = NULL)

png(paste0(home_path, "man_plots/F1S2D_stroma_UMAP_seur_labels.png"), width = 8, 
    height = 5, res = 300, units = "in")
Seurat::DimPlot(object = fibro_plot, 
                group.by = "comb_labels_clusts", 
                pt.size = 0.1,
                reduction = "UMAP",
                cols = stroma_cols, 
                raster = FALSE) + 
  ggplot2::labs(x="Dim 1", y="Dim 2", title="UMAP") + 
  ggplot2::theme(plot.title = element_text(size=6+2),
                 legend.title = element_blank(), 
                 legend.text = element_text(size=6, hjust = 1, angle = 0),
                 legend.key.size = unit(0.25, "line"), 
                 legend.position = "right")
dev.off()

table(pData(fibro)$comb_labels_clusts)

imm <- readRDS(file = paste0(int_path, "all_slides_immune_integrated.RDS")) 
imm_plot <- as.Seurat(imm, assay = NULL)

png(paste0(home_path, "man_plots/F1S2E_immune_UMAP_seur_labels.png"), width = 8, 
    height = 5, res = 300, units = "in")
Seurat::DimPlot(object = imm_plot, 
                group.by = "comb_labels_clusts", 
                pt.size = 0.1,
                reduction = "UMAP",
                cols = imm_cols, 
                raster = FALSE) + 
  ggplot2::labs(x="Dim 1", y="Dim 2", title="UMAP") + 
  ggplot2::theme(plot.title = element_text(size=6+2),
                 legend.title = element_blank(), 
                 legend.text = element_text(size=6, hjust = 1, angle = 0),
                 legend.key.size = unit(0.25, "line"), 
                 legend.position = "right")
dev.off()

table(pData(imm)$comb_labels_clusts)

tumor <- readRDS(file = paste0(int_path, "all_slides_tumor_integrated.RDS"))
tumor_plot <- as.Seurat(tumor, assay = NULL)

png(paste0(home_path, "man_plots/F1S2F_tumor_seur_labels.png"), width = 7, 
    height = 5, res = 300, units = "in")
Seurat::DimPlot(object = tumor_plot, 
                group.by = "comb_labels_clusts", 
                pt.size = 0.1,
                reduction = "UMAP",
                cols = tumor_cols, 
                raster = FALSE) + 
  ggplot2::labs(x="Dim 1", y="Dim 2", title="UMAP") + 
  ggplot2::theme(plot.title = element_text(size=6+2),
                 legend.title = element_blank(), 
                 legend.text = element_text(size=6, hjust = 1, angle = 0),
                 legend.key.size = unit(0.25, "line"), 
                 legend.position = "right")
dev.off()

table(pData(tumor)$comb_labels_clusts)

tcell_sub <- readRDS(file = paste0(int_path, "all_slides_tcell_integrated.RDS"))

tcell_plot <- as.Seurat(tcell_sub, assay = NULL)

png(paste0(home_path, "man_plots/F1S2G_tcell_UMAP_seur_labels.png"), width = 5, 
    height = 4, res = 300, units = "in")
Seurat::DimPlot(object = tcell_plot, 
                       group.by = "aligned_clusts", 
                       pt.size = 0.1,
                       reduction = "UMAP",
                       cols = RColorBrewer::brewer.pal(8,"Set2"), 
                       raster = FALSE) + 
    ggplot2::labs(x="Dim 1", y="Dim 2", title="UMAP") + 
    ggplot2::theme(plot.title = element_text(size=6+2),
                   legend.title = element_blank(), 
                   legend.text = element_text(size=6, hjust = 1, angle = 0),
                   legend.key.size = unit(0.25, "line"), 
                   legend.position = "right")
dev.off()

treg <- c("CD4", "FOXP3", "IL2RA", "TIGIT", "ICOS", "CTLA4", "CCR8")
exhausted <- c("HAVCR2", "TOX", "PDCD1", "LAG3", "CTLA4", "TIGIT", "CD38", "ENTPD1")

tcell_plot <- UCell::AddModuleScore_UCell(tcell_plot, features=list(treg), name="treg")
tcell_plot <- UCell::AddModuleScore_UCell(tcell_plot, features=list(exhausted), name="exhausted")

pData(tcell_sub)$treg_ucell <- tcell_plot$signature_1treg
pData(tcell_sub)$exhausted_ucell <- tcell_plot$signature_1exhausted


for_plot <- data.frame(aligned_clusts = pData(tcell_sub)$aligned_clusts,
                       signature_val_treg = pData(tcell_sub)$treg_ucell)

png(paste0(home_path, "man_plots/F1S2G_dotplot_tcell_subset_treg.png"), width = 5.5, height = 3, 
    res = 300, units = "in")
DotPlot(object = tcell_plot, 
        features = treg[treg %in% row.names(fData(tcell_sub))], 
        group.by = "aligned_clusts") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  theme(legend.title = element_text(size = 10), 
              legend.text  = element_text(size = 10),
              legend.key.size = unit(0.7, "lines")) + coord_flip() 
dev.off()

png(paste0(home_path, "man_plots/tcell_UMAP_treg_sig.png"), width = 2, 
    height = 2, res = 300, units = "in")
plot_cells(tcell_sub, color_cells_by = "treg_ucell", 
                 cell_size = 0.1, alpha = 0.4) + 
  ggtitle("Treg Signature") + theme_classic(base_size = 6) + 
  NoLegend() + scale_color_distiller(palette = "OrRd", direction = 1)
dev.off()

sum(pData(tcell_sub)$aligned_clusts == 6)/nrow(pData(tcell_sub))

pData(all_slides)$patient_label <- factor(pData(all_slides)$tissue)
temp <- setNames(names(patient_labels), patient_labels)
pData(all_slides)$patient_label <- forcats::fct_recode(pData(all_slides)$patient_label, !!!temp)
for_plot <- pData(all_slides)[,c("patient_label", "comb_labels")]
for_plot$comb_labels <- as.character(for_plot$comb_labels)
for_plot$comb_labels[for_plot$comb_labels %in% c("Acinar", "Endocrine")] <- "Acinar / Endocrine"
for_plot$comb_labels[for_plot$comb_labels %in% c("Fibroblast", "Stellate", "Endothelial")] <- "Stroma"
for_plot$comb_labels[for_plot$comb_labels %in% c("B cell/Plasma", "T cell", "Macrophage", "Dendritic", "Mast")] <- "Immune"
for_plot$comb_labels[for_plot$comb_labels %in% c("Low quality")] <- "Other"
for_plot$comb_labels <- factor(for_plot$comb_labels, levels = c("Acinar / Endocrine", "Ductal", "Tumor", "Stroma", "Immune", "Other"))

for_plot <- as.data.frame(table(for_plot[,c("patient_label", "comb_labels")]))

for_plot$patient_label <- factor(for_plot$patient_label, levels = c("C3_D12", "C3_D5", "C3_D6", "C3_D9", "C3_D4", "C2_D1", "C2_D8", "C2_D10", "C2_D6"))

temp_cols <- list("Acinar / Endocrine" = "#153C65",
                  "Ductal" = "#4CB0E1", 
                  "Tumor" = "#9A2626",
                  "Stroma" = "#E5E5AA", 
                  "Immune" = "#DE7126", 
                  "Other" = "lightgray")


png(paste0(home_path, "man_plots/F1S2H_smi_major_types_freq_bar.png"), width = 3, height = 2, 
    res = 300, units = "in")
ggplot(data = for_plot, aes(patient_label, Freq, fill = comb_labels)) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw(base_size = 6) + 
  labs(x = "", y = "Fraction of cells") + 
  monocle3:::monocle_theme_opts() + 
  theme(legend.key.size = unit(5, "point")) + 
  guides(fill = guide_legend(override.aes = list(size=2))) +
  scale_fill_manual(values = temp_cols, name = "Cell type") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

#SingleCell_FOVCount
pData(all_slides)$slide_fov <- paste0(pData(all_slides)$tissue, "_", pData(all_slides)$fov)
as.data.frame(table(unique(pData(all_slides)[,c("patient_label", "slide_fov")])$patient_label))
#SingleCell_Num
as.data.frame(table(pData(all_slides)[,c("patient_label")]))

pData(all_slides)$patient_ROI_label <- factor(pData(all_slides)$tissue)
temp <- setNames(names(patient_ROI_labels), patient_ROI_labels)
pData(all_slides)$patient_ROI_label <- forcats::fct_recode(pData(all_slides)$patient_ROI_label, !!!temp)

pData(tumor)$comb_labels_clusts <- as.character(pData(tumor)$comb_labels_clusts)
pData(all_slides)$tumor_comb_labels <- ""

pData(all_slides)$tumor_comb_labels <- pData(tumor)[row.names(pData(all_slides)),]$comb_labels_clusts

pData(all_slides)$tumor_labels <- ifelse(!is.na(pData(all_slides)$tumor_comb_labels), as.character(pData(all_slides)$tumor_comb_labels), 
                                         as.character(pData(all_slides)$comb_labels_clusts))


t1 <- as.data.frame(pData(all_slides)) %>% group_by(patient_ROI_label, fov) %>% count(tumor_labels) %>% 
  tidyr::pivot_wider(names_from = tumor_labels, values_from = n, values_fill = 0)


t2 <- as.data.frame(pData(all_slides)) %>% group_by(patient_ROI_label, fov) %>% 
  summarize(CellCount = n(), TotalCount = sum(nCount_Nanostring))

t1.2 <- merge(t2, t1, by = c("patient_ROI_label", "fov"))

names(t1.2) <- c("patient_ROI_label", "FOV", "NumCells", "TotalMolCount", "B cell/Plasma_Count", 
                 "Ductal-like_Count", "Endocrine_Count", "Endothelial_Count", "Fibroblast_Count", 
                 "Other_Count", "Macrophage_Count", "Mast_Count", "Stellate_Count", "T cell_Count", 
                 "Basal PDAC_Count", "Classical PDAC_Count", "Proliferative PDAC_Count", 
                 "Mixed PDAC_Count", "Fibroblast high PDAC_Count", "Acinar_Count", 
                 "Dendritic_Count", "Tumor Other_Count")

write.csv(t1.2, paste0(int_path, "TS1_3_singlecell_fov_summary.csv"),
          quote = FALSE, row.names = FALSE)

all_sig_cx <- all_sig[all_sig$Gene_list %in% row.names(fData(all_slides)),]
cx_sig_collapse <- all_sig_cx %>% group_by(Signature_name) %>% 
  summarize(sig_list = paste(Gene_list, collapse = ","))

write.table(cx_sig_collapse, paste0(int_path, "TS1_4_cosmx_sigs.csv"),
          quote = FALSE, row.names = FALSE, col.names = TRUE)

temp_cols <- c(all_cols, imm_cols)

options(future.globals.maxSize = 6.291e+8) # 600MiB
slide4 <- readRDS(paste0(int_path, "slide4.RDS"))

slide_num <- 4
fov_num <- 11
slide_obj <- slide4
fov_met <- slide_obj@meta.data[slide_obj@meta.data$fov == fov_num,]
fov.crop <- Crop(slide_obj[[paste0("Slide", slide_num)]], 
                 y = c(min(fov_met$CenterX_global_px),
                       max(fov_met$CenterX_global_px)), 
                 x = c(min(fov_met$CenterY_global_px), 
                       max(fov_met$CenterY_global_px)))
slide_obj[[paste0("FOV", fov_num)]] <- fov.crop
DefaultBoundary(slide_obj[[paste0("FOV", fov_num)]]) <- "segmentation"

#slide_obj@meta.data$cell_subtype <- NA
#slide_obj@meta.data$cell_subtype <- pData(tumor)[paste0(row.names(slide_obj@meta.data), 
#                                                               "_S", slide_num),]$cell_subtype

png(paste0(home_path, "man_plots/F2E_S", slide_num, "_FOV", fov_num, "_cell_subtype_LTF_CCL19.png"), 
    width = 9, height = 11, res = 300, units = "in")
print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), molecules = c("LTF", "CCL19"), nmols = 100000,
                   alpha = 0.6, border.size = 0.1, border.color = "darkgray", cols = "lightgray", mols.cols = c("#17C3B2", "#EF2D56")) +
        theme(legend.title = element_blank()) +
        ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
        scale_y_reverse() + scale_x_reverse())
dev.off()


slide_obj@meta.data$comb_labels_clusts <- NA
slide_obj@meta.data$comb_labels_clusts <- as.character(pData(all_slides)[paste0(row.names(slide_obj@meta.data), 
                                                                   "_S", slide_num),]$comb_labels_clusts)

slide_obj@meta.data$comb_labels_clusts2 <- NA
slide_obj@meta.data$comb_labels_clusts2 <- as.character(pData(imm)[paste0(row.names(slide_obj@meta.data), 
                                                                   "_S", slide_num),]$comb_labels_clusts)
slide_obj@meta.data$comb_labels_clusts2 <- ifelse(!is.na(slide_obj@meta.data$comb_labels_clusts2), 
                                                  slide_obj@meta.data$comb_labels_clusts2, 
                                                  slide_obj@meta.data$comb_labels_clusts)

png(paste0(home_path, "man_plots/F2E_S", slide_num, "_FOV", fov_num, "_cell_type.png"), 
    width = 9, height = 11, res = 300, units = "in")
print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), #molecules = c("CEACAM6", "KRT17"), nmols = 100000,
                   group.by = "comb_labels_clusts2", cols = temp_cols ,na.value = "lightgrey",
                   alpha = 0.6, border.size = 0.1, border.color = "darkgray") +
        theme(legend.title = element_blank()) +
        ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
        scale_y_reverse() + scale_x_reverse())
dev.off()


slide_num <- 4
fov_num <- 2
slide_obj <- slide4
fov_met <- slide_obj@meta.data[slide_obj@meta.data$fov == fov_num,]
fov.crop <- Crop(slide_obj[[paste0("Slide", slide_num)]], 
                 y = c(min(fov_met$CenterX_global_px),
                       max(fov_met$CenterX_global_px)), 
                 x = c(min(fov_met$CenterY_global_px), 
                       max(fov_met$CenterY_global_px)))
slide_obj[[paste0("FOV", fov_num)]] <- fov.crop
DefaultBoundary(slide_obj[[paste0("FOV", fov_num)]]) <- "segmentation"

png(paste0(home_path, "man_plots/F2E_S", slide_num, "_FOV", fov_num, "_cell_subtype_LTF_CCL19.png"), 
    width = 9, height = 11, res = 300, units = "in")
print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), molecules = c("LTF", "CCL19"), nmols = 100000,
                   alpha = 0.6, border.size = 0.1, border.color = "darkgray", cols = "lightgray", mols.cols = c("#17C3B2", "#EF2D56")) +
        theme(legend.title = element_blank()) +
        ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
        scale_y_reverse() + scale_x_reverse())
dev.off()


slide_obj@meta.data$comb_labels_clusts <- NA
slide_obj@meta.data$comb_labels_clusts <- as.character(pData(all_slides)[paste0(row.names(slide_obj@meta.data), 
                                                                   "_S", slide_num),]$comb_labels_clusts)

slide_obj@meta.data$comb_labels_clusts2 <- NA
slide_obj@meta.data$comb_labels_clusts2 <- as.character(pData(imm)[paste0(row.names(slide_obj@meta.data), 
                                                                   "_S", slide_num),]$comb_labels_clusts)
slide_obj@meta.data$comb_labels_clusts2 <- ifelse(!is.na(slide_obj@meta.data$comb_labels_clusts2), 
                                                  slide_obj@meta.data$comb_labels_clusts2, 
                                                  slide_obj@meta.data$comb_labels_clusts)
png(paste0(home_path, "man_plots/F2E_S", slide_num, "_FOV", fov_num, "_cell_type.png"), 
    width = 9, height = 11, res = 300, units = "in")
print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), #molecules = c("CEACAM6", "KRT17"), nmols = 100000,
                   group.by = "comb_labels_clusts2", cols = temp_cols,na.value = "lightgrey",
                   alpha = 0.6, border.size = 0.1, border.color = "darkgray") +
        theme(legend.title = element_blank()) +
        ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
        scale_y_reverse() + scale_x_reverse())
dev.off()

fibro_plot <- UCell::AddModuleScore_UCell(fibro_plot, features=list(all_sig[all_sig$Signature_name == "selected",]$Gene_list), name="selected")

pData(fibro)$ifna_sig <- fibro_plot@meta.data[row.names(pData(fibro)),]$signature_1BMS.Pathway.IFNa
for_plot <- data.frame(comb_labels_clusts = pData(fibro)$comb_labels_clusts,
                       signature_val = pData(fibro)$ifna_sig)

png(paste0(home_path, "man_plots/F2S1I_IFNa_signature_boxplot_fibro_labels.png"), width = 2, height = 2.5, 
    res = 300, units = "in")
ggplot(for_plot, aes(comb_labels_clusts, signature_val, fill = comb_labels_clusts)) + 
  geom_boxplot(size = 0.2, outlier.size = 0.01) +
  scale_fill_manual(values=stroma_cols) + 
  theme_bw(base_size = 6) + monocle3:::monocle_theme_opts() + 
  theme(legend.title = element_blank()) + 
  labs(x = "", y = "UCell signature value") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position = "none") 
dev.off()

png(paste0(home_path, "man_plots/F2S1_stroma_dotplot_IFNa_markers.png"), width = 4, height = 5.5, 
    res = 300, units = "in")
DotPlot(object = fibro_plot, 
        features = intersect(all_sig[all_sig$Signature_name == "selected",]$Gene_list, row.names(fData(fibro))), 
        group.by = "comb_labels_clusts") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  theme(legend.title = element_text(size = 10), 
              legend.text  = element_text(size = 10),
              legend.key.size = unit(0.7, "lines")) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+
  #scale_color_gradient(low = "lightgrey", high = "blue", oob = scales::squish, limits = c(0, 1))
dev.off()


pData(fibro)$slide_fov <- paste0(pData(fibro)$tissue, "_", pData(fibro)$fov)
ifna_list <- as.data.frame(table(pData(fibro)[,c("slide_fov", "comb_labels_clusts")]))
ifna_list <- ifna_list[ifna_list$comb_labels_clusts == "Fibroblasts/IFNa pathway/ resistance (11)",]
ifna_list <- ifna_list[order(ifna_list$Freq, decreasing = T),]

make_ifna_plots <- function(slide_num, fov_num, slide_obj) {
  fov_met <- slide_obj@meta.data[slide_obj@meta.data$fov == fov_num,]
  fov.crop <- Crop(slide_obj[[paste0("Slide", slide_num)]], 
                   y = c(min(fov_met$CenterX_global_px),
                         max(fov_met$CenterX_global_px)), 
                   x = c(min(fov_met$CenterY_global_px), 
                         max(fov_met$CenterY_global_px)))
  slide_obj[[paste0("FOV", fov_num)]] <- fov.crop
  DefaultBoundary(slide_obj[[paste0("FOV", fov_num)]]) <- "segmentation"
  
  png(paste0(home_path, "man_plots/F2S1_S", slide_num, "_FOV", fov_num, "_cell_subtype_OAS3_MX1_IFIT1_IFI27_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), 
                     molecules = c("OAS3", "MX1", "IFIT1", "IFI27"), 
                     nmols = 100000, dark.background = FALSE,
                     alpha = 0.6, border.size = 0.1, border.color = "darkgray", 
                     cols = "lightgray", 
                     mols.cols = c("#A50104", "#324376", "#3F7D20", "#ED9B40")) +
          theme(legend.title = element_blank()) +
          ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
          scale_y_reverse() + scale_x_reverse())
  dev.off()
  
  slide_obj@meta.data$comb_labels_clusts <- NA
  slide_obj@meta.data$comb_labels_clusts <- pData(fibro)[paste0(row.names(slide_obj@meta.data), 
                                                                "_S", slide_num),]$comb_labels_clusts
  
  png(paste0(home_path, "man_plots/F2S1_S", slide_num, "_FOV", fov_num, 
             "_fibro_subtype_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), 
                     group.by = "comb_labels_clusts", cols = stroma_cols,
                     na.value = "lightgrey", alpha = 0.8, border.size = 0.1, 
                     border.color = "darkgray", dark.background = FALSE) +
          theme(legend.title = element_blank()) +
          ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
          scale_y_reverse() + scale_x_reverse())
  dev.off()
  
  slide_obj@meta.data$ifna_sig <- NA
  slide_obj@meta.data$ifna_sig <- pData(fibro)[paste0(row.names(slide_obj@meta.data),
                                                       "_S", slide_num),]$ifna_sig
  
  slide_obj$ifna_sig[is.na(slide_obj$comb_labels_clusts)] <- NA

  png(paste0(home_path, "man_plots/F2S1_S", slide_num, "_FOV", fov_num, 
             "_fibro_ifna_sig_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageFeaturePlot(slide_obj, fov = paste0("FOV", fov_num), 
                         features = "ifna_sig", alpha = 0.8, border.size = 0.1, 
                         border.color = "darkgray", dark.background = FALSE) +
          theme(legend.title = element_blank()) +
          ggtitle(paste0("IFNA")) + scale_y_reverse() + scale_x_reverse() + 
          scale_fill_gradient(low = "#000000", high = '#FFFF00', 
                              na.value = "lightgrey"))
  dev.off()
}

slide2 <- readRDS(paste0(int_path, "slide2.RDS"))
make_ifna_plots(6, 7, slide6) #S6_P1, C3_D5
#make_ifna_plots(6, 6, slide6) #S6_P1, C3_D5
make_ifna_plots(2, 17, slide2) #S2_P1 C2_D10

tumor_plot <- as.Seurat(tumor, assay = NULL)

tumor_plot <- tumor_plot[,tumor_plot$comb_labels_clusts != "Tumor Other (13)"]
tumor_plot <- Seurat::SCTransform(tumor_plot, assay="originalexp", variable.features.n=988,
                                  return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")
tumor_plot <- ScaleData(object = tumor_plot, assay = "originalexp")

# Downsample to 30k total cells - Seurat limitation
ds <- lapply(unique(tumor_plot$comb_labels_clusts), function(lab) {
  temp <- tumor_plot@meta.data[tumor_plot$comb_labels_clusts == lab,]
  num <- 5000
  if(nrow(temp) > num)
    return(row.names(temp[sample(1:nrow(temp), num, replace = FALSE),]))
  return(row.names(temp))
})

sig_genes <- read.table("/path/to/genelist.txt", header = TRUE)
sig_genes <- c("CLU", "CRP",
               "OLFM4", "DMBT1", "CEACAM1",
               "S100P", "CEACAM5", "CEACAM6", "AGR2", "TFF1", "S100A6", "KRT8", "KRT18",
               "UBE2C", "CENPF", "TOP2A", "MKI67",
               "KRT13", "KRT23",
               "LY6D", "CCL20", "CD9", "TNFRSF11B", "KRT80", "KRT7", "KRT19", "KRT17", "KRT6B", "KRT15", "AREG", "KRT16", "SLC2A1", "NDRG1", "KRT5", "KRT6A", "KRT14",
               "COL1A1", "COL1A2", "FN1")
#sig_genes <- c("GC", "REG1A", "SOD2", "FGFR3", "CCL28", "ANXA4", "CLU", "LEFTY1", "MMP7", "SPP1", "CRP",
#               "AGR2", "CEACAM5", "CEACAM6", "ERBB3", "LGALS4", "S100P", "TFF1", "TFF3", "TSPAN8",
#               "UBE2C", "CENPF", "TOP2A", "MKI67",
#               "KRT13", "KRT23",
#               "AREG", "KRT15", "KRT17", "KRT6A", "KRT6C", "KRT7", "LY6D", "S100A2", "SLC2A1", "TNS4", "KRT16",
#               "COL1A1", "COL1A2", "FN1")

tumor_plot_ds <- tumor_plot[,unlist(ds)]
tumor_plot_ds$comb_labels_clusts <- factor(tumor_plot_ds$comb_labels_clusts, levels = c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"))

png(paste0(home_path, "man_plots/F3B_tumor_heatmap_key_markers.png"), width = 7, height = 10, 
    res = 300, units = "in")
DoHeatmap(object = tumor_plot_ds, 
          features = sig_genes, 
          group.by = "comb_labels_clusts", raster = TRUE, group.colors = unlist(tumor_cols), 
          disp.max = 2, lines.width = 100, draw.lines = TRUE) +
  scale_fill_gradient2(low = '#FF00FF', mid = "#000000", high = '#FFFF00', limits = c(-2, 2)) +
  scale_color_manual(name = "Identity", values = tumor_cols)
dev.off()

pData(tumor)$slide_fov <- paste0(pData(tumor)$tissue, "_", pData(tumor)$fov)
mixed_list <- as.data.frame(table(pData(tumor)[,c("slide_fov", "comb_labels_clusts")]))
mixed_list <- mixed_list[mixed_list$comb_labels_clusts == "Tumor Mixed (10)",]
mixed_list <- mixed_list[order(mixed_list$Freq, decreasing = T),]

tumor_seur <- as.Seurat(tumor, assay = NULL)
markers_l <- all_sig[all_sig$Signature_name == "selected",]$Gene_list
tumor_seur <- UCell::AddModuleScore_UCell(tumor_seur, features=list(markers_l[markers_l %in% row.names(fData(all_slides))]), name="selected")
pData(tumor)$ucell_seurat_prolif_sig <- tumor_seur@meta.data[row.names(tumor_seur@meta.data),]$signature_1Seurat.Proliferation

make_fig3_mixed_plots <- function(slide_num, fov_num, slide_obj, 
                                  low_ntile = 0.1, high_ntile = 0.9) {
  fov_met <- slide_obj@meta.data[slide_obj@meta.data$fov == fov_num,]
  fov.crop <- Crop(slide_obj[[paste0("Slide", slide_num)]], 
                   y = c(min(fov_met$CenterX_global_px),
                         max(fov_met$CenterX_global_px)), 
                   x = c(min(fov_met$CenterY_global_px), 
                         max(fov_met$CenterY_global_px)))
  slide_obj[[paste0("FOV", fov_num)]] <- fov.crop
  DefaultBoundary(slide_obj[[paste0("FOV", fov_num)]]) <- "segmentation"
  
  png(paste0(home_path, "man_plots/F3_S", slide_num, "_FOV", fov_num, 
             "_cell_subtype_CEACAM6_KRT6A_MKI67_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), 
                     molecules = c("CEACAM6", "KRT6A", "MKI67"), nmols = 100000, 
                     dark.background = FALSE, alpha = 0.6, border.size = 0.1, 
                     border.color = "darkgray", cols = "lightgray", 
                     mols.cols = c("#A50104", "#324376", "#3F7D20")) +
          theme(legend.title = element_blank()) +
          ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
          scale_y_reverse() + scale_x_reverse())
  dev.off()
  
    png(paste0(home_path, "man_plots/F3_S", slide_num, "_FOV", fov_num, 
               "_cell_subtype_CEACAM6_KRT6A_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), 
                     molecules = c("CEACAM6", "KRT6A"), nmols = 100000, 
                     dark.background = FALSE, alpha = 0.6, border.size = 0.1, 
                     border.color = "darkgray", cols = "lightgray", 
                     mols.cols = c("#A50104", "#324376")) +
          theme(legend.title = element_blank()) +
          ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
          scale_y_reverse() + scale_x_reverse())
  dev.off()
  
  slide_obj@meta.data$comb_labels_clusts <- NA
  slide_obj@meta.data$comb_labels_clusts <- pData(tumor)[paste0(row.names(slide_obj@meta.data), 
                                                                "_S", slide_num),]$comb_labels_clusts
  
  png(paste0(home_path, "man_plots/F3_S", slide_num, "_FOV", fov_num, 
             "_tumor_subtype_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), 
                     group.by = "comb_labels_clusts", cols = tumor_cols,
                     na.value = "lightgrey", alpha = 0.8, border.size = 0.1, 
                     border.color = "darkgray", dark.background = FALSE) +
          theme(legend.title = element_blank()) +
          ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
          scale_y_reverse() + scale_x_reverse())
  dev.off()
  
  slide_obj@meta.data$classical_sig <- NA
  slide_obj@meta.data$classical_sig <- pData(tumor)[paste0(row.names(slide_obj@meta.data),
                                                           "_S", 
                                                           slide_num),]$ucell_collisson_classical_sig
  slide_obj@meta.data$basal_sig <- NA
  slide_obj@meta.data$basal_sig <- pData(tumor)[paste0(row.names(slide_obj@meta.data),
                                                       "_S", slide_num),]$ucell_moffitt_basal_sig
  slide_obj@meta.data$prol_sig <- NA
  slide_obj@meta.data$prol_sig <- pData(tumor)[paste0(row.names(slide_obj@meta.data),
                                                       "_S", slide_num),]$ucell_seurat_prolif_sig
  
  slide_obj$classical_sig[is.na(slide_obj$comb_labels_clusts)] <- NA
  slide_obj$basal_sig[is.na(slide_obj$comb_labels_clusts)] <- NA
  slide_obj$prol_sig[is.na(slide_obj$comb_labels_clusts)] <- NA

  bas_low <- quantile(pData(tumor)$ucell_moffitt_basal_sig, 
                      c(low_ntile, high_ntile))[[1]]
  bas_high <- quantile(pData(tumor)$ucell_moffitt_basal_sig, 
                       c(low_ntile, high_ntile))[[2]]
  
  slide_obj$basal_sig[slide_obj$basal_sig < bas_low] <- bas_low
  slide_obj$basal_sig[slide_obj$basal_sig > bas_high] <- bas_high
  png(paste0(home_path, "man_plots/F3_S", slide_num, "_FOV", fov_num, 
             "_tumor_basal_sig_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageFeaturePlot(slide_obj, fov = paste0("FOV", fov_num), 
                         features = "basal_sig", alpha = 0.8, border.size = 0.1, 
                         border.color = "darkgray", dark.background = FALSE) +
          theme(legend.title = element_blank()) + 
          ggtitle(paste0("Moffitt Basal")) + scale_y_reverse() + 
          scale_x_reverse() + 
          scale_fill_gradient(low = "#000000", high = '#FFFF00', 
                              na.value = "lightgrey", 
                              limits = c(bas_low,bas_high)))
  dev.off()
  
  class_low <- quantile(pData(tumor)$ucell_collisson_classical_sig, 
                        c(low_ntile, high_ntile))[[1]]
  class_high <- quantile(pData(tumor)$ucell_collisson_classical_sig, 
                         c(low_ntile, high_ntile))[[2]]
  
  slide_obj$classical_sig[slide_obj$classical_sig < class_low] <- class_low
  slide_obj$classical_sig[slide_obj$classical_sig > class_high] <- class_high
  png(paste0(home_path, "man_plots/F3_S", slide_num, "_FOV", fov_num, 
             "_tumor_classical_sig_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageFeaturePlot(slide_obj, fov = paste0("FOV", fov_num), 
                         features = "classical_sig", alpha = 0.8, 
                         border.size = 0.1, border.color = "darkgray", 
                         dark.background = FALSE) +
          theme(legend.title = element_blank()) +
          ggtitle(paste0("Collisson Classical")) + scale_y_reverse() + 
          scale_x_reverse() + 
          scale_fill_gradient(low = "#000000", high = '#FFFF00', 
                              na.value = "lightgrey", 
                              limits = c(class_low,class_high)))
  dev.off()
  
  prof_low <- quantile(pData(tumor)$ucell_seurat_prolif_sig, 
                       c(low_ntile, high_ntile))[[1]]
  prof_high <- quantile(pData(tumor)$ucell_seurat_prolif_sig, 
                        c(low_ntile, high_ntile))[[2]]
  
  slide_obj$prol_sig[slide_obj$prol_sig < prof_low] <- prof_low
  slide_obj$prol_sig[slide_obj$prol_sig > prof_high] <- prof_high
  png(paste0(home_path, "man_plots/F3_S", slide_num, "_FOV", fov_num, 
             "_tumor_proliferation_sig_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageFeaturePlot(slide_obj, fov = paste0("FOV", fov_num), 
                         features = "prol_sig", alpha = 0.8, border.size = 0.1, 
                         border.color = "darkgray", dark.background = FALSE) +
          theme(legend.title = element_blank()) +
          ggtitle(paste0("Proliferative")) + scale_y_reverse() + 
          scale_x_reverse() + 
          scale_fill_gradient(low = "#000000", high = '#FFFF00', 
                              na.value = "lightgrey", 
                              limits = c(prof_low,prof_high)))
  dev.off()
}

options(future.globals.maxSize = 6.291e+8) # 600MiB

slide1 <- readRDS(paste0(int_path, "slide1.RDS"))
slide3 <- readRDS(paste0(int_path, "slide3.RDS"))
slide4 <- readRDS(paste0(int_path, "slide4.RDS"))
#slide5b2 <- readRDS(paste0(int_path, "slide5b2.RDS"))
slide6 <- readRDS(paste0(int_path, "slide6.RDS"))

# Main
make_fig3_mixed_plots(3, 8, slide3, low_ntile = 0, high_ntile = 1) #S3_P1, C3_D4

# Supp
make_fig3_mixed_plots(1, 16, slide1, low_ntile = 0, high_ntile = 1) #S1_P1 C2_D6
make_fig3_mixed_plots(3, 20, slide3, low_ntile = 0, high_ntile = 1) #S3_P2 C3_D4
make_fig3_mixed_plots(4, 19, slide4, low_ntile = 0, high_ntile = 1) #S4_P2 C3_D12
make_fig3_mixed_plots(6, 3, slide6, low_ntile = 0, high_ntile = 1) #S6_P1 C3_D5

#make_fig3_mixed_plots(1, 6, slide1) #S1_P2 C2_D1
#make_fig3_mixed_plots('5b2', 18, slide5b2) #S5b2_P2 C3_D9

tumor_plot <- tumor_seur[,!tumor_seur$comb_labels_clusts %in% c("Tumor Other (13)")]
tumor_plot$comb_labels_clusts <- forcats::fct_rev(tumor_plot$comb_labels_clusts)

png(paste0(home_path, "man_plots/tumor_dotplot_prolif_markers.png"), width = 8, height = 4, 
    res = 300, units = "in")
DotPlot(object = tumor_plot, 
        features = markers_l[markers_l %in% row.names(fData(all_slides))], 
        group.by = "comb_labels_clusts") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))  +
  theme(legend.title = element_text(size = 10), 
              legend.text  = element_text(size = 10),
              legend.key.size = unit(0.7, "lines"))#+
  #scale_color_gradient(low = "lightgrey", high = "blue", oob = scales::squish, limits = c(0, 1))
dev.off()

for_plot <- pData(tumor)[,c("comb_labels_clusts", "ucell_seurat_prolif_sig")]
for_plot <- for_plot[!for_plot$comb_labels_clusts %in% c("Tumor Other (13)"),]
for_plot$comb_labels_clusts <- factor(for_plot$comb_labels_clusts, levels = c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"))

png(paste0(home_path, "man_plots/F3I_tumor_prolif_sig_boxplot.png"), 
        width = 1.5, height = 3, res = 300, units = "in")
ggplot(data = as.data.frame(for_plot), 
       aes(comb_labels_clusts, ucell_seurat_prolif_sig, fill = comb_labels_clusts)) + 
  geom_boxplot(outlier.size = 0.1, size = 0.2) + theme_bw(base_size = 6) + 
  NoLegend() +
  scale_fill_manual(values = tumor_cols) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  labs(x = "", y = "Signature value") 
dev.off()

plot_tme <- function(slide_obj, cell_members, group.by, plot_title = "", 
                     cols = NULL, dark.background = FALSE, buffer = 10, 
                     with_conns = FALSE, edge_df = NULL, slide_suffix = NULL, 
                     linewidth = 0.5) {
  fov_met <- slide_obj@meta.data[cell_members,]
  fov.crop <- Crop(slide_obj[[names(slide_obj@images)[[1]]]], 
                       y = c(min(fov_met$CenterX_global_px) - buffer,
                             max(fov_met$CenterX_global_px) + buffer), 
                       x = c(min(fov_met$CenterY_global_px) - buffer, 
                             max(fov_met$CenterY_global_px) + buffer))
  slide_obj[["for_plot"]] <- fov.crop
  DefaultBoundary(slide_obj[["for_plot"]]) <- "segmentation"
  
  pl <- ImageDimPlot(slide_obj, fov = "for_plot", group.by = group.by, cols = cols,
                     alpha = 0.8, border.size = 0.1, border.color = "darkgray", 
                     dark.background = dark.background) +
            theme(legend.title = element_blank()) +
            ggtitle(plot_title) + 
            scale_y_reverse() + scale_x_reverse()
  
  if (!with_conns) {
    print(pl)
  } else {
    print(pl + ggplot2::geom_segment(data = edge_df[edge_df$from %in% paste0(cell_members, slide_suffix) & 
                                           edge_df$to %in% paste0(cell_members, slide_suffix)], 
                           aes(x = sdimy_begin, y = sdimx_begin,
                                              xend = sdimy_end, yend = sdimx_end),
                           color = "black", linewidth = linewidth, alpha = 0.7, inherit.aes = FALSE))
  }
}

comp_mat <- read.csv(paste0(int_path, "TME_analysis/tme_subtype/composition_full_qc_20.csv"))
meta <- read.csv(paste0(int_path, "TME_analysis/tme_subtype/composition_meta.csv"))
all_edges <- readRDS(file = paste0(int_path, "all_edges.RDS"))

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

pData(tumor)$comb_labels_clusts <- as.character(pData(tumor)$comb_labels_clusts) 
pData(all_slides)$tumor_comb_labels <- NA

pData(all_slides)$tumor_comb_labels <- pData(tumor)[row.names(pData(all_slides)),]$comb_labels_clusts
pData(all_slides)$tumor_labels <- ifelse(!is.na(pData(all_slides)$tumor_comb_labels), as.character(pData(all_slides)$tumor_comb_labels), 
                                         as.character(pData(all_slides)$comb_labels_clusts))

pData(all_slides)$tumor_labels[pData(all_slides)$tumor_labels %in% c("B cell/Plasma (13,24)", "Dendritic (27)", "Macrophage (2,11,28)", "Mast (22)", "T cell (5)")] <- "Immune"

slide5b2 <- readRDS(paste0(int_path, "slide5b2.RDS"))
slide5b2@meta.data$tumor_labels <- NA
slide5b2@meta.data$tumor_labels <- pData(all_slides)[paste0(row.names(slide5b2@meta.data), 
                                                               "_S5b2"),]$tumor_labels
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
temp_cols <- c(temp_cols, all_cols, tumor_cols)

png(paste0(home_path, "man_plots/F3J_empty_tme.png"), width = 8, height = 5, 
    res = 300, units = "in")
plot_tme(slide5b2, gsub("_S5b2", "", comp_sub_long[comp_sub_long$subgraph == "slide5b2_FOV17_1000",]$cell_members), "Slide_name", cols = list("S3" = "lightgray"), with_conns = FALSE, edge_df = all_edges, slide_suffix = "_S5b2", linewidth = 0.8, buffer = 20)
dev.off()

slide5b2$center_cell <- row.names(slide5b2@meta.data) == "1269_17"
png(paste0(home_path, "man_plots/F3J_empty_tme_cc.png"), width = 8, height = 5, 
    res = 300, units = "in")
plot_tme(slide5b2, gsub("_S5b2", "", comp_sub_long[comp_sub_long$subgraph == "slide5b2_FOV17_1000",]$cell_members), "center_cell", cols = list(`FALSE` = "lightgray", `TRUE` = "#186928"), with_conns = FALSE, edge_df = all_edges, slide_suffix = "_S5b2", linewidth = 0.8, buffer = 20)
dev.off()

pdf(paste0(home_path, "man_plots/F3J_empty_tme_with_conns_cc.pdf"), width = 8, height = 5)
plot_tme(slide5b2, gsub("_S5b2", "", comp_sub_long[comp_sub_long$subgraph == "slide5b2_FOV17_1000",]$cell_members), "center_cell", cols = list(`FALSE` = "lightgray", `TRUE` = "#186928"), with_conns = TRUE, edge_df = all_edges, slide_suffix = "_S5b2", linewidth = 0.8, buffer = 20)
dev.off()

png(paste0(home_path, "man_plots/F3J_example_tms_with_conns_labels.png"), width = 8, height = 5, 
    res = 300, units = "in")
plot_tme(slide5b2, gsub("_S5b2", "", comp_sub_long[comp_sub_long$subgraph == "slide5b2_FOV17_1000",]$cell_members), "tumor_labels", cols = temp_cols, with_conns = TRUE, edge_df = all_edges, slide_suffix = "_S5b2", linewidth = 0.8, buffer = 20)
dev.off()

pData(tumor)$comb_labels_clusts <- as.character(pData(tumor)$comb_labels_clusts)
 
comp_sub$center_cell_type <- 
  pData(tumor)[comp_sub$center_cell_id,]$comb_labels_clusts

for_plot <- comp_sub[comp_sub$center_cell_type %in%  c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"),]

for_plot_count <- for_plot %>% group_by(center_cell_type) %>% 
  summarize(High = sum(Fibroblast > 0)/n())

for_plot_count$center_cell_type <- factor(for_plot_count$center_cell_type, levels = c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"))

levels(for_plot_count$center_cell_type) <- c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)")

png(paste0(home_path, "man_plots/F3K_tumor_neighbors_fibro0_bar.png"), width = 2, height = 3, 
    res = 300, units = "in")
ggplot(for_plot_count, aes(center_cell_type, High, fill = center_cell_type)) + 
  geom_bar(stat = "identity") + theme_bw(base_size = 6) + 
  labs(x = "", y = "Fraction of neighborhoods - Fibro > 0") + 
  monocle3:::monocle_theme_opts() + 
  theme(legend.key.size = unit(5, "point")) + 
  guides(fill = guide_legend(override.aes = list(size=2))) +
  scale_fill_manual(values = tumor_cols, name = "Cell type") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  NoLegend()
dev.off()

print(for_plot_count)

slide5 <- readRDS(paste0(int_path, "slide5.RDS"))
slide5b3 <- readRDS(paste0(int_path, "slide5b3.RDS"))
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
temp_cols <- c(temp_cols, all_cols, tumor_cols) 

pData(all_slides)$tumor_comb_labels <- NA
pData(tumor)$comb_labels_clusts <- as.character(pData(tumor)$comb_labels_clusts)
pData(all_slides)$tumor_comb_labels <- pData(tumor)[row.names(pData(all_slides)),]$comb_labels_clusts
pData(all_slides)$tumor_comb_labels <-  as.character(pData(all_slides)$tumor_comb_labels)

pData(all_slides)$tme_type <- as.character(pData(all_slides)$comb_labels)
pData(all_slides)$tme_subtype <- ifelse(is.na(pData(all_slides)$tumor_comb_labels),
                                         as.character(pData(all_slides)$comb_labels),
                                         as.character(pData(all_slides)$tumor_comb_labels))
pData(all_slides)$tme_subtype_ext <- pData(all_slides)$tme_subtype

slide1@meta.data$tme_subtype_ext <- NA
slide1@meta.data$tme_subtype_ext <- pData(all_slides)[paste0(row.names(slide1@meta.data), 
                                                               "_S1"),]$tme_subtype_ext
slide2@meta.data$tme_subtype_ext <- NA
slide2@meta.data$tme_subtype_ext <- pData(all_slides)[paste0(row.names(slide2@meta.data), 
                                                               "_S2"),]$tme_subtype_ext
slide3@meta.data$tme_subtype_ext <- NA
slide3@meta.data$tme_subtype_ext <- pData(all_slides)[paste0(row.names(slide3@meta.data), 
                                                               "_S3"),]$tme_subtype_ext
slide4@meta.data$tme_subtype_ext <- NA
slide4@meta.data$tme_subtype_ext <- pData(all_slides)[paste0(row.names(slide4@meta.data), 
                                                               "_S4"),]$tme_subtype_ext
slide5@meta.data$tme_subtype_ext <- NA
slide5@meta.data$tme_subtype_ext <- pData(all_slides)[paste0(row.names(slide5@meta.data), 
                                                               "_S5"),]$tme_subtype_ext
slide5b2@meta.data$tme_subtype_ext <- NA
slide5b2@meta.data$tme_subtype_ext <- pData(all_slides)[paste0(row.names(slide5b2@meta.data), 
                                                               "_S5b2"),]$tme_subtype_ext
slide5b3@meta.data$tme_subtype_ext <- NA
slide5b3@meta.data$tme_subtype_ext <- pData(all_slides)[paste0(row.names(slide5b3@meta.data), 
                                                               "_S5b3"),]$tme_subtype_ext
slide6@meta.data$tme_subtype_ext <- NA
slide6@meta.data$tme_subtype_ext <- pData(all_slides)[paste0(row.names(slide6@meta.data), 
                                                               "_S6"),]$tme_subtype_ext

slide_params <- list("slide1" = c("_S1", slide1),
                     "slide2" = c("_S2", slide2),
                     "slide3" = c("_S3", slide3),
                     "slide4" = c("_S4", slide4),
                     "slide5" = c("_S5", slide5),
                     "slide5b2" = c("_S5b2", slide5b2),
                     "slide5b3" = c("_S5b3", slide5b3),
                     "slide6" = c("_S6", slide6))

# Classical example with fibro
potent_subgraphs <- comp_sub[!is.na(comp_sub$center_cell_type) & comp_sub$center_cell_type == "Tumor Classical (1,2,3,5,7,8,11,14)" & comp_sub$Fibroblast > 0,]$subgraph
set.seed(20)
for(subg in potent_subgraphs[sample(1:length(potent_subgraphs), size = 6)]) {
  slide_string <- strsplit(subg, "_", )[[1]][1]
  png(paste0(home_path, "man_plots/class_fib_example_tme_", subg, ".png"), width = 8, height = 5, 
    res = 300, units = "in")
  plot_tme(slide_params[[slide_string]][[2]], 
           gsub(slide_params[[slide_string]][[1]], "", comp_sub_long[comp_sub_long$subgraph == subg,]$cell_members), 
           "tme_subtype_ext", cols = temp_cols, with_conns = TRUE, 
           edge_df = all_edges, slide_suffix = slide_params[[slide_string]][[1]], 
           linewidth = 0.8, buffer = 0) + ggtitle(subg)
  dev.off()
}

# Classical example no fibro
potent_subgraphs <- comp_sub[!is.na(comp_sub$center_cell_type) & comp_sub$center_cell_type == "Tumor Classical (1,2,3,5,7,8,11,14)" & comp_sub$Fibroblast == 0,]$subgraph
set.seed(20)
for(subg in potent_subgraphs[sample(1:length(potent_subgraphs), size = 6)]) {
  slide_string <- strsplit(subg, "_", )[[1]][1]
  png(paste0(home_path, "man_plots/class_nofib_example_tme_", subg, ".png"), width = 8, height = 5, 
    res = 300, units = "in")
  plot_tme(slide_params[[slide_string]][[2]], 
           gsub(slide_params[[slide_string]][[1]], "", comp_sub_long[comp_sub_long$subgraph == subg,]$cell_members), 
           "tme_subtype_ext", cols = temp_cols, with_conns = TRUE, 
           edge_df = all_edges, slide_suffix = slide_params[[slide_string]][[1]], 
           linewidth = 0.8, buffer = 0) + ggtitle(subg)
  dev.off()
}

# Basal example with fibro
potent_subgraphs <- comp_sub[!is.na(comp_sub$center_cell_type) & comp_sub$center_cell_type == "Tumor Basal (4,15)" & comp_sub$Fibroblast > 0,]$subgraph
set.seed(20)
for(subg in potent_subgraphs[sample(1:length(potent_subgraphs), size = 6)]) {
  slide_string <- strsplit(subg, "_", )[[1]][1]
  png(paste0(home_path, "man_plots/basal_fib_example_tme_", subg, ".png"), width = 8, height = 5, 
    res = 300, units = "in")
  plot_tme(slide_params[[slide_string]][[2]], 
           gsub(slide_params[[slide_string]][[1]], "", comp_sub_long[comp_sub_long$subgraph == subg,]$cell_members), 
           "tme_subtype_ext", cols = temp_cols, with_conns = TRUE, 
           edge_df = all_edges, slide_suffix = slide_params[[slide_string]][[1]], 
           linewidth = 0.8, buffer = 0) + ggtitle(subg)
  dev.off()
}

# Basal example no fibro
potent_subgraphs <- comp_sub[!is.na(comp_sub$center_cell_type) & comp_sub$center_cell_type == "Tumor Basal (4,15)" & comp_sub$Fibroblast == 0,]$subgraph
set.seed(20)
for(subg in potent_subgraphs[sample(1:length(potent_subgraphs), size = 6)]) {
  slide_string <- strsplit(subg, "_", )[[1]][1]
  png(paste0(home_path, "man_plots/basal_nofib_example_tme_", subg, ".png"), width = 8, height = 5, 
    res = 300, units = "in")
  plot_tme(slide_params[[slide_string]][[2]], 
           gsub(slide_params[[slide_string]][[1]], "", comp_sub_long[comp_sub_long$subgraph == subg,]$cell_members), 
           "tme_subtype_ext", cols = temp_cols, with_conns = TRUE, 
           edge_df = all_edges, slide_suffix = slide_params[[slide_string]][[1]], 
           linewidth = 0.8, buffer = 0) + ggtitle(subg)
  dev.off()
}

rm(slide_params)
gc()

pData(tumor)$slide_fov <- paste0(pData(tumor)$tissue, "_", pData(tumor)$fov)
fib_list <- as.data.frame(table(pData(tumor)[,c("slide_fov", "comb_labels_clusts")]))
fib_list <- fib_list[fib_list$comb_labels_clusts == "Tumor Mixed / Fibro adjacent (12)",]
fib_list <- fib_list[order(fib_list$Freq, decreasing = T),]

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
temp_cols <- c(temp_cols, all_cols, tumor_cols)

make_fig3_fib_plots <- function(slide_num, fov_num, slide_obj, minval = 0.25) {
  fov_met <- slide_obj@meta.data[slide_obj@meta.data$fov == fov_num,]
  fov.crop <- Crop(slide_obj[[paste0("Slide", slide_num)]], 
                   y = c(min(fov_met$CenterX_global_px),
                         max(fov_met$CenterX_global_px)), 
                   x = c(min(fov_met$CenterY_global_px), 
                         max(fov_met$CenterY_global_px)))
  slide_obj[[paste0("FOV", fov_num)]] <- fov.crop
  DefaultBoundary(slide_obj[[paste0("FOV", fov_num)]]) <- "segmentation"
  
  png(paste0(home_path, "man_plots/S", slide_num, "_FOV", fov_num, "_CEACAM6_KRT6A_COL1A1_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), molecules = c("CEACAM6", "KRT6A", "COL1A1"), nmols = 100000, dark.background = FALSE,
                     alpha = 0.6, border.size = 0.1, border.color = "darkgray", cols = "lightgray", mols.cols = c("#A50104", "#324376", "#0E7C7B")) +
          theme(legend.title = element_blank()) +
          ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
          scale_y_reverse() + scale_x_reverse())
  dev.off()
  
  
  slide_obj@meta.data$tumor_labels <- NA
  slide_obj@meta.data$tumor_labels <- pData(all_slides)[paste0(row.names(slide_obj@meta.data), 
                                                               "_S", slide_num),]$tumor_labels
  
  
  png(paste0(home_path, "man_plots/S", slide_num, "_FOV", fov_num, "_tumor_subtype_light.png"), 
      width = 9, height = 11, res = 300, units = "in")
  print(ImageDimPlot(slide_obj, fov = paste0("FOV", fov_num), #molecules = c("CEACAM6", "KRT17"), nmols = 100000,
                     group.by = "tumor_labels", cols = temp_cols, na.value = "lightgrey",
                     alpha = 0.8, border.size = 0.1, border.color = "darkgray", dark.background = FALSE) +
          theme(legend.title = element_blank()) +
          ggtitle(paste0("Slide ", slide_num, ", FOV", fov_num)) + 
          scale_y_reverse() + scale_x_reverse())
  dev.off()
}

options(future.globals.maxSize = 6.291e+8) # 600MiB

slide1 <- readRDS(paste0(int_path, "slide1.RDS"))
slide2 <- readRDS(paste0(int_path, "slide2.RDS"))
slide6 <- readRDS(paste0(int_path, "slide6.RDS"))

make_fig3_fib_plots(1, 18, slide1) #S1_P1 C2_D6
#make_fig3_fib_plots(2, 12, slide2) #S2_P1 C2_D10
make_fig3_fib_plots(6, 15, slide6) #S6_P2 C3_D6

#comp_mat <- read.csv("/path/to/composition_full_qc_20.csv")
#meta <- read.csv("/path/to/tme_subtype/composition_meta.csv")

#comp_sub <- merge(comp_mat, meta, by = c("Subgraph.ID", "FOV"))

#comp_sub$slide <- stringr::str_split_fixed(comp_sub$fov, "_", 2)[,1]
#comp_sub$slide <- gsub('slide', "S", comp_sub$slide)
#comp_sub$center_cell_id <- paste0(comp_sub$center_cell, "_", comp_sub$slide)
#row.names(comp_sub) <- comp_sub$center_cell_id
#comp_sub$subgraph <- paste0(comp_sub$fov, "_", comp_sub$Subgraph.ID)

#pData(tumor)$comb_labels_clusts <- as.character(pData(tumor)$comb_labels_clusts)

#comp_sub$center_cell_type <- 
#  pData(tumor)[comp_sub$center_cell_id,]$comb_labels_clusts
 
for_plot <- comp_sub[comp_sub$center_cell_type %in%  c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"),] 

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

for_plot_count <- for_plot_count[c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"), c("Ductal", "Classical", "Proliferative", "Mixed", "Basal", "Fib.high", "Acinar", "B.cell.Plasma", "Dendritic", "Endocrine", "Endothelial", "Fibroblast", "Macrophage", "Mast", "Stellate", "T.cell")]

colmins <- apply(for_plot_count, 2, min)
colmax <- apply(for_plot_count, 2, max)
colrange <- colmax - colmins
for_plot_count_sc <- t((t(for_plot_count) - colmins)/colrange)

png(paste0(home_path, "man_plots/F3N_tumor_neighbors_heatmap.png"), width = 9, height = 3.5, 
    res = 300, units = "in")
print(pheatmap::pheatmap(for_plot_count_sc, cluster_rows = F, cluster_cols = F, 
                   color = RColorBrewer::brewer.pal(9, "Blues"), 
                   display_numbers = round(for_plot_count * 100, digits = 1), 
                   number_color = "darkgrey"))
dev.off()

pData(tumor)$patient_label <- factor(pData(tumor)$tissue)
temp <- setNames(names(patient_labels), patient_labels)
pData(tumor)$patient_label <- forcats::fct_recode(pData(tumor)$patient_label, !!!temp)

print(for_plot[for_plot$comb_labels_clusts %in% c("Tumor Mixed (10)", "Tumor Basal (4,15)"),])
for_plot <- as.data.frame(table(pData(tumor)[,c("patient_label", "comb_labels_clusts")]))
#for_plot <- for_plot[for_plot$comb_labels_clusts != "Tumor Other (13)",]
#for_plot$label <- unlist(patient_ROI_labels[for_plot$patient_label])
for_plot$comb_labels_clusts <- factor(for_plot$comb_labels_clusts, levels = c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)", "Tumor Other (13)"))

for_plot <- for_plot %>% group_by(patient_label) %>% mutate(total = sum(Freq))
plot_order <- for_plot[for_plot$comb_labels_clusts %in% c("Tumor Basal (4,15)", "Tumor Mixed (10)"),]
plot_order <- plot_order %>% group_by(patient_label) %>% mutate(perc = sum(Freq)/total)
plot_order <- plot_order[plot_order$comb_labels_clusts == "Tumor Basal (4,15)",]
for_plot$patient_label <- factor(for_plot$patient_label, levels = c(plot_order[order(plot_order$perc),]$patient_label))

png(paste0(home_path, "man_plots/F3S2A_smi_tumor_types_freq_bar.png"), width = 2, height = 2, 
    res = 300, units = "in")
ggplot(data = for_plot, aes(patient_label, Freq * 100, fill = comb_labels_clusts)) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw(base_size = 6) + 
  labs(x = "", y = "Fraction of cells") + 
  monocle3:::monocle_theme_opts() + 
  theme(legend.key.size = unit(5, "point")) + 
  guides(fill = guide_legend(override.aes = list(size=2))) +
  scale_fill_manual(values = tumor_cols, name = "Cell type") + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 60, hjust=1)) #vjust = 0.5, 
dev.off()


for_plot_perc <- for_plot[for_plot$comb_labels_clusts %in% c("Tumor Mixed (10)", "Tumor Basal (4,15)"),]
for_plot_perc <- for_plot_perc %>% group_by(patient_label, total) %>% summarize(mix_bas = sum(Freq))
for_plot_perc$perc <- for_plot_perc$mix_bas/for_plot_perc$total * 100

png(paste0(home_path, "man_plots/F3S2A_smi_basal_mixed_colorbar.png"), width = 2, height = 0.7, 
    res = 300, units = "in")
ggplot(data = for_plot_perc, aes(patient_label, 100, fill = perc, label=round(perc, digits = 1))) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw(base_size = 6) + 
  labs(x = "", y = "Fraction of cells") + 
  monocle3:::monocle_theme_opts() + 
  theme(legend.key.size = unit(5, "point")) + 
  guides(fill = guide_legend(override.aes = list(size=2))) +
  scale_fill_gradient(low = "lightgrey", high = '#96257D') +
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  geom_text(aes(y = 0.5), size = 2)
dev.off()

for_plot <- pData(tumor)[,c("comb_labels_clusts", "ucell_collisson_classical_sig")]
for_plot <- for_plot[!for_plot$comb_labels_clusts %in% c("Tumor Other (13)"),]
for_plot$comb_labels_clusts <- factor(for_plot$comb_labels_clusts, levels = c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"))

png(paste0(home_path, "man_plots/F3S2B_tumor_classical_sig_boxplot.png"), 
    width = 1.5, height = 3, res = 300, units = "in")
ggplot(data = as.data.frame(for_plot), 
       aes(comb_labels_clusts, ucell_collisson_classical_sig, fill = comb_labels_clusts)) + 
  geom_boxplot(outlier.size = 0.1, size = 0.2) + theme_bw(base_size = 6) + 
  NoLegend() +
  scale_fill_manual(values = tumor_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  labs(x = "", y = "Signature value") 
dev.off()

for_plot <- pData(tumor)[,c("comb_labels_clusts", "ucell_moffitt_basal_sig")]
for_plot <- for_plot[!for_plot$comb_labels_clusts %in% c("Tumor Other (13)"),]
for_plot$comb_labels_clusts <- factor(for_plot$comb_labels_clusts, levels = c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"))


png(paste0(home_path, "man_plots/F3S2C_tumor_basal_sig_boxplot.png"), 
    width = 1.5, height = 3, res = 300, units = "in")
ggplot(data = as.data.frame(for_plot), 
       aes(comb_labels_clusts, ucell_moffitt_basal_sig, fill = comb_labels_clusts)) + 
  geom_boxplot(outlier.size = 0.1, size = 0.2) + theme_bw(base_size = 6) + 
  NoLegend() +
  scale_fill_manual(values = tumor_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  labs(x = "", y = "Signature value") 
dev.off()

markers_l <- all_sig[all_sig$Signature_name == "Syng.U.State.Hypoxia_Hallmark",]$Gene_list

tumor_plot <- as.Seurat(tumor, assay = NULL)
tumor_plot <- tumor_plot[,!tumor_plot$comb_labels_clusts %in% c("Tumor Other (13)")]
tumor_plot$comb_labels_clusts <- forcats::fct_rev(tumor_plot$comb_labels_clusts)
png(paste0(home_path, "man_plots/tumor_dotplot_hypox_markers.png"), width = 13, height = 4, 
    res = 300, units = "in")
DotPlot(object = tumor_plot, 
        features = markers_l[markers_l %in% row.names(fData(all_slides))], 
        group.by = "comb_labels_clusts") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))  +
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 10),
        legend.key.size = unit(0.7, "lines"))#+
#scale_color_gradient(low = "lightgrey", high = "blue", oob = scales::squish, limits = c(0, 1))
dev.off()

tumor_plot <- as.Seurat(tumor, assay = NULL)
tumor_plot <- UCell::AddModuleScore_UCell(tumor_plot, features=list(markers_l[markers_l %in% row.names(fData(all_slides))]), name="Syng.U.State.Hypoxia_Hallmark")
pData(tumor)$ucell_hypox_prolif_sig <- tumor_plot@meta.data[row.names(tumor_plot@meta.data),]$signature_1Syng.U.State.Hypoxia_Hallmark

for_plot <- pData(tumor)[,c("comb_labels_clusts", "ucell_hypox_prolif_sig")]
for_plot <- for_plot[!for_plot$comb_labels_clusts %in% c("Tumor Other (13)"),]
for_plot$comb_labels_clusts <- factor(for_plot$comb_labels_clusts, levels = c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"))

png(paste0(home_path, "man_plots/F3S2D_tumor_hypox_sig_boxplot.png"), 
    width = 1.5, height = 3, res = 300, units = "in")
ggplot(data = as.data.frame(for_plot), 
       aes(comb_labels_clusts, ucell_hypox_prolif_sig, fill = comb_labels_clusts)) + 
  geom_boxplot(outlier.size = 0.1, size = 0.2) + theme_bw(base_size = 6) + 
  NoLegend() +
  scale_fill_manual(values = tumor_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  labs(x = "", y = "Signature value") 
dev.off()

temp <- tumor[,sample(1:nrow(pData(tumor)), nrow(pData(tumor)), 
                      replace = FALSE)]

minval <- min(c(pData(temp)$ucell_moffitt_basal_sig, 
                pData(temp)$ucell_collisson_classical_sig,
                pData(temp)$ucell_peng_ductal1_sig))
maxval <- max(c(pData(temp)$ucell_moffitt_basal_sig, 
                pData(temp)$ucell_collisson_classical_sig,
                pData(temp)$ucell_peng_ductal1_sig))
s2 <- plot_cells(temp, color_cells_by = "ucell_moffitt_basal_sig", 
                 cell_size = 0.1, alpha = 0.4) + 
  ggtitle("Moffitt Basal Signature") + theme_classic(base_size = 6) + 
  NoLegend() + scale_color_distiller(palette = "OrRd", direction = 1, limits = c(minval, maxval))

s3 <- plot_cells(temp, color_cells_by = "ucell_collisson_classical_sig", 
                 cell_size = 0.1, alpha = 0.4) + 
  ggtitle("Collisson Classical Signature") + theme_classic(base_size = 6) + 
  NoLegend() + scale_color_distiller(palette = "OrRd", direction = 1, limits = c(minval, maxval))

s6 <- plot_cells(temp, color_cells_by = "ucell_peng_ductal1_sig", 
                 cell_size = 0.1, alpha = 0.4) + 
  ggtitle("Peng2019 Ductal1 Signature") + theme_classic(base_size = 6) + 
  scale_color_distiller(palette = "OrRd", direction = 1, limits = c(minval, maxval))

png(paste0(home_path, "man_plots/F3S2E_tumor_only_signatures_ucell_legend.png"), width = 6, 
    height = 2, res = 300, units = "in")
s6 + s3 + s2
dev.off()

png(paste0(home_path, "man_plots/F3S2E_tumor_only_signatures_ucell.png"), width = 6, 
    height = 2, res = 300, units = "in")
(s6 + NoLegend()) + s3 + s2
dev.off()



for_plot <- comp_sub[comp_sub$center_cell_type %in%  c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"),]
for_plot_count <- for_plot %>% group_by(center_cell_type) %>% 
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

for_plot_count <- for_plot_count[c("Ductal (6)", "Tumor Classical (1,2,3,5,7,8,11,14)", "Tumor Classical / Proliferative (9)", "Tumor Mixed (10)", "Tumor Basal (4,15)", "Tumor Mixed / Fibro adjacent (12)"), c("Ductal", "Classical", "Proliferative", "Mixed", "Basal", "Fib.high", "Acinar", "B.cell.Plasma", "Dendritic", "Endocrine", "Endothelial", "Fibroblast", "Macrophage", "Mast", "Stellate", "T.cell")]

png(paste0(home_path, "man_plots/F3S2G_tumor_neighbors_average_heatmap.png"), width = 9, height = 3.5, 
    res = 300, units = "in")
print(pheatmap::pheatmap(for_plot_count, cluster_rows = F, cluster_cols = F, 
                   color = RColorBrewer::brewer.pal(9, "Blues"), 
                   display_numbers = round(for_plot_count * 100, digits = 1), 
                   number_color = "darkgrey"))
dev.off()

sessionInfo()
