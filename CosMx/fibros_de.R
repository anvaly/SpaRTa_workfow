library(tidyverse)
library(monocle3)
library(Seurat)
library(SeuratData)


home_path <- "/out/cosmx_TAP/"
int_path <- "/out/cosmx_TAP/intermediate/"

tumor <- readRDS(file = paste0(int_path, "all_slides_tumor_integrated.RDS"))

fibrosde1=data.frame()
fibrosde2=data.frame()
for (i in c(1:9)) {
    comp_mat <- read.csv(paste0(int_path, "TME_analysis/tme_subtype/composition_full_qc_20_", as.character(i), "hop" ,".csv"))
    meta <- read.csv(paste0(int_path, "TME_analysis/tme_subtype/composition_meta_", as.character(i), "hop" , ".csv"))
    
    
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
    comp_sub_long$center_cell_type <- 
        pData(tumor)[comp_sub_long$center_cell_id,]$comb_labels_clusts
    
    comp_sub_long$center_cell_type[is.na(comp_sub_long$center_cell_type)]='not tumor'
    nearductal=comp_sub_long$cell_members[comp_sub_long$center_cell_type=="Ductal (6)"]
    nearbasal=comp_sub_long$cell_members[comp_sub_long$center_cell_type=="Tumor Basal (4,15)"]
    
    
    ##get all fibros
    fibro <- readRDS(file=paste0(int_path, "all_slides_fibro_integrated.RDS"))
    fibro_seur <- as.Seurat(fibro, data = NULL)
    
    fibids=colnames(fibro_seur)
    ductfibs=as.integer(fibids %in% nearductal)
    basfibs=as.integer(fibids %in% nearbasal)
    fibrobin=ductfibs*2+basfibs ## 0=neither, 1=basal, 2=ductal, 3=both
    
    table(fibrobin)
    
    fibro_seur@meta.data$fibrobin=fibrobin
    fibro_seur=Seurat::SetIdent(fibro_seur, value="fibrobin")
    
    bas_v_duct=FindMarkers(fibro_seur, ident.1=1, ident.2=2,
                           min.pct=0.1,
                           logfc.threshold=0.25,
                           test.use = "wilcox")
    bas_v_duct$name=rownames(bas_v_duct)
    bas_v_duct$comp='basal v ductal'
    
    bas_v_neither=FindMarkers(fibro_seur, ident.1=1, ident.2=0,
                              min.pct=0.1,
                              logfc.threshold=0.25,
                              test.use = "wilcox")
    bas_v_neither$name=rownames(bas_v_neither)
    bas_v_neither$comp='basal v neither'
    
    duct_v_neither=FindMarkers(fibro_seur, ident.1=2, ident.2=0,
                               min.pct=0.1,
                               logfc.threshold=0.25,
                               test.use = "wilcox")
    duct_v_neither$name=rownames(duct_v_neither)
    duct_v_neither$comp='ductal v neither'
    
    all1=rbind(bas_v_duct, bas_v_neither, duct_v_neither)
    all1$logpval=log10(all1$p_val_adj)
    minlogp=min(all1$logpval[all1$logpval!=-Inf])-50
    all1$logpval[all1$logpval==-Inf]=minlogp
    all1$logpval=-all1$logpval
    all1$nhops=i
    
    fibrosde1=rbind(fibrosde1, all1)

    
    neartumor=comp_sub_long$cell_members[!comp_sub_long$center_cell_type=='not tumor']
    farfibs=!fibids %in% neartumor
    
    fibrobin[farfibs]=-1 ## 0=neither, 1=basal, 2=ductal, 3=both
    fibro_seur@meta.data$fibrobin=fibrobin
    fibro_seur=Seurat::SetIdent(fibro_seur, value="fibrobin")
    
    bas_v_far=FindMarkers(fibro_seur, ident.1=1, ident.2=-1,
                          min.pct=0.1,
                          logfc.threshold=0.25,
                          test.use = "wilcox")
    bas_v_far$name=rownames(bas_v_far)
    bas_v_far$comp='basal v no nearby tumor'
    
    duct_v_far=FindMarkers(fibro_seur, ident.1=2, ident.2=-1,
                           min.pct=0.1,
                           logfc.threshold=0.25,
                           test.use = "wilcox")
    duct_v_far$name=rownames(duct_v_far)
    duct_v_far$comp='ductal v no nearby tumor'
    
    all2=rbind(bas_v_duct, bas_v_far, duct_v_far)
    all2$logpval=log10(all2$p_val_adj)
    minlogp=min(all2$logpval[all2$logpval!=-Inf])-50
    all2$logpval[all2$logpval==-Inf]=minlogp
    all2$logpval=-all2$logpval
    all2$nhops=i

    fibrosde2=rbind(fibrosde2, all2)
}



library(RColorBrewer)
library(ggrepel)
colors=c('#000000','#FC8D62')

##fibrosde1=read.csv(paste0(home_path, 'man_plots/diff_exp1_full.csv'))
##fibrosde2=read.csv(paste0(home_path, 'man_plots/diff_exp2_full.csv'))

pdf(paste0(home_path, 'man_plots/diff_exp1.pdf'), w=9, h=3.5)
for (i in c(1:9)) {
    forplot=fibrosde1[fibrosde1$nhops==i,]
    forplot$sig=abs(forplot$avg_log2FC)>=.75
    print(ggplot(forplot, aes(x=avg_log2FC, y=logpval, colour=sig)) +
          geom_point() +
          geom_label_repel(data=forplot[forplot$sig,], aes(avg_log2FC, logpval, label=name, colour=sig), size=2) +
          ##geom_text_repel(data=forplot[forplot$sig,], aes(avg_log2FC, logpval, label=name, colour=sig), size=3, angle=90, hjust=1, nudge_y=-7) +
          scale_colour_manual(values=colors) +
          facet_wrap(~comp) +
          ggtitle(paste0('Diff Exp: Fibros within/without ', as.character(i), ' hop neighborhoods of tumor types')) +
          theme_bw()+
          theme(legend.position="none"))
}
dev.off()


pdf(paste0(home_path, 'man_plots/diff_exp2.pdf'), w=9, h=3.5)
for (i in c(1:9)) {
    forplot=fibrosde2[fibrosde2$nhops==i,]
    forplot$sig=abs(forplot$avg_log2FC)>=.75
    print(ggplot(forplot, aes(x=avg_log2FC, y=logpval, colour=sig)) +
          geom_point() +
          geom_label_repel(data=forplot[forplot$sig,], aes(avg_log2FC, logpval, label=name, colour=sig), size=2) +
          ##geom_text(data=forplot[forplot$sig,], aes(avg_log2FC, logpval, label=name, colour=sig), check_overlap = TRUE, size=3, angle=90, hjust=1, nudge_y=-7) +
          scale_colour_manual(values=colors) +
          facet_wrap(~comp) +
          ggtitle(paste0('Diff Exp: Fibros within/without ', as.character(i), ' hop neighborhoods of tumor types')) +
          theme_bw()+
          theme(legend.position="none"))
}
dev.off()


rownames(fibrosde1) <- NULL
rownames(fibrosde2) <- NULL

write.csv(fibrosde1, paste0(home_path, 'man_plots/diff_exp1_full.csv'))
write.csv(fibrosde2, paste0(home_path, 'man_plots/diff_exp2_full.csv'))

write.csv(fibrosde1[fibrosde1$nhops==3,], paste0(home_path, 'man_plots/diff_exp1_3hop.csv'))
write.csv(fibrosde2[fibrosde1$nhops==3,], paste0(home_path, 'man_plots/diff_exp2_3hop.csv'))

write.csv(fibrosde1[abs(fibrosde1$avg_log2FC)>=.75,], paste0(home_path, 'man_plots/diff_exp1.csv'))
write.csv(fibrosde2[abs(fibrosde2$avg_log2FC)>=.75,], paste0(home_path, 'man_plots/diff_exp2.csv'))
