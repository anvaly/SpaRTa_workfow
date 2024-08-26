library(tidyverse)
library(monocle3)

home_path <- "/out/cosmx_TAP/"
int_path <- "/out/cosmx_TAP/intermediate/"

tumor <- readRDS(file = paste0(int_path, "all_slides_tumor_integrated.RDS"))

plotfile=paste0(home_path, 'man_plots/nb_population.pdf')
pdf(plotfile, w=11, h=8)
for (i in c(1,3)) {
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
    
    counts=as.data.frame(table(comp_sub_long$center_cell_id))
    colnames(counts)=c('center_cell', 'nb_pop')
    counts$center_type <- 
        pData(tumor)[counts$center_cell,]$comb_labels_clusts

    ##overall nb population histo
    maxpop=max(counts$nb_pop)
    forplot=counts[!is.na(counts$center_type),]
    forplot_sum=data.frame(popmean=mean(counts$nb_pop), popmedian=median(counts$nb_pop))
    
    plot=ggplot(counts, aes(x=nb_pop, alpha=.3)) +
        geom_histogram(binwidth=1) +
        geom_vline(aes(xintercept=mean(nb_pop))) +
        geom_text(data=forplot_sum,
                  aes(x=-Inf,
                      y=Inf,
                      hjust=0,
                      vjust=1,
                      label=paste0('mean: ', as.character(popmean), '\nmedian: ', as.character(popmedian))),
                  inherit.aes=FALSE) +
        xlim(0, maxpop+2) +
        ggtitle(paste0(as.character(i), ' hop neighborhood populations')) +
        xlab('Neighborhood population') +
        theme_bw()
    print(plot)

    types=names(table(counts$center_type))

    forfacet=counts[!is.na(counts$center_type),]
    forfacet_sum=as_tibble(forfacet) %>%
        group_by(center_type) %>%
        summarise(popmean=mean(nb_pop), popmedian=median(nb_pop)) %>%
        mutate(popmean=round(popmean, 2))

    plot=ggplot(forfacet, aes(x=nb_pop, alpha=.3)) +
        geom_histogram(binwidth=1) +
        geom_text(data=forfacet_sum,
                  aes(x=-Inf,
                      y=Inf,
                      hjust=0,
                      vjust=1,
                      label=paste0('mean: ', popmean, '\nmedian: ', popmedian)),
                  size=3,
                  inherit.aes=FALSE) +
        facet_wrap(~center_type, scales='free_y') +
        ggtitle(paste0(as.character(i), ' hop neighborhood populations')) +
        xlab('Neighborhood population') +
        xlim(0, maxpop+2) +
        theme_bw()
    print(plot)

}
dev.off()
