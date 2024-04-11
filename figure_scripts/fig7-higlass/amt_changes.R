

setwd("~/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig6-higlass")

library(circlize)
library(ComplexHeatmap)
library(readr)
library(dplyr)

h1_net <- read_tsv("../fig4-all_sankey/create_bed_network_files/network_200bp_window_h1.bed", col_types = 'cccccccccccccc',
                   col_names = c("h1_chr","h1_pos1","h1_pos2","h1_base_region","h1_cluster","h1_cell_type","match",
                                 "2_chr","2_pos1","2_pos2","2_base_region","2_cluster","2_cell_type","overlap"))

endo_net <- read_tsv("../fig4-all_sankey/create_bed_network_files/network_200bp_window_endo.bed", col_types = 'cccccccccccccc',
                   col_names =c("endo_chr","endo_pos1","endo_pos2","endo_base_region","endo_cluster","endo_cell_type","match",
                                  "2_chr","2_pos1","2_pos2","2_base_region","2_cluster","2_cell_type","overlap"))

endo_all <- read_tsv("../fig4-all_sankey/create_bed_network_files/endo_labelled.bed",col_names = F)

hff_net <- read_tsv("../fig4-all_sankey/create_bed_network_files/network_200bp_window_hff.bed", col_types = 'cccccccccccccc',
                   col_names =c("hff_chr","hff_pos1","hff_pos2","hff_base_region","hff_cluster","hff_cell_type","match",
                                 "2_chr","2_pos1","2_pos2","2_base_region","2_cluster","2_cell_type","overlap"))

hela_net <- read_tsv("../fig4-all_sankey/create_bed_network_files/network_200bp_window_hela.bed", col_types = 'cccccccccccccc',
                   col_names =c("hela_chr","hela_pos1","hela_pos2","hela_base_region","hela_cluster","hela_cell_type","match",
                                  "2_chr","2_pos1","2_pos2","2_base_region","2_cluster","2_cell_type","overlap"))

# endo_net <- endo_net %>% subset(`2_cell_type` == 'hff')


# PICK INTERESTING LABELS
# LOOK FOR INTERESTING GENES

saveRDS(endo_net,"endo_net.rds")

ss <- endo_net %>% subset(`2_base_region` == endo_base_region & endo_cluster != `2_cluster` &
                      grepl("(Tss)|(Enh)",endo_cluster) & !(grepl('Wk',endo_cluster)) &
                        grepl("(Biv)|(Het)",`2_cluster`))

ss <- endo_net %>% subset(`2_base_region` == endo_base_region & endo_base_region == 'TssA' & endo_cluster != `2_cluster` &
                            grepl("(Tss_)|(Enh_)",endo_cluster) & !(grepl('Wk',endo_cluster)) &
                            grepl("Low",`2_cluster`) & `2_cell_type` == 'h1') 

write_tsv(ss, "interesting_regions.bed",col_names = F)
#bedtools intersect -a ../fig5-gene_sankey/gene_exp/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500b.bed -b interesting_regions.bed -wo > interesting_regions_genes.bed


ss_wgenes <- read_tsv("interesting_regions_genes.bed",col_names = F)

gene_exp_endo_hff <- read_tsv("../fig5-gene_sankey/gene_exp/dec_hff_de.tsv")
gene_exp_endo_h1 <- read_tsv("../fig5-gene_sankey/gene_exp/dec_h1_de.tsv")

gene_exp <- gene_exp_endo_hff[gene_exp_endo_hff$gene %in% unique(ss_wgenes$X4),]





# # amount of changes
# matched <- (endo_net[endo_net$`endo_cluster` == endo_net$`2_cluster`,])
# mismatched <- (endo_net[endo_net$`endo_cluster` != endo_net$`2_cluster`,])
# 
# 
# 
# #h1
# mismatched_h1 <- mismatched[mismatched$`2_cell_type` == "h1",]
# 
# # 40.7% have same chromHMM, different CCS cluster
# sum(mismatched_h1$`endo_base_region` == mismatched_h1$`2_base_region`) /
#   (nrow(mismatched_h1))
# 
# #endo
# mismatched_endo <- mismatched[mismatched$`2_cell_type` == "endo",]
# 
# # 40.7% have same chromHMM, different CCS cluster
# sum(mismatched_endo$`endo_base_region` == mismatched_endo$`2_base_region`) /
#   (nrow(mismatched_endo))
# 
# 
# #hff
# mismatched_hff <- mismatched[mismatched$`2_cell_type` == "hff",]
# 
# # 33.3% have same chromHMM, different CCS cluster
# sum(mismatched_hff$`endo_base_region` == mismatched_hff$`2_base_region`) /
#   (nrow(mismatched_hff))
# 
# #hela
# mismatched_hela <- mismatched[mismatched$`2_cell_type` == "hela",]
# 
# # 28.7% have same chromHMM, different CCS cluster
# sum(mismatched_hela$`endo_base_region` == mismatched_hela$`2_base_region`) /
#   (nrow(mismatched_hela))
# 
# 
# # UPDATED 11/10/22
# #Hela
# # endo 1, h1 2, hff3
# # % regions changed clusters =
# #  28.7% of Hela regions in endo
# #  25.1% of Hela regions in h1
# #  32.6% of Hela regions in hff
# 
# #H1
# # endo 1 hela 2 hff 3
# #  40.7% h1 regions in endo
# #  25.1% h1 regions in hela
# #  28.0% h1 regions in hff
# 
# #HFF
# # endo 1, h1 2, hela 3
# #  33.3% of hff regions in endo
# #  28.0 % of hff regions in h1
# #  32.6% of hff regions in hela
# 
# # Endo
# # h1 1, hela 2, hff 3
# #  40.7% of endo regions in h1
# #  28.7% of endo regions in hela
# #  33.3% of endo regions in hff
# 
# 
# table(matched$match) / total_regions
# table(mismatched$match) / total_regions
# 
