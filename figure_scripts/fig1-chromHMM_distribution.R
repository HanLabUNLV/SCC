setwd("~/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig0-qc_chromhmm/")

library(readr)
library(dplyr)
library(ggplot2)

endo_bed <- read_tsv("../../create_cluster_matrix_100perc_chromHMM_in_window/chromHMM_epimap/endoderm_200bp_window.bed",
                     col_names = c("chr","pos1","pos2","label","score","strand","x","y","z"))
h1_bed <- read_tsv("../../create_cluster_matrix_100perc_chromHMM_in_window/chromHMM_epimap/h1_200bp_window.bed",
                   col_names = c("chr","pos1","pos2","label","score","strand","x","y","z"))
hff_bed <- read_tsv("../../create_cluster_matrix_100perc_chromHMM_in_window/chromHMM_epimap/HFF_200bp_window.bed",
                    col_names = c("chr","pos1","pos2","label","score","strand","x","y","z"))
hela_bed <- read_tsv("../../create_cluster_matrix_100perc_chromHMM_in_window/chromHMM_epimap/HeLa-S3_200bp_window.bed",
                     col_names = c("chr","pos1","pos2","label","score","strand","x","y","z"))

endo <- cbind.data.frame(data.frame(table(endo_bed$label)),cell = 'Endoderm')
hff <- cbind.data.frame(data.frame(table(hff_bed$label)),cell = 'HFF')
h1 <- cbind.data.frame(data.frame(table(h1_bed$label)),cell = 'H1')
hela <- cbind.data.frame(data.frame(table(hela_bed$label)),cell = 'HeLa')

counts <- rbind.data.frame(endo,hff,h1,hela)


ss_counts <- counts[grepl("Quies",counts$Var1),]
png("ChromHMM_dist_justQuies.png",height = 600,width = 1200,res = 200)
png("ChromHMM_dist_justQuies.png",height = 550,width = 460,res = 200)

ggplot(ss_counts,aes(x = Var1,y = Freq,fill = cell)) +
  geom_bar(stat = "identity",position = position_dodge()) + 
  labs(x = "ChromHMM label", y = "Num 200bp windows",fill = "Cell Type",title = "ChromHMM 200bp window labels") +
  theme(axis.text.x = element_text(angle = 60, hjust=1))
dev.off()

sum(endo$Freq)
sum(hff$Freq)
sum(hela$Freq)
sum(h1$Freq)


# #look at base regions:
# 
# endo_all <- read_tsv("../fig4-all_sankey/create_bed_network_files/endo_labelled_200bp_window.bed",col_names = F)
# hff_all <- read_tsv("../fig4-all_sankey/create_bed_network_files/hff_labelled_200bp_window.bed",col_names = F)
# h1_all <- read_tsv("../fig4-all_sankey/create_bed_network_files/h1_labelled_200bp_window.bed",col_names = F)
# hela_all <- read_tsv("../fig4-all_sankey/create_bed_network_files/hela_labelled_200bp_window.bed",col_names = F)


