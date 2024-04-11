
setwd("~/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig4-all_sankey//")
library(tidyverse)
library(ggsignif)
library(circlize)
library(ComplexHeatmap)

final_net <- readRDS("create_bed_network_files/final_net.rds")

tot_endo <- sum(!is.na(final_net$endo_base_region))
tot_h1 <- sum(!is.na(final_net$h1_base_region))
tot_hff <- sum(!is.na(final_net$hff_base_region))
tot_hela <- sum(!is.na(final_net$hela_base_region))


# ENDO

# Endo - H1 40.7%
sum(final_net$endo_base_region == final_net$h1_base_region & 
      final_net$endo_cluster != final_net$h1_cluster,na.rm = T)/tot_endo

# H1 - Endo 19.9%
sum(final_net$endo_base_region == final_net$h1_base_region & 
      final_net$endo_cluster != final_net$h1_cluster,na.rm = T)/tot_h1

# Endo - hff 26.2%
sum(final_net$endo_base_region == final_net$hff_base_region & 
      final_net$endo_cluster != final_net$hff_cluster,na.rm = T)/tot_endo

# hff - Endo 15.8%
sum(final_net$endo_base_region == final_net$hff_base_region & 
      final_net$endo_cluster != final_net$hff_cluster,na.rm = T)/tot_hff

# Endo - hela 26.9%
sum(final_net$endo_base_region == final_net$hela_base_region & 
      final_net$endo_cluster != final_net$hela_cluster,na.rm = T)/tot_endo

# hela - Endo 12.5%
sum(final_net$endo_base_region == final_net$hela_base_region & 
      final_net$endo_cluster != final_net$hela_cluster,na.rm = T)/tot_hela

## H1

# h1 - hff 18.3%
sum(final_net$h1_base_region == final_net$hff_base_region & 
      final_net$h1_cluster != final_net$hff_cluster,na.rm = T)/tot_h1

# hff - h1 22.5%
sum(final_net$h1_base_region == final_net$hff_base_region & 
      final_net$h1_cluster != final_net$hff_cluster,na.rm = T)/tot_hff

# h1 - hela 16.4%
sum(final_net$h1_base_region == final_net$hela_base_region & 
      final_net$h1_cluster != final_net$hela_cluster,na.rm = T)/tot_h1

# hela - h1 15.6%
sum(final_net$h1_base_region == final_net$hela_base_region & 
      final_net$h1_cluster != final_net$hela_cluster,na.rm = T)/tot_hela


## HFF

# hff - hela 23.4%
sum(final_net$hff_base_region == final_net$hela_base_region & 
      final_net$hff_cluster != final_net$hela_cluster,na.rm = T)/tot_hff

# hela - hff 18.1%
sum(final_net$hff_base_region == final_net$hela_base_region & 
      final_net$hff_cluster != final_net$hela_cluster,na.rm = T)/tot_hela
