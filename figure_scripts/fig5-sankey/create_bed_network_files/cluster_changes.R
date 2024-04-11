

setwd("~/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig5-all_sankey/create_bed_network_files/")
library(tidyverse)
library(ggsignif)
library(circlize)
library(ComplexHeatmap)

# # bash - # cut -f 1,2,3 *200bp_window.bed | sort | uniq > all_regions.txt
# all_regions <- read_tsv("all_regions_uniq.txt",col_names = c("chr","pos1","pos2"))
#   
# 
# h1_net <- read_tsv("h1_labelled_200bp_window.bed", col_types = 'cccccc',
#                    col_names = c("chr","pos1","pos2","h1_base_region","h1_cluster","h1_cell_type"))
# endo_net <- read_tsv("endo_labelled_200bp_window.bed", col_types = 'cccccc',
#                      col_names = c("chr","pos1","pos2","endo_base_region","endo_cluster","endo_cell_type"))
# hff_net <- read_tsv("hff_labelled_200bp_window.bed", col_types = 'cccccc',
#                     col_names = c("chr","pos1","pos2","hff_base_region","hff_cluster","hff_cell_type"))
# hela_net <- read_tsv("hela_labelled_200bp_window.bed", col_types = 'cccccc',
#                      col_names = c("chr","pos1","pos2","hela_base_region","hela_cluster","hela_cell_type"))
# 
# 
# 
# all_h1_net <- merge(all_regions,h1_net,
#                     by = c("chr","pos1","pos2"),all.x = T)
# 
# all_h1_endo_net <- merge(all_h1_net,endo_net,
#                     by = c("chr","pos1","pos2"),all.x = T)
# 
# all_h1_endo_hff_net <- merge(all_h1_endo_net,hff_net,
#                     by = c("chr","pos1","pos2"),all.x = T)
# 
# all_all_net <- merge(all_h1_endo_hff_net,hela_net,
#                     by = c("chr","pos1","pos2"),all.x = T)
# 
# 
# final_net <- all_all_net %>% select(chr,pos1,pos2,
#                                     endo_base_region,endo_cluster,h1_base_region,h1_cluster,
#                                     hela_base_region,hela_cluster,hff_base_region,hff_cluster)
# 
# saveRDS(final_net,"final_net.rds")
final_net <- readRDS("final_net.rds")

final_net_cluster <- final_net %>% select(h1_cluster, endo_cluster, hff_cluster,hela_cluster)


# subcategories:
# TX Contact: cont_Tx_EnhG1, cont_Tx_EnhG2,
# cont_TxWk_Enh, cont_TxWk_EnhWk
# TSS Contact: 	cont_Tss_Enh, cont_Tss_EnhG,
# cont_Tss_noEnh, cont_Tss_EnhWk
# Enhancer Contact: cont_EnhA, cont_EnhA_EnhWk
# Bivalent Contact: cont_Biv
# ReprPC Contact: cont_ReprPC_Biv, cont_ReprPC
# Het Contact: cont_Het_Rpts, cont_Het_Rpts_strong, cont_Het_Rpts_strongest
# Low Contact: cont_Low_TssBiv, cont_Low_EnhWk
# c(
#   7 = NULL
#   6 low contact = "cont_Low_TssBiv", "cont_Low_EnhWk"
#   5 het contact = "cont_Het_Rpts", "cont_Het_Rpts_strong", "cont_Het_Rpts_strongest"
#   4 ReprPC contact =  "cont_ReprPC_Biv","cont_ReprPC"
#   3 biv contact = "cont_Biv"
#   2 enh contact = "cont_EnhA","cont_EnhA_EnhWk"
#   1 tss contact = "cont_Tss_Enh","cont_Tss_EnhG","cont_Tss_noEnh","cont_Tss_EnhWk"
#   0 tx contact = "cont_Tx_EnhG1","cont_Tx_EnhG2","cont_TxWk_Enh","cont_TxWk_EnhWk"
# )



final_net_cluster$endo_cluster <- ifelse(
  is.na(final_net_cluster$endo_cluster),7,ifelse(
  final_net_cluster$endo_cluster %in% c("cont_Low_TssBiv", "cont_Low_EnhWk"),6,ifelse(
  final_net_cluster$endo_cluster %in% c("cont_Het_Rpts",
                                        "cont_Het_Rpts_strong",
                                        "cont_Het_Rpts_strongest"),5,ifelse(
    final_net_cluster$endo_cluster %in% c("cont_ReprPC_Biv","cont_ReprPC"),4, ifelse(
      final_net_cluster$endo_cluster %in% c("cont_Biv"),3, ifelse(
        final_net_cluster$endo_cluster %in% c("cont_EnhA","cont_EnhA_EnhWk"),2, ifelse(
          final_net_cluster$endo_cluster %in% c("cont_Tss_Enh",
                                                "cont_Tss_EnhG",
                                                "cont_Tss_noEnh",
                                                "cont_Tss_EnhWk"),1, 0
        )
      )
    )
  )
)))

final_net_cluster$hff_cluster <- ifelse(
  is.na(final_net_cluster$hff_cluster),7,ifelse(
    final_net_cluster$hff_cluster %in% c("cont_Low_TssBiv", "cont_Low_EnhWk"),6,ifelse(
      final_net_cluster$hff_cluster %in% c("cont_Het_Rpts",
                                            "cont_Het_Rpts_strong",
                                            "cont_Het_Rpts_strongest"),5,
      ifelse(
        final_net_cluster$hff_cluster %in% c("cont_ReprPC_Biv","cont_ReprPC"),4, ifelse(
          final_net_cluster$hff_cluster %in% c("cont_Biv"),3, ifelse(
            final_net_cluster$hff_cluster %in% c("cont_EnhA","cont_EnhA_EnhWk"),2, ifelse(
              final_net_cluster$hff_cluster %in% c("cont_Tss_Enh",
                                                    "cont_Tss_EnhG",
                                                    "cont_Tss_noEnh",
                                                    "cont_Tss_EnhWk"),1, 0
            )
          )
        )
      )
    )))

final_net_cluster$hela_cluster <- ifelse(
  is.na(final_net_cluster$hela_cluster),7,ifelse(
    final_net_cluster$hela_cluster %in% c("cont_Low_TssBiv", "cont_Low_EnhWk"),6,ifelse(
      final_net_cluster$hela_cluster %in% c("cont_Het_Rpts",
                                           "cont_Het_Rpts_strong",
                                           "cont_Het_Rpts_strongest"),5,
      ifelse(
        final_net_cluster$hela_cluster %in% c("cont_ReprPC_Biv","cont_ReprPC"),4, ifelse(
          final_net_cluster$hela_cluster %in% c("cont_Biv"),3, ifelse(
            final_net_cluster$hela_cluster %in% c("cont_EnhA","cont_EnhA_EnhWk"),2, ifelse(
              final_net_cluster$hela_cluster %in% c("cont_Tss_Enh",
                                                   "cont_Tss_EnhG",
                                                   "cont_Tss_noEnh",
                                                   "cont_Tss_EnhWk"),1, 0
            )
          )
        )
      )
    )))

final_net_cluster$h1_cluster <- ifelse(
  is.na(final_net_cluster$h1_cluster),7,ifelse(
    final_net_cluster$h1_cluster %in% c("cont_Low_TssBiv", "cont_Low_EnhWk"),6,ifelse(
      final_net_cluster$h1_cluster %in% c("cont_Het_Rpts",
                                            "cont_Het_Rpts_strong",
                                            "cont_Het_Rpts_strongest"),5,
      ifelse(
        final_net_cluster$h1_cluster %in% c("cont_ReprPC_Biv","cont_ReprPC"),4, ifelse(
          final_net_cluster$h1_cluster %in% c("cont_Biv"),3, ifelse(
            final_net_cluster$h1_cluster %in% c("cont_EnhA","cont_EnhA_EnhWk"),2, ifelse(
              final_net_cluster$h1_cluster %in% c("cont_Tss_Enh",
                                                    "cont_Tss_EnhG",
                                                    "cont_Tss_noEnh",
                                                    "cont_Tss_EnhWk"),1, 0
            )
          )
        )
      )
    )))

# final_net_cluster$h1_cluster <- factor(final_net_cluster$h1_cluster,levels = c(0,1,2,3,4,5,6))

final_net_cluster$endo_cluster <- final_net_cluster$endo_cluster + 8
# final_net_cluster$endo_cluster <- factor(final_net_cluster$endo_cluster,levels = c(7,8,9,10,11,12,13))


final_net_cluster$hff_cluster <- final_net_cluster$hff_cluster + 16
# final_net_cluster$hff_cluster <- factor(final_net_cluster$hff_cluster,levels = c(14,15,16,17,18,19,20))


final_net_cluster$hela_cluster <- final_net_cluster$hela_cluster + 24
# final_net_cluster$hela_cluster <- factor(final_net_cluster$hela_cluster,levels = c(21,22,23,24,25,26,27))
colnames(final_net_cluster) <- c('x','x','x','x')

network <- rbind.data.frame(table(final_net_cluster[,1:2]),table(final_net_cluster[,2:3]),table(final_net_cluster[,3:4]))
# merged_table <- rbind.data.frame(final_net_cluster[,1:2],final_net_cluster[,2:3])

# merged_table$x <- factor(merged_table$x,levels = 0:13)
# merged_table$x.1 <- factor(merged_table$x.1,levels = 7:27)


# network <- as_tibble(table(merged_table))

colnames(network) <- c("source", "target", "value")

# colnames(network) <- c("target", "source", "value")
network$source <- as.numeric(network$source)
network$target <- as.numeric(network$target)
network$value <- as.numeric(network$value)

# # MUST BE ZERO-INDEXED
# network$target <- network$target + 7

network$group <- 'col'

#sm_net <- network[network$n > 200 & network$n < 500,]
# sm_net$target <- sm_net$target + 15


nodes = data.frame("name" = 
                     c("TX Contact",
                       "TSS Contact",
                       "Enhancer Contact",
                       "Bivalent Contact",
                       "ReprPC Contact",
                       "Het Contact",
                       "Low Contact",
                       "No Label",
                       
                       "TX Contact",
                       "TSS Contact",
                       "Enhancer Contact",
                       "Bivalent Contact",
                       "ReprPC Contact",
                       "Het Contact",
                       "Low Contact",
                       "No Label",
                       
                       "TX Contact",
                       "TSS Contact",
                       "Enhancer Contact",
                       "Bivalent Contact",
                       "ReprPC Contact",
                       "Het Contact",
                       "Low Contact",
                       "No Label",
                       
                       "TX Contact",
                       "TSS Contact",
                       "Enhancer Contact",
                       "Bivalent Contact",
                       "ReprPC Contact",
                       "Het Contact",
                       "Low Contact",
                       "No Label"),
                   
                   # "group" = 
                   #   c(rep("Endoderm",7),rep("HFF",7))
                   "group" = c(
                     'TX',"TSS","Enh","Biv","ReprPC","Het","Low","NA",
                     'TX',"TSS","Enh","Biv","ReprPC","Het","Low","NA",
                     'TX',"TSS","Enh","Biv","ReprPC","Het","Low","NA",
                     'TX',"TSS","Enh","Biv","ReprPC","Het","Low","NA"
                   )
)


names(network) = c("source", "target", "value",'group')

# my_color <- 'd3.scaleOrdinal() .domain(["upreg", "downreg",
                                        # "TX","TSS","Enh","WkEnh","Biv","Het","Low"]) 

my_color <- 'd3.scaleOrdinal() .domain(["col",
                                        "TX","TSS","Enh","Biv","ReprPC","Het","Low","NA"]) 
                               .range(["rgba(93,93,93,.25)",
                                       "#00CD00","#FFFF00","#FF0000","#A020F0","#964B00","#000000","#F0F0F0"])'

# Plot
library(networkD3)
sn <- sankeyNetwork(Links = data.frame(network), Nodes = nodes,
                    NodeGroup = 'group',LinkGroup = 'group',
                    colourScale=my_color,
                    Source = "source", Target = "target",
                    Value = "value", NodeID = "name",fontSize= 14, nodeWidth = 7,iterations = 0)
sn
# sn <- htmlwidgets::prependContent(sn, htmltools::tags$h1("HERVH Candidate Loci"))
sn <- htmlwidgets::prependContent(sn, htmltools::HTML('<h3 style="font-family:Arial, sans-serif;font-weight:bold;margin-bottom:0;text-align:center">CIS Changes Between Cell Types</h3><div style="position:relative;padding: 10px 20px 0px 20px;margin-bottom:50px;">
<p style="font-family:Arial, sans-serif;font-weight:bold;position:absolute;left:0;">H1</p>
<p style="font-family:Arial, sans-serif;font-weight:bold;position:absolute;left:29%;">Endoderm</p>
<p style="font-family:Arial, sans-serif;font-weight:bold;position:absolute;right:32%;">HFF</p>
<p style="font-family:Arial, sans-serif;font-weight:bold;position:absolute;right:0;">HeLa</p>
</div>'))

# htmlwidgets::sizingPolicy(padding = 10, browser.fill = TRUE)
# sn$sizingPolicy$viewer$fill <- FALSE
sn$sizingPolicy$browser$fill <- F
sn$sizingPolicy$browser$defaultWidth <- 1100

# you save it as an html
saveNetwork(sn, "total_regions_200bp_window.html")

library(webshot2)
# you convert it as png
# HAVE TO OPEN GOOGLE CHROME
webshot2::webshot("total_regions_200bp_window.html",
                  "total_regions_200bp_window.png", vwidth =1100, vheight = 800,zoom = 2)





#### SUB SANKEY ChromHMM the same, labels change

cell_type_1 <- 'hff'
cell_type_2 <- 'hela'


final_net_ss <- final_net_cluster[final_net[[paste0(cell_type_1,"_base_region")]] == final_net[[paste0(cell_type_2,"_base_region")]] & 
                                    !is.na(final_net[[paste0(cell_type_1,"_base_region")]]) &
                                    !is.na(final_net[[paste0(cell_type_2,"_base_region")]]),]
final_net_ss <- final_net_ss[,c(3,4)]

network <- rbind.data.frame(table(final_net_ss))
colnames(network) <- c("source", "target", "value")

# colnames(network) <- c("target", "source", "value")
network$source <- as.numeric(network$source)
network$target <- as.numeric(network$target)
network$value <- as.numeric(network$value)

# # MUST BE ZERO-INDEXED
# network$target <- network$target + 7

network$group <- 'col'

#sm_net <- network[network$n > 200 & network$n < 500,]
# sm_net$target <- sm_net$target + 15


nodes = data.frame("name" = 
                     c("TX Contact",
                       "TSS Contact",
                       "Enhancer Contact",
                       "Bivalent Contact",
                       "ReprPC Contact",
                       "Het Contact",
                       "Low Contact",

                       "TX Contact",
                       "TSS Contact",
                       "Enhancer Contact",
                       "Bivalent Contact",
                       "ReprPC Contact",
                       "Het Contact",
                       "Low Contact"),
                       
                   
                   # "group" = 
                   #   c(rep("Endoderm",7),rep("HFF",7))
                   "group" = c(
                     'TX',"TSS","Enh","Biv","ReprPC","Het","Low",
                     'TX',"TSS","Enh","Biv","ReprPC","Het","Low"
                     
                   )
)


names(network) = c("source", "target", "value",'group')

network$target <- network$target - min(network[,1]) - 1  # NO NO LABEL
network$source <- network$source - min(network[,1])   # NO NO LABEL


# my_color <- 'd3.scaleOrdinal() .domain(["upreg", "downreg",
# "TX","TSS","Enh","WkEnh","Biv","Het","Low"]) 

my_color <- 'd3.scaleOrdinal() .domain(["col",
                                        "TX","TSS","Enh","Biv","ReprPC","Het","Low"]) 
                               .range(["rgba(93,93,93,.25)",
                                       "#00CD00","#FFFF00","#FF0000",
                                       "#A020F0","#964B00","#000000","#F0F0F0"])'






# Plot
library(networkD3)
sn <- sankeyNetwork(Links = data.frame(network), Nodes = nodes,
                    NodeGroup = 'group',LinkGroup = 'group',
                    colourScale=my_color,
                    Source = "source", Target = "target",
                    Value = "value", NodeID = "name",fontSize= 14, nodeWidth = 7,iterations = 0)
sn
# sn <- htmlwidgets::prependContent(sn, htmltools::tags$h1("HERVH Candidate Loci"))
sn <- htmlwidgets::prependContent(sn, htmltools::HTML('<h3 style="font-family:Arial, sans-serif;font-weight:bold;margin-bottom:0;text-align:center">HFF - HeLa CCS Changes with Constant ChromHMM label</h3><div style="position:relative;padding: 10px 20px 0px 20px;margin-bottom:50px;">
<p style="font-family:Arial, sans-serif;font-weight:bold;position:absolute;left:0;">HFF</p>
<p style="font-family:Arial, sans-serif;font-weight:bold;position:absolute;right:0;">HeLa</p>
</div>'))

# htmlwidgets::sizingPolicy(padding = 10, browser.fill = TRUE)
# sn$sizingPolicy$viewer$fill <- FALSE
sn$sizingPolicy$browser$fill <- F
sn$sizingPolicy$browser$defaultWidth <- 400

# you save it as an html
saveNetwork(sn, "constant_chromhmm/hff_hela_200bp_window.html")

library(webshot2)
# you convert it as png
# HAVE TO OPEN GOOGLE CHROME
webshot2::webshot("constant_chromhmm/hff_hela_200bp_window.html",
                  "constant_chromhmm/hff_hela_200bp_window.png", vwidth = 400, vheight = 800,zoom = 2)










# 
# 
# 
# 
# 
# 
# ############
# 
# 
# types = unique(total_int$`2_cell_type`)
# types <- types[!grepl('\\.',types)]
# 
# # nm <- c("1/TssBiv","2/EnhA2","3/TssFlnkU","4/TssA",
# #         "5/ReprPC","6/EnhWk1","7/TxWk","8/EnhA1","9/EnhBiv",
# #         "10/TssFlnkD",'11/EnhWk2',"12/TssFlnk")
# 
# order = c(2,4,6,10,12,8,7,5,13,9,14,1,11,3)
# h_list = list()
# for (t in types) {
#   
#   ss <- total_int[total_int$`2_cell_type`==t,]
#   
#   tabled <- table(ss$`1_cluster`,ss$`2_cluster`)
#   
#   tabled <- tabled[order,order]
#   
#   col_fun = colorRamp2(c(0, 1000), c("gray90", "red3"))
#   
#   h <- Heatmap(tabled,cluster_rows = F,cluster_columns = F,name = t,column_title = t,col=col_fun)
#   h_list[[t]] = h
# }
# 
# # png(paste0(cell_type,"_cluster_changes.png"),width = 2400,height = 900,res = 200)
# draw(h_list[[1]] + h_list[[2]] + h_list[[3]], column_title=paste0(cell_type," cluster changes by cell type"),
#                                                                   row_title = paste0(cell_type," cluster"))
# # dev.off()
# 
# 
# 
# 
# 
# # table(total_int$match)
# 
#  # 26% have same chromHMM, different CCS cluster
# sum(mismatched$`1_base_region` == mismatched$`2_base_region`) /
#   (nrow(mismatched))
# 
# sum(matched$`1_base_region` == matched$`2_base_region`) /
#   (nrow(matched))
# 
# sum(matched$`1_base_region` != matched$`2_base_region`) /
#   (nrow(matched))
# 
# sum(matched$`1_base_region` == matched$`2_base_region`) + sum(matched$`1_base_region` != matched$`2_base_region`)
# 
# mismatch_cluster__match_base_region <- mismatched[mismatched$`1_base_region` == mismatched$`2_base_region`,]
# 
# enha1 <- mismatch_cluster__match_base_region[grepl("EnhA",mismatch_cluster__match_base_region$`1_base_region`) & 
#                                                mismatch_cluster__match_base_region$`2_cell_type` == 'hela' &
#                                                mismatch_cluster__match_base_region$`1_cluster` == 5 &
#                                                mismatch_cluster__match_base_region$`2_cluster` == 9,]
# 
# 
# 
# table(enha1$`1_cluster`,enha1$`2_cluster`)
# 
# 
# total_regions <- nrow(cbind(total_int$`1_chr`,total_int$`1_pos1`,total_int$`1_pos2`))
# 
# 
# 
# 
# 
# # amount of changes
# matched <- (total_int[total_int$`1_cluster` == total_int$`2_cluster`,])
# mismatched <- (total_int[total_int$`1_cluster` != total_int$`2_cluster`,])
# 
# 
# 
# 
# #h1
# mismatched_h1 <- mismatched[mismatched$`2_cell_type` == "h1",]
# 
# # 40.7% have same chromHMM, different CCS cluster
# sum(mismatched_h1$`1_base_region` == mismatched_h1$`2_base_region`) /
#   (nrow(mismatched_h1))
# 
# #endo
# mismatched_endo <- mismatched[mismatched$`2_cell_type` == "endo",]
# 
# # 40.7% have same chromHMM, different CCS cluster
# sum(mismatched_endo$`1_base_region` == mismatched_endo$`2_base_region`) /
#   (nrow(mismatched_endo))
# 
# 
# #hff
# mismatched_hff <- mismatched[mismatched$`2_cell_type` == "hff",]
# 
# # 33.3% have same chromHMM, different CCS cluster
# sum(mismatched_hff$`1_base_region` == mismatched_hff$`2_base_region`) /
#   (nrow(mismatched_hff))
# 
# #hela
# mismatched_hela <- mismatched[mismatched$`2_cell_type` == "hela",]
# 
# # 28.7% have same chromHMM, different CCS cluster
# sum(mismatched_hela$`1_base_region` == mismatched_hela$`2_base_region`) /
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
