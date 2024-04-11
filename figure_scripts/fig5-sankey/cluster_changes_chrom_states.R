

setwd("~/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig5-all_sankey/create_bed_network_files/")

library(tidyverse)
library(ggsignif)
library(circlize)
library(ComplexHeatmap)

final_net <- readRDS("final_net.rds")

final_net_cluster <- final_net %>% select(h1_base_region, endo_base_region, hff_base_region,hela_base_region)


# 6 = "No Label"
# 5 repressed  = "ReprPCWk","Het","ZNF-Rpts","ReprPC"
# 4 bivalent = "EnhBiv", "TssBiv"
# 3 weak enhancer = "EnhWk" 
# 2 enhancer = "EnhA1"  "EnhA2"  "EnhG1" "EnhG2" 
# 1 TSS = "TssFlnk" "TssA"   "TssFlnkU" "TssFlnkD"
# 0 Tx =   "Tx"  "TxWk" 

final_net_cluster$h1_base_region <- ifelse(
  is.na(final_net_cluster$h1_base_region),6,ifelse(
  final_net_cluster$h1_base_region %in% c("ReprPCWk","Het","ZNF-Rpts","ReprPC"),5,ifelse(
  final_net_cluster$h1_base_region %in% c("EnhBiv", "TssBiv"),4,ifelse(
    final_net_cluster$h1_base_region == "EnhWk",3, ifelse(
      final_net_cluster$h1_base_region %in% c("EnhA1","EnhA2","EnhG1","EnhG2"),2, ifelse(
        final_net_cluster$h1_base_region %in% c("TssFlnk", "TssA","TssFlnkU", "TssFlnkD"),1,0
      )
    )
  )
)))

final_net_cluster$endo_base_region <- ifelse(
  is.na(final_net_cluster$endo_base_region),6,ifelse(
    final_net_cluster$endo_base_region %in% c("ReprPCWk","Het","ZNF-Rpts","ReprPC"),5,ifelse(
      final_net_cluster$endo_base_region %in% c("EnhBiv", "TssBiv"),4,ifelse(
        final_net_cluster$endo_base_region == "EnhWk",3, ifelse(
          final_net_cluster$endo_base_region %in% c("EnhA1","EnhA2","EnhG1","EnhG2"),2, ifelse(
            final_net_cluster$endo_base_region %in% c("TssFlnk", "TssA","TssFlnkU", "TssFlnkD"),1,0
          )
        )
      )
    )))

final_net_cluster$hff_base_region <- ifelse(
  is.na(final_net_cluster$hff_base_region),6,ifelse(
    final_net_cluster$hff_base_region %in% c("ReprPCWk","Het","ZNF-Rpts","ReprPC"),5,ifelse(
      final_net_cluster$hff_base_region %in% c("EnhBiv", "TssBiv"),4,ifelse(
        final_net_cluster$hff_base_region == "EnhWk",3, ifelse(
          final_net_cluster$hff_base_region %in% c("EnhA1","EnhA2","EnhG1","EnhG2"),2, ifelse(
            final_net_cluster$hff_base_region %in% c("TssFlnk", "TssA","TssFlnkU", "TssFlnkD"),1,0
          )
        )
      )
    )))


final_net_cluster$hela_base_region <- ifelse(
  is.na(final_net_cluster$hela_base_region),6,ifelse(
    final_net_cluster$hela_base_region %in% c("ReprPCWk","Het","ZNF-Rpts","ReprPC"),5,ifelse(
      final_net_cluster$hela_base_region %in% c("EnhBiv", "TssBiv"),4,ifelse(
        final_net_cluster$hela_base_region == "EnhWk",3, ifelse(
          final_net_cluster$hela_base_region %in% c("EnhA1","EnhA2","EnhG1","EnhG2"),2, ifelse(
            final_net_cluster$hela_base_region %in% c("TssFlnk", "TssA","TssFlnkU", "TssFlnkD"),1,0
          )
        )
      )
    )))


# final_net_cluster$h1_cluster <- factor(final_net_cluster$h1_cluster,levels = c(0,1,2,3,4,5,6))

final_net_cluster$endo_base_region <- final_net_cluster$endo_base_region + 7
# final_net_cluster$endo_cluster <- factor(final_net_cluster$endo_cluster,levels = c(7,8,9,10,11,12,13))


final_net_cluster$hff_base_region <- final_net_cluster$hff_base_region + 14
# final_net_cluster$hff_cluster <- factor(final_net_cluster$hff_cluster,levels = c(14,15,16,17,18,19,20))


final_net_cluster$hela_base_region <- final_net_cluster$hela_base_region + 21
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

# 5 repressed  = "ReprPCWk","Het","ZNF-Rpts","ReprPC"
# 4 bivalent = "EnhBiv", "TssBiv"
# 3 weak enhancer = "EnhWk" 
# 2 enhancer = "EnhA1"  "EnhA2"  "EnhG1" "EnhG2"  # CHANGE THIS
# 1 TSS = "TssFlnk" "TssA"   "TssFlnkU" "TssFlnkD"
# 0 Tx =   "Tx"  "TxWk" 

nodes = data.frame("name" = 
                     c("TX",
                       "TSS",
                       "Enhancer",
                       "Weak Enhancer",
                       "Bivalent",
                       "Repressed",
                       "No Label",
                       
                       "TX",
                       "TSS",
                       "Enhancer",
                       "Weak Enhancer",
                       "Bivalent",
                       "Repressed",
                       "No Label",
                       
                       
                       "TX",
                       "TSS",
                       "Enhancer",
                       "Weak Enhancer",
                       "Bivalent",
                       "Repressed",
                       "No Label",
                       
                       
                       "TX",
                       "TSS",
                       "Enhancer",
                       "Weak Enhancer",
                       "Bivalent",
                       "Repressed",
                       "No Label"),
                   
                   # "group" = 
                   #   c(rep("Endoderm",7),rep("HFF",7))
                   "group" = c(
                     'TX',"TSS","Enh","WkEnh","Biv","Het","NA",
                     'TX',"TSS","Enh","WkEnh","Biv","Het","NA",
                     'TX',"TSS","Enh","WkEnh","Biv","Het","NA",
                     'TX',"TSS","Enh","WkEnh","Biv","Het","NA"
                   )
)


names(network) = c("source", "target", "value",'group')

# my_color <- 'd3.scaleOrdinal() .domain(["upreg", "downreg",
                                        # "TX","TSS","Enh","WkEnh","Biv","Het","Low"]) 

my_color <- 'd3.scaleOrdinal() .domain(["col",
                                        "TX","TSS","Enh","WkEnh","Biv","Het","Low","NA"]) 
                               .range(["rgba(93,93,93,.25)",
                                       "#00CD00","#FFFF00","#FF0000","#FFC1C1","#A020F0","#000000","#F0F0F0"])'

# Plot
sn <- sankeyNetwork(Links = data.frame(network), Nodes = nodes,
                    NodeGroup = 'group',LinkGroup = 'group',
                    colourScale=my_color,
                    Source = "source", Target = "target",
                    Value = "value", NodeID = "name",fontSize= 14, nodeWidth = 7,iterations = 0)

sn
# sn <- htmlwidgets::prependContent(sn, htmltools::tags$h1("HERVH Candidate Loci"))
sn <- htmlwidgets::prependContent(sn, htmltools::HTML('<h3 style="font-family:Arial, sans-serif;font-weight:bold;margin-bottom:0;text-align:center">Chromatin State Changes Between Cell Types</h3><div style="position:relative;padding: 10px 20px 0px 20px;margin-bottom:50px;">
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
saveNetwork(sn, "total_regions_chromatin_200bp_window.html")

library(webshot2)
# you convert it as png
# HAVE TO OPEN GOOGLE CHROME
webshot2::webshot("total_regions_chromatin_200bp_window.html",
                  "total_regions_chromatin_200bp_window.png", vwidth =1100, vheight = 800,zoom = 2)

