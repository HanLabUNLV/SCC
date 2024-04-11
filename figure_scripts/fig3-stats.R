# connection stats

setwd("/Users/coripenrod/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig3-raw_heatmaps/")

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms

library(clustMixType)
library(RColorBrewer)
library(circlize)

# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)


# mat_raw <- read_tsv("../matrix_connections.tsv.gz",col_names = T)
# saveRDS(mat_raw,"../data/mat_raw.rds")

mat_raw <- readRDS("../data/mat_raw_mean.rds")
# 
# bed_file <- mat %>% separate(region,into = c("chr","cell_type","pos1","pos2"),remove = F) %>% select("chr","pos1","pos2","cell_type","region")
# write_tsv(bed_file,"../mat.bed",col_names = F)

# 
# library(circlize)
# col_fun = colorRamp2(c(0, 0.5), c("gray95","darkgoldenrod3"))
# col_fun2 = colorRamp2(c(0, 0.05), c("gray95","blue4"))

mat_raw$cell_type <- factor(ifelse(grepl('_endo',mat_raw[,1]),'Endoderm',
                                                                  ifelse(grepl('_hff',mat_raw[,1]),'HFF',
                                                                         ifelse(grepl('_h1',mat_raw[,1]),'H1','HeLa')
                                                                         )
                                                                  )
                                                               )
mat_raw$base_region <- factor(mat_raw$base_region, levels = c("EnhA1","EnhA2","EnhBiv", "EnhG1",
                                                              "EnhG2", "EnhWk","Het", "Quies",
                                                              "ReprPC", "ReprPCWk",
                                                              "TssA", "TssBiv", 
                                                              "TssFlnkD","TssFlnk","TssFlnkU", 
                                                              "Tx", "TxWk", "ZNF-Rpts"))


row_ann <- HeatmapAnnotation("Base Label" = mat_raw$base_region, 
                             "Cell Type" = mat_raw$cell_type,
                             col = list("Base Label" = c("EnhA1"="red","EnhA2" = "red4",
                                                         "EnhBiv" = "purple", "EnhG1" = "palevioletred2",
                                                         "EnhG2" = "palevioletred", "EnhWk" = "rosybrown1",
                                                         "Het" = "gray90", "Quies" = "gray48","ReprPC" = "gray", "ReprPCWk" = "gray60",
                                                         "TssA" = "yellow", "TssBiv" = "plum", 
                                                         "TssFlnkD" ="lightyellow","TssFlnk" = "orange","TssFlnkU" = "khaki", 
                                                         "Tx" = "green3", "TxWk" = "lightgreen", "ZNF-Rpts" = "black"),
                                          "Cell Type" = c("H1" = "lightskyblue",
                                                        "Endoderm" = "darkgoldenrod1",
                                                        "HeLa" = "salmon",
                                                        "HFF" = "mediumpurple1")),
                             which = 'row')
# install.packages("seriation")
library(seriation)
# list_seriation_methods("matrix")
# o = seriate(as.matrix(mat_raw[,3:20]), method = "PCA")
col_fun = colorRamp2(c(0,0.05, 0.1), c("mediumblue","white","firebrick1"))

colnames(mat_raw) <- sub('_contact','',colnames(mat_raw))


png("raw_mean.png",width = 1800,height = 1600,res = 250)
Heatmap(as.matrix(mat_raw[,3:20]),cluster_rows = F,cluster_columns = F, col = col_fun,
        column_title = "Summed contact values for single regions annotated by ChromHMM",
                  row_title = "Labelled Regions", left_annotation = row_ann,
                  # row_order = get_order(o,1), column_order = get_order(o,2),
                  row_order = order(mat_raw$base_region,
                                    mat_raw$cell_type),
                  heatmap_legend_param = list(title = 'Sum Contact Score'),name = "contact",
                  show_row_dend = F,show_column_names = T,show_column_dend = F,show_row_names = F,
                  use_raster = F)
dev.off()




# Degree of connectedness

degree <- apply(mat_raw[,3:20],1,function(x) sum(x!= 0))

connectedness <- data.frame(degree = degree, cell_type = mat_raw$cell_type,base_region = mat_raw$base_region)

summary(connectedness$degree[connectedness$cell_type == 'H1'])
summary(connectedness$degree[connectedness$cell_type == 'HeLa'])
summary(connectedness$degree[connectedness$cell_type == 'Endoderm'])
summary(connectedness$degree[connectedness$cell_type == 'HFF'])

summary(connectedness$degree[connectedness$cell_type == 'h1'])
summary(connectedness$degree[connectedness$cell_type == 'hela'])
summary(connectedness$degree[connectedness$cell_type == 'endo'])
summary(connectedness$degree[connectedness$cell_type == 'hff'])


png("connectedness_mean.png",width = 1400,height = 500,res = 225)
ggplot(connectedness,aes(x = degree,fill = cell_type)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~cell_type,nrow = 1) +
  labs(title = "# of Unique ChromHMM Labels in Contact with Single Region",
       x = "# of unique ChromHMM labels",
       y = "Count of Regions") + theme(legend.position="none")
dev.off()

