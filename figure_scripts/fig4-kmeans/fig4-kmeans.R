
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

library(RColorBrewer)
library(circlize)

library(ComplexHeatmap)
library(fpc)


# mat_raw <- read_tsv("../data/matrix_connections_mean.tsv.gz",col_names = T)
# 
# mat_raw <- as.data.frame(mat_raw)
# mat_raw$base_region <- factor(mat_raw$base_region)
# mat_raw[,3:20] <- apply(mat_raw[,3:20],2,as.numeric)
# mat_raw$cell_type <- factor(ifelse(grepl('_endo',mat_raw[,1]),'endo',
#                                ifelse(grepl('_hff',mat_raw[,1]),'hff',
#                                       ifelse(grepl('_h1',mat_raw[,1]),'h1','hela')
#                                       )
#                                )
#                             )
# # mat_raw <- mat_raw[apply(mat_raw[,3:20],1,function(x) sum(x!= 0)) != 0,]
# saveRDS(mat_raw,"../data/mat_raw_mean.rds")

mat_raw <- readRDS("../data/mat_raw_mean.rds")

 
data <- mat_raw[,c(3:9,11:20)]
data <- data[apply(data,1,function(x) sum(x > 0)) > 1,] # get rid of only Quies contact

mat_raw_ss <- mat_raw[apply(mat_raw[,c(3:9,11:20)],1,function(x) sum(x > 0)) > 1,]


i = 18
mat <- data

clus = readRDS(paste0("kmeans_k-",i,".rds"))

mat$cluster <- factor(clus$cluster,levels = 0:max(clus$cluster))
mat$cell_type <- factor(ifelse(grepl('_endo',mat_raw_ss[,1]),'endo',
                               ifelse(grepl('_hff',mat_raw_ss[,1]),'hff',
                                      ifelse(grepl('_h1',mat_raw_ss[,1]),'h1','hela')
                               )
)
)
mat$base_region <- factor(mat_raw_ss$base_region,levels = 
                            c("EnhA1","EnhA2", "EnhBiv", "EnhG1", "EnhG2", "EnhWk", 
                              "TssA","TssBiv","TssFlnk","TssFlnkD", "TssFlnkU", "Tx","TxWk",
                              "ReprPC", "ReprPCWk","Het","ZNF-Rpts"))

mat$cluster <- factor(mat$cluster, levels = as.character(
  c(11,8,4,12,13,7,6,16,10,2,3,17,18,1,5,9,15,14)),
                        labels =
                        c(`11` = "cont_Tss_Enh",
                            `8` = "cont_Tss_EnhG",
                            `4` = "cont_Tss_noEnh",
                            `12` = "cont_Tss_EnhWk",
                            `13` = "cont_TxWk_Enh",
                            `7` = "cont_TxWk_EnhWk",
                            `6` = "cont_Tx_EnhG1",
                            `16` = "cont_Tx_EnhG2",
                            `10` = "cont_EnhA",
                            `2` = "cont_EnhA_EnhWk",
                            `3` = "cont_Biv",
                            `17` = "cont_ReprPC_Biv",
                            `18` = "cont_ReprPC",
                            `1` = "cont_Low_TssBiv",
                            `5` = "cont_Low_EnhWk",
                            `9` = "cont_Het_Rpts",
                            `15` = "cont_Het_Rpts_strong",
                            `14` ="cont_Het_Rpts_strongest"))



class_clus_freq <- mat %>% group_by(base_region,cluster) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

class_clus_freq$base_region <- factor(class_clus_freq$base_region,levels = 
                                        c("EnhA1","EnhA2", "EnhWk", "TxWk",
                                          "TssA","TssFlnkU","TssFlnk","TssFlnkD", 
                                          "EnhG1", "EnhG2", "Tx",
                                          "EnhBiv", "TssBiv", "ReprPC", "ReprPCWk",
                                          "Het","ZNF-Rpts"))


heatmap_df <- pivot_wider(class_clus_freq,names_from = base_region,
                          values_from = freq,values_fill = 0,-n)

heatmap_df <- heatmap_df[,match(c("cluster","EnhA1","EnhA2", "EnhWk", "TxWk",
                                              "TssA","TssFlnkU","TssFlnk","TssFlnkD", 
                                              "EnhG1", "EnhG2", "Tx",
                                              "EnhBiv", "TssBiv", "ReprPC", "ReprPCWk",
                                              "Het","ZNF-Rpts"),colnames(heatmap_df))]
heatmap_df <- heatmap_df[order(heatmap_df$cluster),]

heatmap_mat <- as.matrix(heatmap_df[,-1])
rownames(heatmap_mat) <- heatmap_df$cluster


plot_regions <- gather(mat[,1:19],
                       key = "contact_region_type",
                       value = 'score',-c('cell_type','cluster'))

plot_regions$contact_region_type <- factor(plot_regions$contact_region_type,levels = 
                                             c("EnhA1_contact","EnhA2_contact",
                                               "EnhWk_contact", "TxWk_contact",
                                               "TssA_contact","TssFlnkU_contact",
                                               "TssFlnk_contact","TssFlnkD_contact",
                                               "EnhG1_contact", "EnhG2_contact","Tx_contact",
                                               "EnhBiv_contact","TssBiv_contact",
                                               "ReprPC_contact", "ReprPCWk_contact",
                                               # "Het_contact","ZNF-Rpts_contact"),
                                           "Het_contact", "Quies_contact","ZNF-Rpts_contact"),

                                           labels = c("EnhA1","EnhA2", "EnhWk", "TxWk",
                                                      "TssA","TssFlnkU","TssFlnk","TssFlnkD", 
                                                      "EnhG1", "EnhG2", "Tx",
                                                      "EnhBiv", "TssBiv", "ReprPC", "ReprPCWk",
                                                      # "Het","ZNF-Rpts"))
                                                      
                                                      "Het", "Quies","ZNF-Rpts"))


total_reg <- plot_regions %>% group_by(contact_region_type) %>% 
  summarise(mean_type_score = mean(score))
              
plot_regions_med <- plot_regions %>% group_by(contact_region_type,cluster) %>% 
  summarise(mean_cluster = mean(score)) %>% left_join(.,total_reg,by = "contact_region_type") %>%
  mutate(score = mean_cluster/mean_type_score) %>% select(contact_region_type,score,cluster)



med_df <- pivot_wider(plot_regions_med,names_from = contact_region_type,
                      values_from = score,values_fill = NA)
med_mat <- as.matrix(med_df[,-1])
rownames(med_mat) <- med_df$cluster

# annotations

clust_stats <- mat %>% group_by(cluster,cell_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

clust_stats_df <- pivot_wider(clust_stats,names_from = cell_type,cluster, values_from = freq)
rownames(clust_stats_df) <- clust_stats_df$cluster
clust_stats_df <- clust_stats_df[,-1]
row_ann <- HeatmapAnnotation("% of Endo" = clust_stats_df$endo,
                             "% of H1" = clust_stats_df$h1,
                             "% of HFF" = clust_stats_df$hff,
                             "% of HeLa" = clust_stats_df$hela,
                             
                             col = list("% of Endo" = colorRamp2(c(0, 1), c("white",
                                                                            "palevioletred")),
                                        "% of H1" = colorRamp2(c(0, 1), c("white", "skyblue3")),
                                        "% of HFF" = colorRamp2(c(0, 1), c("white", "orangered2")),
                                        "% of HeLa" = colorRamp2(c(0, 1), c("white", "olivedrab4"))),
                             
                             which = 'row')

library(circlize)
col_fun = colorRamp2(c(0, 0.3), c("gray95","darkgoldenrod3"))
col_fun2 = colorRamp2(c(-2, 0, 2), c("cornflowerblue","white","red"))

contact = Heatmap(log2(med_mat),cluster_rows = T,cluster_columns = F,col = col_fun2,
                  column_title = "Contact Region",
                  row_title = "Cluster",
                  heatmap_legend_param = list(title = 'log2 ratio\nof contact'),name = "median",
                  show_row_dend = T,clustering_method_rows = "complete",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.1f", log2(med_mat[i, j])), x, y, gp = gpar(fontsize = 9))
                    })

base = Heatmap(heatmap_mat,cluster_rows = F,cluster_columns = F,col = col_fun,
               column_title = "Base Region",right_annotation = row_ann,
               heatmap_legend_param = list(title = 'state %'), name = 'base',
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (heatmap_mat[i, j]*100 >= 1){
                 grid.text(sprintf("%.0f%%", heatmap_mat[i, j]*100), x, y, gp = gpar(fontsize = 9))
               }
               })

png(paste0("kmeans_k-",i,".png"),width = 3000,height = 1100,res = 200)
draw(contact + base)
dev.off()


png(paste0("kmeans_k-",i,"just_contact.png"),width = 1700,height = 1100,res = 200)
draw(contact)
dev.off()



# GET LABELS:

con_mat <- mat_raw_ss %>% separate(region,c("chr","type","pos1","pos2"),sep = '[:_-]')

# scaled_h <- apply(model$h, 2, function(x) x/sum(x))

all_labelled <- data.frame(type = con_mat$type,
                           chr = con_mat$chr,
                           pos1 = con_mat$pos1,
                           pos2 = con_mat$pos2,
                           base_region = con_mat$base_region,
                           cluster = mat$cluster)

# all_labelled <- all_labelled[all_labelled$top_region > 0.75,]

h1_labelled <- subset(all_labelled, type == 'h1')[,c(2:6)]
h1_labelled$cell_type <- 'h1'

endo_labelled <- subset(all_labelled, type == 'endo')[,c(2:6)]
endo_labelled$cell_type <- 'endo'

hela_labelled <- subset(all_labelled, type == 'hela')[,c(2:6)]
hela_labelled$cell_type <- 'hela'

hff_labelled <- subset(all_labelled, type == 'hff')[,c(2:6)]
hff_labelled$cell_type <- 'hff'


write_tsv(all_labelled,"../data/clustered_SCC_matrix.tsv.gz")

write_tsv(h1_labelled, "../fig5-all_sankey/create_bed_network_files/h1_labelled.bed", col_names = F)
write_tsv(endo_labelled, "../fig5-all_sankey/create_bed_network_files/endo_labelled.bed", col_names = F)
write_tsv(hela_labelled, "../fig5-all_sankey/create_bed_network_files/hela_labelled.bed", col_names = F)
write_tsv(hff_labelled, "../fig5-all_sankey/create_bed_network_files/hff_labelled.bed", col_names = F)

