
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

library(RColorBrewer)
library(circlize)

# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(fpc)

# generate mat_raw data

# mat_raw <- read_tsv("matrix_connections.tsv.gz",col_names = T)
# 
# mat_raw <- as.data.frame(mat_raw)
# mat_raw$base_region <- factor(mat_raw$base_region)
# mat_raw[,3:20] <- apply(mat_raw[,3:20],2,as.numeric)
# # mat_raw$cell_type <- factor(ifelse(grepl('_endo',mat_raw[,1]),'endo',
# #                                ifelse(grepl('_hff',mat_raw[,1]),'hff',
# #                                       ifelse(grepl('_h1',mat_raw[,1]),'h1','hela')
# #                                       )
# #                                )
# #                             )
# # mat_raw <- mat_raw[apply(mat_raw[,3:20],1,function(x) sum(x!= 0)) != 0,]
# saveRDS(mat_raw,"mat_raw.rds")

mat_raw <- readRDS("mat_raw.rds")

data <- mat_raw[,c(3:9,11:20)] # remove quiescent col
not_only_quies <- apply(data,1,function(x) sum(x > 0)) > 1

data <- data[not_only_quies,] # get rid of rows with only quiescent contact
mat_raw_ss <- mat_raw[not_only_quies,]

# # # Elbow Method for finding the optimal number of clusters
# # # https://www.r-bloggers.com/2017/02/finding-optimal-number-of-clusters/

# set.seed(123)
# k.max <- 25
# 
# # wss <- sapply(10:k.max,
# #               function(k){kmeans(data, k, nstart=5)$tot.withinss})
# # 
# # 
# # takes a long time
# # png(paste0("elbow_plot_noquies.png"),width = 500,height = 500,res = 150)
# # plot(10:k.max, wss,
# #      type="b", pch = 19, frame = FALSE,
# #      xlab="Number of clusters K",
# #      ylab="Total within-clusters sum of squares",
# #      main="kmeans sum squared")
# # dev.off()


#### If looking for best clustering, loop to run kmeans, save models and heatmaps:
# ***  must create directories: models/  and  heatmaps/ for saving objects  ***


# for (i in 10:25) {
i = 18 # set up for 18 clusters

    mat <- data
    
    
    ### RUN kmeans
    clus <- kmeans(mat,centers = i,nstart = 100)
    saveRDS(clus,paste0("models/kmeans_k-",i,".rds"))
    
    # clus = readRDS(paste0("models/kmeans_k-",i,".rds"))
    
    
    mat$cluster <- factor(clus$cluster,levels = 0:max(clus$cluster))
    # mat$cell_type <- factor(ifelse(grepl('_endo',data[,1]),'endo',
    #                                ifelse(grepl('_hff',data[,1]),'hff',
    #                                       ifelse(grepl('_h1',data[,1]),'h1','hela')
    #                                )
    # )
    # )
    mat$base_region <- factor(mat_raw_ss$base_region,levels = 
                                c("EnhA1","EnhA2", "EnhBiv", "EnhG1", "EnhG2", "EnhWk", 
                                  "TssA","TssBiv","TssFlnk","TssFlnkD", "TssFlnkU", "Tx","TxWk",
                                  "ReprPC", "ReprPCWk","Het","ZNF-Rpts"))
    
    
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
                              values_from = freq,values_fill = 0,id_cols = -n)
    
    heatmap_df <- heatmap_df[,match(c("cluster","EnhA1","EnhA2", "EnhWk", "TxWk",
                                                  "TssA","TssFlnkU","TssFlnk","TssFlnkD", 
                                                  "EnhG1", "EnhG2", "Tx",
                                                  "EnhBiv", "TssBiv", "ReprPC", "ReprPCWk",
                                                  "Het","ZNF-Rpts"),colnames(heatmap_df))]
    heatmap_df <- heatmap_df[order(heatmap_df$cluster),]
    
    heatmap_mat <- as.matrix(heatmap_df[,-1])
    rownames(heatmap_mat) <- heatmap_df$cluster
   
   
    plot_regions <- gather(mat,
                           key = "contact_region_type",
                           value = 'score',-c('cluster',base_region))
    
    plot_regions$contact_region_type <- factor(plot_regions$contact_region_type,levels = 
                                                 c("EnhA1_contact","EnhA2_contact",
                                                   "EnhWk_contact", "TxWk_contact",
                                                   "TssA_contact","TssFlnkU_contact",
                                                   "TssFlnk_contact","TssFlnkD_contact",
                                                   "EnhG1_contact", "EnhG2_contact","Tx_contact",
                                                   "EnhBiv_contact","TssBiv_contact",
                                                   "ReprPC_contact", "ReprPCWk_contact",
                                                   "Het_contact", "Quies_contact","ZNF-Rpts_contact"),
    
                                               labels = c("EnhA1","EnhA2", "EnhWk", "TxWk",
                                                          "TssA","TssFlnkU","TssFlnk","TssFlnkD", 
                                                          "EnhG1", "EnhG2", "Tx",
                                                          "EnhBiv", "TssBiv", "ReprPC", "ReprPCWk",
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
    
    # clust_stats <- mat %>% group_by(cluster,cell_type) %>% 
    #   summarise(n = n()) %>%
    #   mutate(freq = n / sum(n))
    # 
    # clust_stats_df <- pivot_wider(clust_stats,names_from = cell_type,cluster, values_from = freq)
    # rownames(clust_stats_df) <- clust_stats_df$cluster
    # clust_stats_df <- clust_stats_df[,-1]
    # 
    # for cell type annotations
    # row_ann <- HeatmapAnnotation("% of Endo" = clust_stats_df$endo,
    #                              "% of H1" = clust_stats_df$h1,
    #                              "% of HFF" = clust_stats_df$hff,
    #                              "% of HeLa" = clust_stats_df$hela,
    #                              
    #                              col = list("% of Endo" = colorRamp2(c(0, 1), c("white",
    #                                                                             "palevioletred")),
    #                                         "% of H1" = colorRamp2(c(0, 1), c("white", "skyblue3")),
    #                                         "% of HFF" = colorRamp2(c(0, 1), c("white", "orangered2")),
    #                                         "% of HeLa" = colorRamp2(c(0, 1), c("white", "olivedrab4"))),
    #                              
    #                              which = 'row')
    
    library(circlize)
    col_fun = colorRamp2(c(0, 0.3), c("gray95","darkgoldenrod3"))
    col_fun2 = colorRamp2(c(-5, 0, 5), c("cornflowerblue","white","red"))
    
    contact = Heatmap(log2(med_mat +1e-20),cluster_rows = T,cluster_columns = F,col = col_fun2,
                      column_title = "Contact Region",
                      row_title = "Cluster",
                      heatmap_legend_param = list(title = 'log2 ratio\nof contact'),name = "median",
                      # show_row_dend = T,clustering_method_rows = "complete",
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%.1f", log2(med_mat[i, j])), x, y, gp = gpar(fontsize = 9))
                        })
    
    base = Heatmap(heatmap_mat,cluster_rows = F,cluster_columns = F,col = col_fun,
                   column_title = "Base Region",
                   # right_annotation = row_ann,
                   heatmap_legend_param = list(title = 'state %'), name = 'base',
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if (heatmap_mat[i, j]*100 >= 1){
                     grid.text(sprintf("%.0f%%", heatmap_mat[i, j]*100), x, y, gp = gpar(fontsize = 9))
                   }
                   })
    
    png(paste0("heatmaps/kmeans_k-",i,".png"),width = 3000,height = 1100,res = 200)
    draw(contact + base)
    dev.off()
    
    
    # png(paste0("kmeans_k-",i,"just_contact.png"),width = 1700,height = 1100,res = 200)
    # draw(contact)
    # dev.off()
    
  # }




# GET LABELS:

con_mat <- mat_raw %>% separate(region,c("chr","type","pos1","pos2"),sep = '[:_-]')

# scaled_h <- apply(model$h, 2, function(x) x/sum(x))

all_labelled <- data.frame(type = con_mat$type,
                           chr = con_mat$chr,
                           pos1 = con_mat$pos1,
                           pos2 = con_mat$pos2,
                           base_region = con_mat$base_region,
                           cluster = mat$cluster)

write_tsv(all_labelled,"clustered_SCC_matrix.tsv.gz")
