
# Fig 5 heatmap of TSS labels by gene
setwd("~/Documents/UNLV/Year5/new_kmeans_majority_chromhmm_in_window/fig5-gene_sankey/")

library(readr)
library(ComplexHeatmap)
library(dplyr)

ccs_colors = c("Tss_noEnh_cont" ="#FFFF00",
               "Tss_Enh_cont"="#FFFF00",
               "TssFlnkU_cont"="#FFFF00",
               "TssFlnkD_cont"="#FFFF00",
               
               "EnhWk_cont2"="#FF0000",
               "EnhA1_cont"="#FF0000",
               "Enh_cont"= "#FF0000",
               "EnhWk_cont1" = "#FF0000",
               
               "Tx_cont1" = "#00CD00",
               "Tx_cont2" = "#00CD00",
               
               "Biv_cont" = "#FFC1C1",
               "ReprPC_cont"= "#FFC1C1",
               "Het_cont1"= "#A020F0",
               "Het_cont2" = "#A020F0",
               
               "Low_cont" = "#000000",
               "NA" = 'white'
)

     
chrom_colors = c("TssA" ="#FFFF00",
                 "TssFlnk"="#FFFF00",
                 "TssFlnkD"="#FFFF00",
                 "TssFlnkU"="#FFFF00",
                 
                 "EnhA1"="#FF0000",
                 "EnhA2"="#FF0000",
                 "EnhWk"= "#FF0000",
                 "EnhG1" = "#FF0000",
                 "EnhG2" = "#FF0000",
                 
                 "Tx" = "#00CD00",
                 "TxWk" = "#00CD00",
                 
                 "TssBiv" = "#FFC1C1",
                 "EnhBiv"= "#FFC1C1",
                 
                 "ReprPCWk"= "#FFC1C1",
                 "ReprPC" = "#FFC1C1",
                 
                 "ZNF-Rpts" = "#A020F0",
                 "Het" = "#A020F0",
                 
                 "NA" = 'white'
)

# GENE ORDER

h1_genes <- read_tsv("../fig8-enrichment/expressed_genes/h1_pone.0192625.s007.tsv",col_types = 'ccccccccccccccccccccc') %>% 
        select(HUES1_TPM,ensembl_id,hgnc_symbol) %>% subset(!is.na(hgnc_symbol))
colnames(h1_genes) <- c('expr','ensembl_id','hgnc_symbol')
h1_genes$expr <- as.numeric(h1_genes$expr)
gene_order_h1 <- h1_genes$hgnc_symbol[order(h1_genes$expr,decreasing = T)]


endo_genes <- read_tsv("../fig8-enrichment/expressed_genes/endo_Table_S8C-Table_1.tsv",col_types = 'cccccccc') %>%
        select(AveExpr,Gene,`Gene Name`)
colnames(endo_genes) <- c('expr','ensembl_id','hgnc_symbol')
endo_genes <- endo_genes %>% subset(!is.na(hgnc_symbol))
endo_genes$expr <- as.numeric(endo_genes$expr)
gene_order_endo <- endo_genes$hgnc_symbol[order(endo_genes$expr,decreasing = T)]


hff_genes <- read_tsv("../fig8-enrichment/expressed_genes/hff_GSM2448852_0h_3-aligned.tsv",col_types = 'cccccc') %>% 
        select(`0h_3_FPKM`,gene_id,gene_name)
colnames(hff_genes) <- c('expr','ensembl_id','hgnc_symbol')
hff_genes <- hff_genes %>% subset(!is.na(hgnc_symbol))
hff_genes$expr <- as.numeric(hff_genes$expr)
gene_order_hff <- hff_genes$hgnc_symbol[order(hff_genes$expr,decreasing = T)]


###


h1_ccs <- read_tsv("h1_tss_overlap.ccs_500bp.bed",col_names = F,na = character())

h1_ccs_mat <- as.matrix(h1_ccs[,-1])
rownames(h1_ccs_mat) <- h1_ccs$X1

# all_gene_h1 <- c(gene_order_h1,rownames(h1_ccs_mat))
all_gene_h1 <- gene_order_h1

all_gene_h1 <- all_gene_h1[!duplicated(all_gene_h1)]
all_gene_h1 <- all_gene_h1[all_gene_h1 %in% rownames(h1_ccs_mat)]

ss <- h1_ccs_mat[all_gene_h1,]
# png("h1_ccs_tss_heatmap.png",height = 1500,width = 900,res = 200)
h1_hm_ccs <- Heatmap(ss,use_raster = F,cluster_columns = F,cluster_rows = T,col = ccs_colors,
        show_heatmap_legend = T,show_row_names = F,show_column_names = F,row_title = "Genes",
        column_title = "H1 CCS labels",name = "CCS Label",
        row_order = all_gene_h1)
# dev.off()



h1_chrom <- read_tsv("h1_tss_overlap.chrom_500bp.bed",col_names = F,na = character())

h1_chrom_mat <- as.matrix(h1_chrom[,-1])
rownames(h1_chrom_mat) <- h1_chrom$X1

ss <- h1_chrom_mat[all_gene_h1,]

# png("h1_chrom_tss_heatmap.png",height = 1500,width = 900,res = 200)
h1_hm_chrom <- Heatmap(ss,use_raster = F,cluster_columns = F,cluster_rows = T,col = chrom_colors,
        show_heatmap_legend = T,show_row_names = F,show_column_names = F,row_title = "Genes",
        column_title = "H1 ChromHMM labels",name ="ChromHMM Label",
        row_order = all_gene_h1)
# dev.off()

png("h1_tss_heatmap.png",height = 1800,width = 1600,res = 200)
draw(h1_hm_ccs + h1_hm_chrom)
dev.off()

##########

# 
# 
# h1_ccs <- read_tsv("hela_tss_overlap.ccs_500bp.bed",col_names = F,na = character())
# 
# h1_ccs_mat <- as.matrix(h1_ccs[,-1])
# rownames(h1_ccs_mat) <- h1_ccs$X1
# 
# all_gene_h1 <- c(gene_order_hela,rownames(h1_ccs_mat))
# all_gene_h1 <- all_gene_h1[!duplicated(all_gene_h1)]
# all_gene_h1 <- all_gene_h1[all_gene_h1 %in% rownames(h1_ccs_mat)]
# 
# # png("h1_ccs_tss_heatmap.png",height = 1500,width = 900,res = 200)
# h1_hm_ccs <- Heatmap(h1_ccs_mat,use_raster = F,cluster_columns = F,cluster_rows = T,col = ccs_colors,
#                      show_heatmap_legend = T,show_row_names = F,show_column_names = F,row_title = "Genes",
#                      column_title = "HeLa CCS labels",name = "CCS Label",
#                      row_order = all_gene_h1)
# # dev.off()
# 
# 
# 
# h1_chrom <- read_tsv("hela_tss_overlap.chrom_500bp.bed",col_names = F,na = character())
# 
# h1_chrom_mat <- as.matrix(h1_chrom[,-1])
# rownames(h1_chrom_mat) <- h1_chrom$X1
# 
# # png("h1_chrom_tss_heatmap.png",height = 1500,width = 900,res = 200)
# h1_hm_chrom <- Heatmap(h1_chrom_mat,use_raster = F,cluster_columns = F,cluster_rows = T,col = chrom_colors,
#                        show_heatmap_legend = T,show_row_names = F,show_column_names = F,row_title = "Genes",
#                        column_title = "HeLa ChromHMM labels",name ="ChromHMM Label",
#                        row_order = all_gene_h1)
# # dev.off()
# 
# png("hela_tss_heatmap.png",height = 1800,width = 1600,res = 200)
# draw(h1_hm_ccs + h1_hm_chrom)
# dev.off()
# 
# 



############

h1_ccs <- read_tsv("hff_tss_overlap.ccs_500bp.bed",col_names = F,na = character())

h1_ccs_mat <- as.matrix(h1_ccs[,-1])
rownames(h1_ccs_mat) <- h1_ccs$X1

# all_gene_h1 <- c(gene_order_hff,rownames(h1_ccs_mat))
all_gene_h1 <- c(gene_order_hff)
all_gene_h1 <- all_gene_h1[!duplicated(all_gene_h1)]
all_gene_h1 <- all_gene_h1[all_gene_h1 %in% rownames(h1_ccs_mat)]


ss <- h1_ccs_mat[all_gene_h1,]
# png("h1_ccs_tss_heatmap.png",height = 1500,width = 900,res = 200)
h1_hm_ccs <- Heatmap(ss,use_raster = F,cluster_columns = F,cluster_rows = T,col = ccs_colors,
                     show_heatmap_legend = T,show_row_names = F,show_column_names = F,row_title = "Genes",
                     column_title = "HFF CCS labels",name = "CCS Label",
                     row_order = all_gene_h1)
# dev.off()



h1_chrom <- read_tsv("hff_tss_overlap.chrom_500bp.bed",col_names = F,na = character())

h1_chrom_mat <- as.matrix(h1_chrom[,-1])
rownames(h1_chrom_mat) <- h1_chrom$X1

ss <- h1_chrom_mat[all_gene_h1,]

# png("h1_chrom_tss_heatmap.png",height = 1500,width = 900,res = 200)
h1_hm_chrom <- Heatmap(ss,use_raster = F,cluster_columns = F,cluster_rows = T,col = chrom_colors,
                       show_heatmap_legend = T,show_row_names = F,show_column_names = F,row_title = "Genes",
                       column_title = "HFF ChromHMM labels",name ="ChromHMM Label",
                       row_order = all_gene_h1)
# dev.off()

png("hff_tss_heatmap.png",height = 1800,width = 1600,res = 200)
draw(h1_hm_ccs + h1_hm_chrom)
dev.off()


##############

h1_ccs <- read_tsv("endo_tss_overlap.ccs_500bp.bed",col_names = F,na = character())

h1_ccs_mat <- as.matrix(h1_ccs[,-1])
rownames(h1_ccs_mat) <- h1_ccs$X1


# all_gene_h1 <- c(gene_order_endo,rownames(h1_ccs_mat))
all_gene_h1 <- gene_order_endo

all_gene_h1 <- all_gene_h1[!duplicated(all_gene_h1)]
all_gene_h1 <- all_gene_h1[all_gene_h1 %in% rownames(h1_ccs_mat)]

ss <- h1_ccs_mat[all_gene_h1,]

# png("h1_ccs_tss_heatmap.png",height = 1500,width = 900,res = 200)
h1_hm_ccs <- Heatmap(ss,use_raster = F,cluster_columns = F,cluster_rows = T,col = ccs_colors,
                     show_heatmap_legend = T,show_row_names = F,show_column_names = F,row_title = "Genes",
                     column_title = "Endoderm CCS labels",name = "CCS Label",
                     row_order = all_gene_h1)
# dev.off()



h1_chrom <- read_tsv("endo_tss_overlap.chrom_500bp.bed",col_names = F,na = character())

h1_chrom_mat <- as.matrix(h1_chrom[,-1])
rownames(h1_chrom_mat) <- h1_chrom$X1

ss <- h1_chrom_mat[all_gene_h1,]

# png("h1_chrom_tss_heatmap.png",height = 1500,width = 900,res = 200)
h1_hm_chrom <- Heatmap(ss,use_raster = F,cluster_columns = F,cluster_rows = T,col = chrom_colors,
                       show_heatmap_legend = T,show_row_names = F,show_column_names = F,row_title = "Genes",
                       column_title = "Endoderm ChromHMM labels",name ="ChromHMM Label",
                       row_order = all_gene_h1)
# dev.off()

png("endo_tss_heatmap.png",height = 1800,width = 1600,res = 200)
draw(h1_hm_ccs + h1_hm_chrom)
dev.off()

